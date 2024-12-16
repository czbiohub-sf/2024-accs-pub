import ast
import csv
from dataclasses import dataclass
import datetime
import fnmatch
import json
import os
from pathlib import Path
import re
import string
from typing import Dict, List, Optional, TextIO, Union


class DataLoaderError(Exception):
    pass


class TextLoader:
    filename: str

    @classmethod
    def from_path(cls, path: Union[Path, str], **kwargs):
        path = Path(path)
        with path.open("r") as f:
            return cls(f=f, filename=path.name, **kwargs)

    def __init__(self, f: TextIO, filename: Optional[str] = None):
        if hasattr(f, 'name'):
            filename = os.path.basename(f.name)
        if filename is None:
            filename = "<string>"
        self.filename = filename


class JsonLoader(TextLoader):
    def __init__(self, f: TextIO, filename: Optional[str] = None):
        super().__init__(f, filename)
        self.data = json.load(f)


class CciCsv(TextLoader):
    def __init__(self, f: TextIO, filename: Optional[str] = None):
        super().__init__(f, filename)
        reader = csv.reader(f)
        self.data: Dict[str, float] = {}
        plate_col_no = 1
        for csv_row in reader:
            if not csv_row:
                continue
            for csv_col_idx, col_val in enumerate(csv_row):
                well_name = \
                    f"{string.ascii_uppercase[csv_col_idx]}{plate_col_no}"
                self.data[well_name] = float(col_val)
            plate_col_no += 1
        if fnmatch.fnmatch(self.filename, "????????-??????-*.csv"):
            ts_str = " ".join(self.filename.split("-", 2)[:2])
            self.start_ts: Optional[datetime.datetime] = (
                datetime.datetime.strptime(ts_str, "%Y%m%d %H%M%S"))
        else:
            self.start_ts = None


class NpfProtocolScript(TextLoader):
    def __init__(self, f: TextIO, filename: Optional[str] = None):
        super().__init__(f, filename)
        self.body = ast.parse(f.read()).body

    def eval_toplevel_assign_value(self, target_name: str):
        for stmt in self.body:
            if not isinstance(stmt, ast.Assign):
                continue
            for target in stmt.targets:
                if not isinstance(target, ast.Name):
                    continue
                if target.id != target_name:
                    continue
                code = compile(
                    ast.Expression(stmt.value),
                    filename=self.filename,
                    mode="eval")
                return eval(code)
        else:
            raise DataLoaderError(f"No assign to name {target_name} found "
                                  f"in top level of {self.filename!r}")

    def get_run_config(self):
        return self.eval_toplevel_assign_value("run_config")


class NpfRunLog(TextLoader):
    def __init__(self, f: TextIO, filename: Optional[str] = None):
        super().__init__(f, filename)
        self.run_id = None
        self.robot_name = None
        self.start_time = None
        self.finish_time = None
        self.aliquot_vols_intended: Dict[str, Dict[str, float]] = {}
        self.aliquot_vols_actual: Dict[str, Dict[str, float]] = {}
        self.cci_readings: Dict[str, float] = {}
        self.cci_valid_areas: Dict[str, float] = {}

        line_start_pat = re.compile(
            r"^\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2}\.\d+\s+[A-Z]+:\s+")
        aliq_tgt_pat = re.compile(
            r"Plate '(?P<plate_name>.*)' well (?P<well_name>\w+\d+): "
            r"target = (?P<tgt_cell_count>\d+\.?\d*)k "
            r"(?:\(x (?P<tgt_multiplier>\d+\.?\d*)\) )?cells, "
            r"aliquot = (?P<aliquot_vol>\d+\.?\d*) μL")
        oor_pat = re.compile(
            r"Plate '(?P<plate_name>.*)' well (?P<well_name>\w+\d+): "
            r"seeding target \((?P<tgt_cell_count>\d+\.?\d*)k cells -> "
            r"(?P<tgt_aliquot_vol>\d+\.?\d*) μL\) out of pipetting range - "
            r"aliquoting (?P<max_or_min>max|min) volume "
            r"\((?P<clip_aliquot_vol>\d+\.?\d*) μL\)")
        cci_rdg_pat = re.compile(
            r"Well (?P<well_name>\w+\d+): "
            r"density = (?P<cci_reading>\d+\.?\d*) cells/μL, "
            r"valid_area = (?P<cci_valid_area>\d+\.?\d*)%")

        self.tgt_multiplier = None
        for raw_line in (" ".join(x.strip().split()) for x in f.readlines()):
            if not line_start_pat.match(raw_line):
                continue
            left, right = (x.strip() for x in raw_line.split(": ", 1))
            datetime_str, msg_type = left.rsplit(None, 1)
            timestamp = (
                datetime.datetime.strptime(datetime_str,
                                           "%Y-%m-%d %H:%M:%S.%f")
                .replace(tzinfo=datetime.timezone.utc)
                )
            if msg_type in ["INFO", "NOTICE"]:
                if right.startswith("Protocol starting:"):
                    if "(" in right:
                        self.run_id = right.split(
                            "(", 1)[-1].strip().rstrip(")")
                    self.start_time = timestamp
                    continue
                elif right.startswith("Protocol complete."):
                    self.finish_time = timestamp
                    continue
            if msg_type == "WARNING":
                match = oor_pat.match(right)
                if match is not None:
                    plate_name = match.group('plate_name')
                    well_name = match.group('well_name')
                    tgt_vol = float(match.group('tgt_aliquot_vol'))
                    if plate_name not in self.aliquot_vols_intended:
                        self.aliquot_vols_intended[plate_name] = {}
                    self.aliquot_vols_intended[plate_name][well_name] = tgt_vol
                    continue
            elif msg_type == "INFO":
                if right.startswith("Robot info: name:"):
                    self.robot_name = right.rsplit(":", 1)[-1].strip()
                    continue
                match = cci_rdg_pat.match(right)
                if match is not None:
                    well_name = match.group('well_name')
                    cci_reading = float(match.group('cci_reading'))
                    cci_valid_area = (
                        float(match.group('cci_valid_area')) / 100.)
                    self.cci_readings[well_name] = cci_reading
                    self.cci_valid_areas[well_name] = cci_valid_area
                    continue
                match = aliq_tgt_pat.match(right)
                if match is not None:
                    plate_name = match.group('plate_name')
                    well_name = match.group('well_name')
                    tgt_count = float(match.group('tgt_cell_count'))
                    aliquot_vol = float(match.group('aliquot_vol'))
                    if match.group('tgt_multiplier') is not None:
                        self.tgt_multiplier = float(
                            match.group('tgt_multiplier'))
                    if plate_name not in self.aliquot_vols_actual:
                        self.aliquot_vols_actual[plate_name] = {}
                    self.aliquot_vols_actual[
                        plate_name][well_name] = aliquot_vol
                    continue

        for plate_name in self.aliquot_vols_actual:
            if plate_name not in self.aliquot_vols_intended:
                self.aliquot_vols_intended[plate_name] = {}
            for well_name, aliquot_vol in self.aliquot_vols_actual[
                    plate_name].items():
                if well_name not in self.aliquot_vols_intended[plate_name]:
                    self.aliquot_vols_intended[plate_name][
                        well_name] = aliquot_vol


@dataclass
class FovScoreLogRecord:
    well_name: str
    site_id: int
    timestamp: datetime.datetime
    score: Optional[float]


class FovScoreLogCsv(TextLoader):
    MISSING_SCORE_VAL = None

    def __init__(self, f: TextIO, filename: Optional[str] = None,
                 missing_score_val: Optional[float] = None):
        if missing_score_val is None:
            missing_score_val = self.MISSING_SCORE_VAL
        super().__init__(f, filename)
        reader = csv.DictReader(f)
        scores_by_well: Dict[str, List[Optional[float]]] = {}
        self.records: Dict[str, Dict[int, FovScoreLogRecord]] = {}
        self.missing_well_names = False
        first_fov_ts = None
        for row in reader:
            ts = datetime.datetime.strptime(
                row['timestamp'], "%Y-%m-%d %H:%M:%S")
            if first_fov_ts is None or ts < first_fov_ts:
                first_fov_ts = ts
            if 'position_well_id' not in row:
                self.missing_well_names = True
                well_name = f"unk-{int(row['position_ind']):08d}"
            else:
                well_name = row['position_well_id']
            if row['score']:
                score: Optional[float] = float(row['score'])
            else:
                score = missing_score_val
            site_id = int(row['position_site_num'])
            record = FovScoreLogRecord(
                well_name=well_name,
                site_id=site_id,
                timestamp=ts,
                score=score
                )
            if well_name not in self.records:
                self.records[well_name] = {}
            self.records[well_name][site_id] = record
            if well_name not in scores_by_well:
                scores_by_well[well_name] = []
            scores_by_well[well_name].append(score)
        self.first_fov_ts = first_fov_ts
        self.scores_by_well = scores_by_well


class DfaExperimentMetadataJson(JsonLoader):
    pass
