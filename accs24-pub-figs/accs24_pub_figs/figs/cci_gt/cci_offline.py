import dataclasses
import hashlib
import importlib
import itertools
import json
import logging
from pathlib import Path
import sys

import numpy as np


CCI_LIB_DIR_BASE = "cci_libs"
STD_CONFIG_NAME = "cci25"
LOCAL_CONFIG_NAME = "cci25_local"


logger = logging.getLogger(__name__)


def import_cci_package(ver_name):
    return import_package_from_path(
        cci_lib_path(ver_name),
        mod_name=ver_name
        )


def cci_ver_base_path(ver_name):
    return Path(CCI_LIB_DIR_BASE).joinpath(ver_name)


def cci_lib_path(ver_name):
    return cci_ver_base_path(ver_name).joinpath("cell_counting_imager")


def cci_std_config_path(ver_name):
    return cci_lib_path(ver_name)\
        .joinpath("configs")\
        .joinpath(f"{STD_CONFIG_NAME}.json")


def cci_local_config_path(ver_name):
    return cci_ver_base_path(ver_name).joinpath(f"{LOCAL_CONFIG_NAME}.json")


def hash_digest_for_path(path):
    with Path(path).open("rb") as f:
        hash_ = hashlib.sha1(f.read())
    return hash_.hexdigest()


def import_package_from_path(path, mod_name=None):
    # yuck
    path = Path(path)
    if mod_name is None:
        mod_name = path.stem

    sys.path.insert(0, str(path.parent.resolve()))
    spec = importlib.util.find_spec(path.stem)
    sys.path.pop(0)

    spec.name = mod_name
    spec.loader.name = mod_name
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


class CciOfflineAnalysis:
    def __init__(self, lib_path, config_paths=(), cci_ver_name=None):
        lib_path = Path(lib_path)
        self.cci_ver_name = (
            cci_ver_name if cci_ver_name is not None
            else lib_path.parent.name)
        self._cci_lib = import_package_from_path(lib_path)
        self._cci_configr = self._cci_lib.CciConfigurator()
        for config_path in config_paths:
            self._cci_configr.load_config_from_file(config_path)
        ch_names = sorted({
            ch_name 
            for (ch_name, ch_pos)
            in self._cci_configr.config.hw_params['lane_positions']})
        self.cell_counters = {
            ch_name: self._cci_configr.get_counter()
            for ch_name in ch_names
            }
        self._bg_hashes = {}
        self._last_bg_hash = {ch_name: None for ch_name in ch_names}

    def analyze_image_pair(self, ch_name, fg_path, bg_path):
        fg_path = Path(fg_path).resolve()
        bg_path = Path(bg_path).resolve()
        counter = self.cell_counters[ch_name]
        if bg_path not in self._bg_hashes:
            self._bg_hashes[bg_path] = hash_digest_for_path(bg_path)
        if self._bg_hashes[bg_path] != self._last_bg_hash[ch_name]:
            bg_id = bg_path.stem.split("-", 1)[-1]
            logger.info(f"Analyzing new BG: {bg_id}")
            counter.set_bg_image_from_path(bg_path)
            self._last_bg_hash[ch_name] = self._bg_hashes[bg_path]
        cc_result = counter.process_fg_image_from_path(fg_path)
        return CciOfflineResult(
            fg_name=fg_path.stem,
            bg_name=bg_path.stem,
            cci_ver=self.cci_ver_name,
            cells_per_ml=cc_result.cells_per_ml,
            n_cells=len(cc_result.cell_locations),
            n_valid_pixels=np.count_nonzero(cc_result.feature_mask),
            n_total_pixels=cc_result.feature_mask.size
            )

    def analyze_all_runs(self, top_path, **kwargs):
        return self._analyze_dir(top_path, level_no=2, **kwargs)

    def analyze_run_dir(self, run_dir_path, **kwargs):
        return self._analyze_dir(run_dir_path, level_no=1, **kwargs)

    def _analyze_dir(self, dir_path, level_no=1, **kwargs):
        dir_path = Path(dir_path)
        if level_no:
            return list(itertools.chain(*(
                self._analyze_dir(sub_dir, level_no-1, **kwargs)
                for sub_dir in sorted(dir_path.glob("????????-??????"))
                if sub_dir.is_dir())))
        return self.analyze_scan_dir(dir_path, **kwargs)

    def analyze_scan_dir(self, dir_path, results_dir=None,
                         skip_processed=True):
        results = []
        dir_path = Path(dir_path)
        for fg_path in sorted(dir_path.glob("image-input-*-*-*.tif")):
            meas_id = fg_path.stem.split("-", 2)[-1]
            scan_id, ch_name = meas_id.rsplit("-", 1)
            bg_path = fg_path.parent.joinpath(f"background-{meas_id}.tif")
            if results_dir is not None:
                json_path = Path(results_dir).joinpath(f"{meas_id}.json")
                if json_path.exists() and skip_processed:
                    logger.info(f"Using existing result for FG: {meas_id}")
                    results.append(CciOfflineResult.read_from_file(json_path))
                    continue
            logger.info(f"Analyzing FG: {meas_id}")
            result = self.analyze_image_pair(ch_name, fg_path, bg_path)
            if results_dir is not None:
                result.write_to_file(json_path)
                logger.debug(f"Wrote result to {str(json_path)!r}")
            results.append(result)
            logger.debug(repr(result))
        return results


@dataclasses.dataclass(frozen=True)
class CciOfflineResult:
    fg_name: str
    bg_name: str
    cci_ver: str
    cells_per_ml: float
    n_cells: int
    n_valid_pixels: int
    n_total_pixels: int

    def cells_per_pixel(self):
        return n_cells / n_valid_pixels

    def meas_id(self):
        return self.fg_name.split("-", 2)[-1]

    def ch_name(self):
        return self.fg_name.rsplit("-", 1)[-1]

    def scan_id(self):
        return self.meas_id().rsplit("-", 1)[0]

    def write_to_file(self, path):
        path = Path(path)
        with path.open("w") as f:
            json.dump(dataclasses.asdict(self), f)

    @classmethod
    def read_from_file(cls, path):
        path = Path(path)
        with path.open("r") as f:
            data = json.load(f)
        return cls(**data)
