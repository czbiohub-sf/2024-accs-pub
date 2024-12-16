import csv
from dataclasses import dataclass
import logging
from pathlib import Path
import string
from typing import Dict, Iterable, List, Mapping, Optional, TextIO, Tuple

import numpy as np

from ...analysis_common import PathOrStr
from ...data_loaders import TextLoader
from ...figure import PubFigure, TextOutput
from ...global_config import BRAND_COLORS, XYLABEL_FONT_SIZE


logger = logging.getLogger(__name__)


DESCRIPTION = "Hands-on time profiling for manual method vs ACCS runs"
FIG_WIDTH_IN = 4.3
FIG_HEIGHT_IN = 2.4

DATA_DIR = "input_data"
CAT_MAP = {
    'setup': ["setup_cci", "setup_labware", "setup_reagents", "setup_sw",
              "setup_trashbin", "annotation", "setup_bsc", "setup_prerun"],
    # 'annotation': ["annotation"],
    'arraying': ["arraying"],
    'passaging': ["wash", "trypsinize", "quench", "seeding", "inspection"],
    'walk_away': ["walk_away"],
    'nocount': ["pause"],
    'cleanup': ["cleanup_cci", "cleanup_bsc", "cleanup_plates",
                "cleanup_plate", "cleanup_labware", "cleanup_trashbin"],
    }
PLOT_DSET_ORDER_MANUAL = [
    "2023-11-14 manual SV",
    "2023-11-16 manual RB",
    "2023-11-28 manual CJ",
    ]
PLOT_DSET_ORDER_ACCS = [
    "2023-11-22 ACCS GC",
    "2023-11-30 ACCS SV",
    ]
PLOT_CAT_ORDER = [
    "setup", "arraying", "passaging", "walk_away", "cleanup"]
PLOT_CAT_DESCS = {
    'setup': "Setup / annotation",
    # 'annotation': "Annotation",
    'arraying': "Pipette tip arraying",
    'walk_away': "Walk-away",
    'passaging': "Manual passaging",
    'cleanup': "Teardown / cleanup",
    }
GROUP_LABEL_ACCS = "Using\nautomation"
GROUP_LABEL_MANUAL = "Manual\nmethod"
PLOT_CAT_STYLES = {
    'setup': {'facecolor': BRAND_COLORS['orangish']},
    'arraying': {'facecolor': BRAND_COLORS['rustish']},
    'walk_away': {'facecolor': "#FFFFFF00"},
    'passaging': {'facecolor': BRAND_COLORS['lavendish']},
    'cleanup': {'facecolor': BRAND_COLORS['sf_cyan_bright']},
    }
PLOT_BAR_HEIGHT = 0.33


@dataclass
class TimeProfSegment:
    name: str
    start_time: float
    end_time: float

    def duration(self):
        if self.start_time > self.end_time:
            raise ValueError("Negative duration?")
        return self.end_time - self.start_time


class TimeProfDataset(TextLoader):
    def __init__(self, f: TextIO, filename: Optional[str] = None):
        super().__init__(f, filename)
        self._fname_fmtd = f"{filename!r}" if filename is not None else "??"
        reader = csv.DictReader(f)
        self.segments = []
        open_segments: Dict[str, float] = {}
        for row_idx, csv_row in enumerate(reader):
            csv_row = {k: v.strip() for (k, v) in csv_row.items()}
            if not csv_row.get('timestamp', None):
                continue
            timestamp = float(csv_row['timestamp']) / 1e3
            start_tag = csv_row['start']
            stop_tag = csv_row['stop']
            if stop_tag:
                if stop_tag not in open_segments:
                    raise Exception(
                        f"{self._fname_fmtd} row {row_idx+1}: segment "
                        f"{stop_tag!r} not open"
                        )
                self.segments.append(TimeProfSegment(
                    name=stop_tag,
                    start_time=open_segments[stop_tag],
                    end_time=timestamp
                    ))
                del open_segments[stop_tag]
            if start_tag:
                if start_tag in open_segments:
                    raise Exception(
                        f"{self._fname_fmtd} row {row_idx+1}: segment "
                        f"{start_tag!r} already open"
                        )
                open_segments[start_tag] = timestamp
        if open_segments:
            raise Exception(
                "Segment(s) not closed: "
                + ", ".join(f"{name!r}" for name in open_segments)
                )

    def get_durations(self):
        durations = {}
        for segment in self.segments:
            if segment.name not in durations:
                durations[segment.name] = 0.
            durations[segment.name] += segment.duration()
        return durations

    def get_endtoend_duration(self) -> float:
        if not self.segments:
            return 0.
        return (
            max(s.end_time for s in self.segments)
            - min(s.start_time for s in self.segments)
            )

    def _no_x_occurrences(self, tag: str):
        raise Exception(f"No {tag!r} occurrences in {self._fname_fmtd}")

    def get_duration(self, tag: str, or_zero: bool = False):
        sub_durs = [s.duration() for s in self.segments if s.name == tag]
        if not sub_durs:
            if or_zero:
                return 0.
            self._no_x_occurrences(tag)
        return sum(sub_durs)

    def merge_durations(self, cat_map: Mapping[str, Iterable[str]],
                        nonexistent_ok: bool = True,
                        include_zeros: bool = False):
        durations = self.get_durations()
        result = {}
        included = set()
        for cat_name, cat_tags in cat_map.items():
            cat_tags = set(cat_tags)
            result[cat_name] = sum(
                self.get_duration(tag, or_zero=nonexistent_ok)
                for tag in cat_tags
                )
            included.update(cat_tags)
        left_out = set(durations.keys()).difference(included)
        if left_out:
            logger.warning(
                f"Tags left out of totals from {self._fname_fmtd}: "
                + ", ".join(f"{tag_name!r}" for tag_name in left_out)
                )
        return {
            k: v for (k, v) in result.items()
            if v > 0. or include_zeros
            }

    def get_time_between(self, tag1: str, tag2: str):
        tag1_evts = sorted(
            (s for s in self.segments if s.name == tag1),
            key=lambda s: s.start_time
            )
        if not tag1_evts:
            self._no_x_occurrences(tag1)
        tag2_evts = sorted(
            (s for s in self.segments if s.name == tag2),
            key=lambda s: s.end_time,
            )
        if not tag2_evts:
            self._no_x_occurrences(tag2)
        return tag2_evts[-1].end_time - tag1_evts[0].start_time


def load_all_dsets(data_dir: Optional[PathOrStr] = None
                   ) -> Dict[str, TimeProfDataset]:
    if data_dir is None:
        data_dir = Path(__file__).parent.joinpath(DATA_DIR)
    else:
        data_dir = Path(data_dir)
    dsets = {}
    for csv_path in sorted(data_dir.glob("*.events.csv")):
        dset_name = csv_path.name.rsplit(".", 2)[0]
        dset = TimeProfDataset.from_path(csv_path)
        logger.info(f"Dataset {dset_name!r}: {len(dset.segments)} events, "
                    f"{len(dset.get_durations())} tags")
        dsets[dset_name] = dset
    return dsets


def mkoutputs_time_prof(fig_id: str
                        ) -> Tuple[List[PubFigure], List[TextOutput]]:
    dsets = load_all_dsets()
    logger.info(f"Loaded {len(dsets)} dataset(s)")

    text_lines = []
    cat_durations = {}
    for dset_idx, (dset_name, dset) in enumerate(dsets.items()):
        text_lines.append(f"Dataset: {dset_name!r}")
        endtoend_s = dset.get_endtoend_duration()
        handson_s = endtoend_s - dset.get_duration('hands_off', or_zero=True)
        durations = dset.get_durations()
        uaf_s = endtoend_s - sum(durations.values())
        if uaf_s > 0.:
            logger.warning(f"Dataset {dset_name!r} has {uaf_s:.1f} s "
                           "of unaccounted-for time")
        text_lines.append(
            f"  End-to-end time:        {endtoend_s / 60.:6.1f} min")
        text_lines.append(
            f"  Hands-on time:          {handson_s / 60.:6.1f} min")
        text_lines.append(
            f"  Unaccounted-for time:   {uaf_s / 60.:6.1f} min")

        def print_durs(durations):
            dur_order = sorted(
                durations.items(), key=lambda x: x[1], reverse=True)
            for name, dur_s in dur_order:
                text_lines.append(f"    {name!r:20s}: {dur_s / 60.:6.1f} min")
        text_lines.append("  Durations by tag:")
        print_durs(durations)
        cat_durations[dset_name] = dset.merge_durations(CAT_MAP)
        text_lines.append("  Durations by category:")
        print_durs(cat_durations[dset_name])
        if dset_idx != len(dsets) - 1:
            text_lines.append("")

    plot_dset_order = PLOT_DSET_ORDER_MANUAL + PLOT_DSET_ORDER_ACCS
    plot_total_durs = {
        dset_name: sum(
            dset_cat_durs.get(cat_name, 0.) for cat_name in PLOT_CAT_ORDER)
        for (dset_name, dset_cat_durs) in cat_durations.items()}
    fig = PubFigure(
        fig_id="fig_3_panel_a",
        width_in=FIG_WIDTH_IN,
        height_in=FIG_HEIGHT_IN,
        )
    ax = fig[0][0]
    lefts = np.zeros(len(plot_dset_order))
    y_vals = list(reversed(range(len(plot_dset_order))))
    walkaway_bars = {}
    walkaway_durs_min = {}
    for cat_name in PLOT_CAT_ORDER:
        widths = np.array([
            cat_durations[dset_name].get(cat_name, 0.) / 60.
            for dset_name in plot_dset_order
            ])
        bars = ax.barh(
            y=y_vals,
            width=widths,
            left=lefts,
            label=(
                PLOT_CAT_DESCS[cat_name] if cat_name != 'walk_away'
                else None
                ),
            height=PLOT_BAR_HEIGHT,
            **PLOT_CAT_STYLES[cat_name],
            )
        if cat_name == 'walk_away':
            walkaway_dset_idxs = [
                i for (i, dset_name) in enumerate(plot_dset_order)
                if 'walk_away' in cat_durations[dset_name]
                ]
            walkaway_bars = {
                plot_dset_order[i]: bars[i] for i in walkaway_dset_idxs
                }
            walkaway_durs_min = {
                plot_dset_order[i]: widths[i] for i in walkaway_dset_idxs
                }
        lefts += widths
    if walkaway_bars is not None:
        assert walkaway_durs_min is not None
        for dset_name, bar in walkaway_bars.items():
            for is_left in True, False:
                annotate_x = bar.get_x() if is_left else bar.get_bbox().max[0]
                ax.annotate(
                    f"<keepsize>({walkaway_durs_min[dset_name]:.1f}' walk-away time)",
                    xy=(annotate_x, bar.get_center()[1]),
                    xytext=bar.get_center(),
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=XYLABEL_FONT_SIZE,
                    color="#FFFFFF00" if is_left else None,
                    arrowprops={
                        'arrowstyle': "->",
                        'lw': 1.0,
                        'color': "#969da0",
                        },
                    )
    for dset_idx, y_pos in enumerate(y_vals):
        dset_name = plot_dset_order[dset_idx]
        total_dur_min = plot_total_durs[dset_name] / 60.
        handson_min = (
            total_dur_min - cat_durations[dset_name].get('walk_away', 0.) / 60.
            )
        total_desc = "<keepsize>" + (
            f"{total_dur_min:.1f}' total, "
            if 'walk_away' in cat_durations[dset_name]
            else ""
            )
        ax.annotate(
            total_desc + f"{handson_min:.1f}' hands-on",
            xy=(total_dur_min, y_pos + PLOT_BAR_HEIGHT / 2.),
            xytext=(-10., 10.),
            textcoords="offset points",
            fontsize=XYLABEL_FONT_SIZE,
            ha='right',
            va='center_baseline',
            arrowprops={
                'arrowstyle': "->",
                'lw': 0.75,
                'connectionstyle': "angle,angleA=0,angleB=90",
                },
            )
    ax.set_ylim(-0.5, len(plot_dset_order))
    ax.set_xlabel("Time (minutes)")
    ax.set_yticks(y_vals)
    ax.set_yticklabels(string.ascii_uppercase[:len(ax.get_yticks())])
    ygroup_setup = [
        (0, len(PLOT_DSET_ORDER_ACCS), GROUP_LABEL_ACCS),
        (len(PLOT_DSET_ORDER_ACCS), len(PLOT_DSET_ORDER_MANUAL),
            GROUP_LABEL_MANUAL),
        ]
    for start_idx, num_yvals, label in ygroup_setup:
        y_bottom = start_idx - 0.25
        y_top = start_idx + num_yvals - 0.75
        y_mid = 0.5 * (y_bottom + y_top)
        ax.annotate(
            "",
            xy=(-25., y_bottom),
            xytext=(-25., y_top),
            xycoords=('axes points', 'data'),
            textcoords=('axes points', 'data'),
            arrowprops={
                'arrowstyle': "-",
                'lw': 0.75,
                },
            )
        ax.annotate(
            label,
            xy=(-30., y_mid),
            xytext=(-30., y_mid),
            ha='right',
            va='center',
            xycoords=('axes points', 'data'),
            textcoords=('axes points', 'data'),
            arrowprops={
                'color': 'none',
                },
            )
    ax.legend(title="Hands-on time", alignment='left')

    text_output = TextOutput(fig_id=fig_id, text_lines=text_lines)

    return [fig], [text_output]
