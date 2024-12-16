import itertools
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, List, Set, Tuple

import matplotlib as mpl
import matplotlib.axes
import matplotlib.ticker
import numpy as np

from ...figure import PubFigure, TextOutput
from ...data_loaders import NpfRunLog
from ...global_config import TITLE_FONT_SIZE, LEGEND_FONT_SIZE, BRAND_COLORS

from . import softmax_xml


logger = logging.getLogger(__name__)


DESCRIPTION = "ACCS-vs-manual Normalization performance " \
              "comparison by ATP assay"
FIG_WIDTH_IN = 5.25
FIG_HEIGHT_IN = 2.4
N_ROWS = 2
N_COLS = 3

DATA_DIR = "input_data"
MAX_SAFE_ALIQUOT_VOL = 175.
START_TIME = 59.9 * 60.
SUBFIG_DSET_NAMES = {
    'a': 'AC31121MB sv hand norm nom12k5',
    'b': 'AC31205MC hand norm 12k5 camille',
    'c': 'AC31205MR hand norm 12k5 rodrigo',
    'd': 'AC31116AR cci norm fcp 12k5',
    'e': 'AC31121AR cci norm nom12k5',
    'f': 'AC31205AR cci norm 12k5',
    }
SUBFIG_FACECOLORS = dict(
    [(k, BRAND_COLORS['sf_cyan_bright']) for k in ["d", "e", "f"]]
    #+ [(k, "#969da0") for k in ["A", "B", "C"]]
    + [(k, BRAND_COLORS['sf_cyan_bright']) for k in ["a", "b", "c"]]
    )
EXCLUDE_WELLS = {
    'AC31116AR cci norm fcp 12k5': [
        # 'G2',  # TODO: Check why this was
        ],
    'AC31128AR cci norm nom12k5': [
        'A1',  # Tip crash on reservoir
        ],
    }
ROTATED_PLATES = [
    'AC31121AR cci norm nom12k5',
    ]
GROUP_DEFS = {
    'neg': ['C3', 'F3', 'C7', 'F7'],
    'ref20': ['A5', 'C5', 'E5', 'G5', 'B10', 'D10', 'F10', 'H10'],
    'ref10': ['B5', 'D5', 'F5', 'H5', 'A10', 'C10', 'E10', 'G10'],
    'no_cells': ['F3', 'C7'],
    'no_reagent': ['C3', 'F7'],
    }
HIST_N_BINS = 32


Dset = Dict[str, List[softmax_xml.LumDataPt]]


def load_all_dsets() -> Dict[str, Dset]:
    dsets = {}
    overrange_aliquot_wells: Dict[str, Set[str]] = {}
    data_dir_path = Path(__file__).parent.joinpath(DATA_DIR)
    for dir_path in sorted(data_dir_path.iterdir()):
        runlog_paths = list(dir_path.glob("*.run_log.txt"))
        if len(runlog_paths) == 1:
            run_log = NpfRunLog.from_path(runlog_paths[0])
            aliquot_vols = run_log.aliquot_vols_actual['stock_plate']
        elif not runlog_paths:
            logger.warning(f"No run log in input dir {dir_path.name!r}")
            aliquot_vols = None
        else:
            raise Exception(f"Expected 0 or 1 *.run_log.txt files, "
                            f"got {len(runlog_paths)} in {dir_path!r}")
        for xml_path in sorted(dir_path.glob("*.xml")):
            with xml_path.open("rb") as f:
                data_by_section = softmax_xml.load_softmax_xml(f)
                for section_name in data_by_section:
                    if section_name in dsets:
                        raise Exception(
                            f"Duplicate plate name: {section_name}")
                dsets.update(data_by_section)
                logger.info(
                    f"Loaded from {xml_path.name!r}: "
                    + ", ".join(f"{x!r}" for x in data_by_section.keys()))
            for expt_name in data_by_section:
                if aliquot_vols is None:
                    overrange_aliquot_wells[expt_name] = set()
                    continue
                overrange_aliquot_wells[expt_name] = {
                    well_name
                    for (well_name, aliquot_vol) in aliquot_vols.items()
                    if aliquot_vol > MAX_SAFE_ALIQUOT_VOL
                    }
                if len(overrange_aliquot_wells[expt_name]):
                    logger.warning(f"Weak SP wells for {expt_name!r}: "
                                   f"{overrange_aliquot_wells[expt_name]}")
    return dsets


def rotate_well_name(well_name: str) -> str:
    col = int(well_name[1:])
    row = well_name[0]
    opp_row = "HGFEDCBA"["ABCDEFGH".index(row)]
    return f"{opp_row}{13 - col}"


def all_timestamp_values(dset: Dset) -> List[float]:
    return sorted({dp.time for dp in itertools.chain(*dset.values())})


def get_lum_data_at_timestamp(dset: Dset, timestamp: float) -> Dict[str, int]:
    return {
        well_name: [dp.lum for dp in dpoints if dp.time == timestamp][0]
        for (well_name, dpoints) in dset.items()
        }


def group_lum_values(
        lum_data: Dict[str, int],
        rotate_group_map: bool = False,
        exclude_wells: Iterable[str] = (),
        ) -> Dict[str, List[int]]:
    group_defs = {
        k: [rotate_well_name(w) if rotate_group_map else w for w in v]
        for (k, v) in GROUP_DEFS.items()
        }
    group_defs['test'] = [
        w for w in [f'{r}{c}' for r in "ABCDEFGH" for c in range(1, 13)]
        if w not in itertools.chain(*group_defs.values())
        ]
    return {
        group_name: [
            lum for (well_name, lum) in lum_data.items()
            if (well_name in group_wells) and (well_name not in exclude_wells)
            ]
        for (group_name, group_wells) in group_defs.items()
        }


def plot_ctg_dset(ax: mpl.axes.Axes, dset: Dset, subfig_name: str,
                  rotate_group_map: bool = False):
    dset_name = next(
        dp.ps_name for (well_name, dpoints) in dset.items() for dp in dpoints)
    all_timestamps = all_timestamp_values(dset)
    try:
        timestamp = [x for x in all_timestamps if x >= START_TIME][0]
    except IndexError:
        raise Exception(f"No data found at or after time {START_TIME:f}")
    lum_data = get_lum_data_at_timestamp(
        dset=dset, timestamp=timestamp)
    lum_data_by_group = group_lum_values(
        lum_data,
        rotate_group_map=rotate_group_map,
        exclude_wells=EXCLUDE_WELLS.get(dset_name, [])
        )
    hist_data = np.array(lum_data_by_group['test'], dtype='float64') / 1e6
    lum_mean = np.average(hist_data)
    lum_std = np.std(hist_data)
    lum_cv = lum_std / lum_mean
    logger.info(f"Stats for {dset_name!r} (subfigure {subfig_name}): "
                f"n={len(hist_data)}, mean={lum_mean:.2e}, "
                f"CV={100. * lum_cv:.1f}%")
    hist_range = (0., 2.*lum_mean)
    ax.text(
        lum_mean + lum_std,
        0.43,
        f"<keepsize>  CV={100. * lum_cv:.1f}%\n  (n={len(hist_data)})",
        transform=ax.get_xaxis_transform(),
        fontsize=LEGEND_FONT_SIZE,
        )
    ax.hist(
        hist_data,
        range=(0, 2. * lum_mean),
        # alpha=0.75,
        bins=HIST_N_BINS,
        fc=SUBFIG_FACECOLORS[subfig_name]
        # ec="black",
        )
    # ax.axvline(lum_mean, color='red', lw=1., label="mean")
    ax.axvline(lum_mean - lum_std, color='black', lw=1., ls="dotted")
    ax.axvline(
        lum_mean + lum_std,
        color='black', lw=1., ls="dotted",
        label="mean±σ"
        )
    ax.set_xlim(*hist_range)
    return {
        'dset_name': dset_name,
        'n_samples': len(hist_data),
        'mean': lum_mean,
        'cv': lum_cv,
        }


def mkoutputs_cci_vs_hand_norm(fig_id: str
        ) -> Tuple[List[PubFigure], List[TextOutput]]:
    dsets = load_all_dsets()
    subfig_dsets = {
        subfig_name: dsets[dset_name]
        for (subfig_name, dset_name) in SUBFIG_DSET_NAMES.items()
        }

    fig = PubFigure(
        fig_id="fig_2_panel_d",
        width_in=FIG_WIDTH_IN,
        height_in=FIG_HEIGHT_IN,
        n_rows=N_ROWS,
        n_cols=N_COLS,
        share_y=True,
        )
    text_lines = ["subfig, dset, n, CV"]
    for ax_idx, (subfig_name, dset) in enumerate(subfig_dsets.items()):
        row_idx = ax_idx // N_COLS
        col_idx = ax_idx % N_COLS
        ax = fig[row_idx][col_idx]
        info = plot_ctg_dset(
            ax=ax,
            dset=dset,
            subfig_name=subfig_name,
            rotate_group_map=SUBFIG_DSET_NAMES[subfig_name] in ROTATED_PLATES
            )
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1.))
        if col_idx == 0:
            ax.set_ylabel("# of samples")
            row_title = {
                0: "Manual method",
                1: "Using automation",
                }.get(row_idx, None)
            ax.text(
                #s=f"<keepsize>{row_title}",
                s=row_title,
                x=0., y=1.07,
                transform=ax.transAxes,
                #fontsize=TITLE_FONT_SIZE
                )
        if row_idx == N_ROWS - 1:
            ax.set_xlabel("Luminance (RLU)")
        if col_idx == N_COLS - 1 and row_idx == N_ROWS - 1:
            ax.legend(
                loc="upper right",
                ncols=1,
                handlelength=1.,
                handletextpad=0.15,
                frameon=True,
                borderpad=0.15,
                borderaxespad=0.15,
                # bbox_to_anchor=(0.5, -0.2)
                )
        ax.set_title(
            f"{subfig_name}", verticalalignment="top", x=0.075, y=0.8
            ).set_fontweight("bold")
        text_lines.append(f"{subfig_name}, {info['dset_name']:32s}, "
                          f"{info['n_samples']:2d}, {100. * info['cv']:.1f}%")

    text_output = TextOutput(fig_id=f"{fig_id}_stats", text_lines=text_lines)

    return [fig], [text_output]
