import itertools
import logging
from pathlib import Path
import string
from typing import Dict, Tuple, List

import numpy as np

from ...data_loaders import FovScoreLogCsv
from ...figure import PubFigure, TextOutput
from ...global_config import BRAND_COLORS


logger = logging.getLogger(__name__)


DESCRIPTION = ""
FIG_WIDTH_IN = 2.6
FIG_HEIGHT_IN = 2.4

DATA_DIR = "input_data"
FOV_LOG_FILENAMES = ["fov-score-log.csv", "fov-classification-log.csv"]
N_TOP_SCORES_PER_WELL = 4
SCORE_THRESHOLDS = [-1.1, 0., 0.33, 0.66]
ACCS_SEEDED_RUNS = ["PML0439", "PML0440", "PML0441", "PML0442"]

GROUP_LABEL_ACCS = "Plates seeded\nusing automation"
GROUP_LABEL_MANUAL = "Plates seeded\nby hand"
BAR_COLORS = [
    BRAND_COLORS['net_gray_700'],
    BRAND_COLORS['orangish'],
    BRAND_COLORS['lavendish'],
    BRAND_COLORS['sf_cyan_bright'],
    "C1", "C2", "C3"
    ]


def load_all_dsets() -> Dict[str, FovScoreLogCsv]:
    dsets = {}
    data_dir = Path(__file__).parent.joinpath(DATA_DIR)
    dset_dirs = [x for x in data_dir.iterdir() if x.is_dir()]
    for dset_dir in dset_dirs:
        for fname in FOV_LOG_FILENAMES:
            file_path = dset_dir.joinpath(fname)
            if file_path.is_file():
                dset = FovScoreLogCsv.from_path(
                    file_path, missing_score_val=-1)
                n_scores = sum(
                    len(scores) for scores in dset.scores_by_well.values())
                logger.info(
                    f"Loaded {n_scores} FOV scores from {dset_dir.name!r}")
                dsets[dset_dir.name] = dset
                break
        else:
            raise Exception(
                f"No FOV score log file found in dataset dir {dset_dir.name!r}"
                )
    return dsets


def mkoutputs_fov_scores(fig_id: str
                         ) -> Tuple[List[PubFigure], List[TextOutput]]:
    dsets = load_all_dsets()
    if len(set(
            len(well_scores)
            for dset in dsets.values()
            for well_scores in dset.scores_by_well.values()
            )) > 1:
        raise Exception(
            f"Datasets don't have consistent number of sites per well")
    topnx96_scores = {}
    for dset_name, dset in dsets.items():
        topn_scores_per_well = {
            well_name: sorted(well_scores, reverse=True)[
                :N_TOP_SCORES_PER_WELL]
            for (well_name, well_scores) in dset.scores_by_well.items()
            }
        topnx96_scores[dset_name] = list(
            itertools.chain(*topn_scores_per_well.values()))

    fig = PubFigure(
        fig_id="fig_3_panel_b",
        width_in=FIG_WIDTH_IN,
        height_in=FIG_HEIGHT_IN,
        )
    ax = fig[0][0]
    dset_order = sorted(topnx96_scores.keys())
    bottoms = np.zeros(len(dset_order))
    for thresh_idx, thresh in enumerate(SCORE_THRESHOLDS):
        upper_lim = (
            SCORE_THRESHOLDS[thresh_idx + 1]
            if thresh_idx < len(SCORE_THRESHOLDS) - 1
            else 1.0
            )
        heights = np.array([
            sum(
                1 for score in topnx96_scores[dset_name]
                if ((score > thresh) if thresh_idx >= 2 else (score >= thresh))
                and (
                    (score <= upper_lim) if thresh_idx > 0
                    else (score < upper_lim)
                    )
                )
            for dset_name in dset_order
            ])
        if thresh <= -1.:
            bar_label = f"<{upper_lim:.2f}"
        else:
            bar_label = (
                ("(" if thresh_idx >= 2 else "[")
                + f"{thresh:0.2f}, {upper_lim:0.2f}]"
                )
        ax.bar(
            x=range(len(dset_order)),
            height=heights,
            bottom=bottoms,
            label=bar_label,
            fc=BAR_COLORS[thresh_idx],
            )
        bottoms += heights
    ax.set_xticks(range(len(dset_order)))
    accs_tick_counter = itertools.count(1)
    man_tick_counter = itertools.count(1)
    ax.set_xticklabels(
        (
            (
                f"A{next(accs_tick_counter)}"
                if dset_name in ACCS_SEEDED_RUNS
                else f"M{next(man_tick_counter)}"
                )
            for dset_name in dset_order
            ),
        rotation='vertical',
        rotation_mode='anchor',
        ha='right',
        va='center'
        )
    n_accs_ticks = next(accs_tick_counter) - 1
    n_man_ticks = next(man_tick_counter) - 1
    ax.set_xlim(-1., len(dset_order))
    ax.legend(
        title="FOV score range",
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncols=2,
        columnspacing=0.8,
        #reverse=True
        )
    ax.set_yticks([96 * N_TOP_SCORES_PER_WELL])
    ax.set_ylim(0, 1.05 * ax.get_yticks()[-1])
    ax.set_ylabel("# of sites in score range")
    ax_r = ax.twinx()
    ax_r.set_ylabel("% of sites in score range")
    ax_r.set_yticks([0., 25., 50., 75., 100.])
    ax_r.set_ylim(0, 1.05 * ax_r.get_yticks()[-1])

    xgroup_setup = [
        (0, n_man_ticks, GROUP_LABEL_MANUAL),
        (n_man_ticks, n_accs_ticks, GROUP_LABEL_ACCS),
        ]
    for group_idx, (start_idx, num_xvals, label) in enumerate(xgroup_setup):
        x_left = start_idx - 0.25
        x_right = start_idx + num_xvals - 0.75
        x_mid = 0.5 * (x_left + x_right)
        ax.annotate(
            "",
            xy=(x_left, -25.),
            xytext=(x_right, -25.),
            xycoords=('data', 'axes points'),
            textcoords=('data', 'axes points'),
            arrowprops={
                'arrowstyle': "-",
                'lw': 0.75,
                },
            )
        x = x_left if group_idx else x_right
        ax.annotate(
            label,
            xy=(x, -30.),
            xytext=(x, -30.),
            ha='left' if group_idx else 'right',
            va='top',
            xycoords=('data', 'axes points'),
            textcoords=('data', 'axes points'),
            arrowprops={
                'color': 'none',
                },
            )

    return [fig], []
