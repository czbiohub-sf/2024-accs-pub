from dataclasses import dataclass
import itertools
import logging
import math
from pathlib import Path
from typing import Collection, Dict, List, Sequence, Tuple

import matplotlib as mpl
import numpy as np
import numpy.typing

from ...analysis_common import (
    PathOrStr, build_gating_strategy, process_fcs_files, count_sample,
    check_fcs_time_data)
from ...data_loaders import NpfProtocolScript, CciCsv
from ...figure import PubFigure, TextOutput
from ...global_config import BRAND_COLORS, LEGEND_FONT_SIZE
#from . import cci_norm_fc_sim


logger = logging.getLogger(__name__)


DESCRIPTION = "CCI normalization performance measurement by FC"
SUMMARY_FIG_WIDTH_IN = 1.8
SUMMARY_FIG_HEIGHT_IN = 2.25
GROUPS_FIG_WIDTH_IN = 2.5
GROUPS_FIG_HEIGHT_IN = 2.25
#MC_FIG_WIDTH_IN = 3.25
#MC_FIG_HEIGHT_IN = 2.5
STRIPPLOT_RAND_SEED = 499

DATA_DIR = "input_data"
SUBFIG_DSET_NAMES = {
    'A': "2023-03-14 crazyplate test",
    'B': "2023-05-17 crazyplate test",
    #'C': "2023-10-31 crazyplate test",
    'C': "2023-12-14 crazyplate test",
    }
GROUPS_FIG_USE_DSET = "2023-03-14 crazyplate test"
BEAD_WELLS = {
    dset_name: (
        [f"{r}4" for r in "ABCDEFGH"]
        if dset_name != '2023-05-17 crazyplate test'
        else ([f"{r}4" for r in "ABCDEFG"] + ["H5"])
        # B4 didn't get correct amount but we're not analyzing beads here
        )
    for dset_name in SUBFIG_DSET_NAMES.values()
    }
EXCLUDE_WELLS = {
    '2023-03-14 crazyplate test': [
        #"H3", "H9", "B1", "D3", "F11",  # Bad time data in FCS
        "B1", "E1", "B2",  # Image obstructions
        ],
    '2023-05-17 crazyplate test': [
        #"F1", "F3", "G2", "G12",  # Bad time data in FCS
        ],
    '2023-10-31 crazyplate test': [
        #"D12", "E7", "F11", "G2",  # Bad time data in FCS
        ],
    '2023-12-14 crazyplate test': [
        #"D3", "D11", "F7", "H11", "A1",  # Bad time data in FCS
        "A1", "B1", "H1"  # Image obstructions
        ],
    }
FC_SAMPLE_VOL_UL = 60.
FC_SAMPLE_LOADING_TIME_S = 20.
FC_TIME_TOL_S = 1.
SEEDING_TARGET_CELLS = 12.5e3
MAX_RELIABLE_ALIQUOT_UL = 180.
HIST_N_BINS = 48
HIST_X_MAX = 300.
BAR_COLORS = [
    BRAND_COLORS['rustish'],
    BRAND_COLORS['lavendish'],
    BRAND_COLORS['sf_cyan'],
    BRAND_COLORS['orangish'],
    "C1", "C2", "C3"
    ]


@dataclass
class Dataset:
    crazyplate_conc: Dict[str, float]
    splate_cci_data: Dict[str, float]
    oplate_fc_data: Dict[str, float]
    n_bad_fcs: int


def load_dataset(dset_name: str, dset_dir: PathOrStr) -> Dataset:
    dset_dir_path = Path(dset_dir)
    script_dir = dset_dir_path.joinpath("crzp_script")
    try:
        script_path = next(script_dir.glob("*.script_contents.[tT][xX][tT]"))
    except StopIteration:
        raise Exception(f"Missing protocol script for {dset_dir_path.name!r}")
    cci_csv_dir = dset_dir_path.joinpath("inp_cci")
    try:
        cci_csv_path = next(
            cci_csv_dir.glob("????????-??????-OT2Counts.[cC][sS][vV]"))
    except StopIteration:
        raise Exception(f"Missing CCI CSV for {dset_dir_path.name!r}")
    fcs_dir = dset_dir_path.joinpath("outp_fcs")
    fcs_paths = list(fcs_dir.glob("*_*_*.[fF][cC][sS]"))
    if not fcs_paths:
        raise Exception(f"Missing FCS files for {dset_dir_path.name!r}")
    excl_well_names = set(
        BEAD_WELLS[dset_name] + EXCLUDE_WELLS.get(dset_name, []))

    proto_script = NpfProtocolScript.from_path(script_path)
    run_config = proto_script.get_run_config()
    split_factors = (
        run_config['cell_splitting']['target_split_factors']['stock_plate'])
    crazyplate_conc = {
        well_name: (1. / split_factor) if split_factor != 0. else 0.
        for (well_name, split_factor) in split_factors.items()
        }
    logger.info(f"Loaded crazyplate map for {dset_dir_path.name!r} "
                f"from {str(script_path)!r}")

    cci_data = CciCsv.from_path(cci_csv_path)
    logger.info(f"Loaded crazyplate CCI readings for {dset_dir_path.name!r} "
                f"from {str(cci_csv_path)!r}")

    gating_strategy = build_gating_strategy()
    fcs_results = process_fcs_files(
        fcs_paths,
        cb_fn=lambda s: count_sample(
            sample=s,
            gating_strategy=gating_strategy,
            skip_wells=excl_well_names,
            sample_loading_time_s=FC_SAMPLE_LOADING_TIME_S,
            time_tol_s=FC_TIME_TOL_S
            ),
        )
    bad_well_names = [
        result['well_name'] for result in fcs_results.values()
        if result['suspicious']
        ]
    logger.info(
        "Excluding wells with suspect FC data: " + ", ".join(bad_well_names)
        )
    oplate_fc_data = {
        result['well_name']: result['n_cells'] / FC_SAMPLE_VOL_UL
        for result in fcs_results.values()
        if result['analyzed'] and not result['suspicious']
        }
    logger.info(f"Processed {len(fcs_results)} FCS files yielding "
                f"{len(oplate_fc_data)} results for "
                f"{dset_dir_path.name!r} from {str(fcs_dir)!r}")
    return Dataset(
        oplate_fc_data=oplate_fc_data,
        n_bad_fcs=len(bad_well_names),
        splate_cci_data=cci_data.data,
        crazyplate_conc=crazyplate_conc
        )


def load_all_datasets():
    data_dir_path = Path(__file__).parent.joinpath(DATA_DIR)
    return {
        dset_name: load_dataset(dset_name, data_dir_path.joinpath(dset_name))
        for dset_name in SUBFIG_DSET_NAMES.values()
        }


def make_violin_plot(ax: mpl.axes.Axes,
                     dataset: Sequence[np.typing.ArrayLike],
                     *args, raincloud: bool = True, clip_hbars: bool = True,
                     scatter_alpha: float = 0.5, **kwargs):
    violinplot_colls = ax.violinplot(dataset, *args, **kwargs)
    for coll in [
            coll for (coll_name, coll) in violinplot_colls.items()
            if coll_name in ['cmeans', 'cmins', 'cmaxes', 'cbars', 'cmedians']
            ]:
        coll.set_linewidth(1.)
    violinplot_colls['cbars'].set_edgecolor(BRAND_COLORS['sf_cyan'])
    vbar_segs = violinplot_colls['cbars'].get_segments()
    vbar_xs = [x0 for ((x0, y0), (x1, y1)) in vbar_segs]
    ax.set_xlim(-0.5, len(vbar_xs) - 0.5)
    if not raincloud:
        return violinplot_colls
    for coll_name in ['cmeans', 'cmaxes', 'cmins']:
        hbar_coll = violinplot_colls[coll_name]
        hbar_coll.set_edgecolor(BRAND_COLORS['sf_cyan'])
        if clip_hbars:
            hbar_segs = hbar_coll.get_segments()
            hbar_coll.set_segments([
                ((x0, y0), (vbar_xs[violin_idx]+0.05, y1))
                for (violin_idx, ((x0, y0), (x1, y1))) in enumerate(hbar_segs)
                ])
    vbar_color = violinplot_colls['cbars'].get_colors()[0]
    for violin_idx, violin_x in enumerate(vbar_xs):
        body = violinplot_colls['bodies'][violin_idx]
        body.set_facecolor(BRAND_COLORS['sf_cyan'])
        verts = body.get_paths()[0].vertices
        new_verts = [(min(x, violin_x), y) for (x, y) in verts]
        body.set_paths([new_verts])
        scat_ys = dataset[violin_idx]
        scat_xs = np.random.default_rng(STRIPPLOT_RAND_SEED).uniform(
            violin_x + 0.1, violin_x + 0.3, len(scat_ys))
        ax.scatter(
            scat_xs, scat_ys,
            ec='none', fc=vbar_color, alpha=scatter_alpha, s=8.
            )
    return violinplot_colls


def group_data_by_crazyplate_group(
        data: Dict[str, float],
        seeding_density_groups: Dict[str, List[str]],
        excl_wells: Collection[str]
        ):
    group_names = list(seeding_density_groups.keys())
    group_data = [
        [
            data[well_name]
            for well_name in seeding_density_groups[group_name]
            if well_name in data and well_name not in excl_wells
            ]
        for group_name in group_names
        ]
    if len(set(len(x) for x in group_data)) > 1:
        # This is a dumb hack to get avoid MPL triggering a numpy deprecation warning
        group_data = [np.asarray(x, dtype="float") for x in group_data]
        group_data = np.asarray(group_data, dtype="object")
    return group_data


def get_sp_weak_wells(cci_data: Dict[str, float], critical_vol: float,
                      seeding_target_cells: float):
    aliquot_targets = {
        well_name: (
            (seeding_target_cells / cci_cell_conc)
            if cci_cell_conc > 0
            else None
            )
        for (well_name, cci_cell_conc) in cci_data.items()
        }
    return {
        well_name for (well_name, aliquot_vol) in aliquot_targets.items()
        if (aliquot_vol is None) or (aliquot_vol > critical_vol)
        }


def mkoutputs_cci_norm_fc(fig_id: str
                          ) -> Tuple[List[PubFigure], List[TextOutput]]:
    dsets = load_all_datasets()
    subfig_dsets = {
        subfig_name: dsets[dset_name]
        for (subfig_name, dset_name) in SUBFIG_DSET_NAMES.items()
        }

    overview_fig = PubFigure(
        fig_id=f"fig_2_panel_b",
        width_in=SUMMARY_FIG_WIDTH_IN,
        height_in=SUMMARY_FIG_HEIGHT_IN
        )
    ax = overview_fig[0][0]

    text_lines = ["subfig, dset, n, mean, CV"]
    plot_data = []
    plot_cats = []
    plot_infos = []
    for cat_idx, (subfig_name, dset) in enumerate(subfig_dsets.items()):
        dset_name = SUBFIG_DSET_NAMES[subfig_name]
        man_excl_well_names = EXCLUDE_WELLS.get(dset_name, [])
        excl_well_names = set(BEAD_WELLS[dset_name] + man_excl_well_names)
        sp_weak_wells = get_sp_weak_wells(
            dset.splate_cci_data,
            critical_vol=MAX_RELIABLE_ALIQUOT_UL,
            seeding_target_cells=SEEDING_TARGET_CELLS,
            )
        if sp_weak_wells:
            logger.info(
                f"Excluding wells that required "
                f">{MAX_RELIABLE_ALIQUOT_UL:.1f} μL aliquots: "
                + ", ".join(sp_weak_wells)
                )
            excl_well_names.update(sp_weak_wells)
        conc_data = [
            val for (well_name, val) in dset.oplate_fc_data.items()
            if well_name not in excl_well_names
            ]
        conc_mean = np.average(conc_data)
        conc_std = np.std(conc_data)
        conc_cv = conc_std / conc_mean
        text_lines.append(f"{subfig_name}, {dset_name!r:32s}, "
                          f"{len(conc_data):2d}, "
                          f"{conc_mean:.1f}, {100. * conc_cv:.1f}%")
        plot_data.append(conc_data)
        plot_cats.append(subfig_name)
        plot_infos.append({
            'subfig_name': subfig_name,
            'x_pos': cat_idx,
            'dset_name': dset_name,
            'n_pts': len(conc_data),
            'n_bad_fcs': dset.n_bad_fcs,
            'n_man_excl': len(man_excl_well_names),
            'n_weak_src': len(sp_weak_wells),
            'conc_mean': conc_mean,
            'conc_cv': conc_cv,
            })

    # TODO: pull this out into analysis_common and have cci_gt use it too?
    make_violin_plot(ax, dataset=plot_data, raincloud=True,
                     positions=range(len(plot_data)),
                     showmeans=True)
    for plot_info in plot_infos:
        ax.text(
            plot_info['x_pos'],
            0.05,
            f"<keepsize>CV={100. * plot_info['conc_cv']:.1f}%"
                f" (n={plot_info['n_pts']})"
                #f"\n(n={plot_info['n_pts']})\n\n"
                #"Excluded:\n"
                #f"{plot_info['n_weak_src']} weak src\n"
                #f"{plot_info['n_man_excl']} bad img\n"
                #f"{plot_info['n_bad_fcs']} bad FCS"
                ,
            ha='center',
            transform=ax.get_xaxis_transform(),
            fontsize=LEGEND_FONT_SIZE,
            rotation=90
            )

    ax.set_xticks(range(len(plot_data)))
    ax.set_xticklabels(plot_cats)
    ax.set_ylabel("FC-measured conc. (cells/μL)")
    ax.set_ylim(0., 1.1 * max(itertools.chain(*plot_data)))

    groups_fig = PubFigure(
        fig_id=f"fig_2_panel_c",
        width_in=GROUPS_FIG_WIDTH_IN,
        height_in=GROUPS_FIG_HEIGHT_IN,
        n_rows=2,
        share_x=True
        )
    cci_splate_ax = groups_fig[0][0]
    fc_oplate_ax = groups_fig[1][0]
    dset = dsets[GROUPS_FIG_USE_DSET]
    excl_wells = (
        BEAD_WELLS[GROUPS_FIG_USE_DSET] + EXCLUDE_WELLS[GROUPS_FIG_USE_DSET])
    seeding_density_groups = {
        seeding_strength: [
            well_name
            for (well_name, well_seeding_strength)
            in dset.crazyplate_conc.items()
            if seeding_strength == well_seeding_strength
            ]
        for seeding_strength in sorted(
            set(dset.crazyplate_conc.values()), reverse=True)
        }
    splate_cci_group_data = group_data_by_crazyplate_group(
        data=dset.splate_cci_data,
        seeding_density_groups=seeding_density_groups,
        excl_wells=excl_wells
        )
    cci_splate_ax.hist(
        splate_cci_group_data,
        range=(0., HIST_X_MAX),
        bins=HIST_N_BINS,
        edgecolor='none',
        color=BAR_COLORS[:len(splate_cci_group_data)],
        stacked=True
        )
    oplate_fc_group_data = group_data_by_crazyplate_group(
        data=dset.oplate_fc_data,
        seeding_density_groups=seeding_density_groups,
        excl_wells=excl_wells
        )
    fc_oplate_ax.hist(
        oplate_fc_group_data,
        label=[
            f"{100. * x:.0f}%"
            for x in seeding_density_groups.keys()
            ],
        range=(0., HIST_X_MAX),
        bins=HIST_N_BINS,
        color=BAR_COLORS[:len(oplate_fc_group_data)],
        stacked=True
        )
    cci_splate_ax.set_xlabel("Input plate CCI measurement (cells/μL)")
    #cci_splate_ax.set_ylabel("# of samples")
    cci_splate_ax.yaxis.set_major_locator(
        mpl.ticker.MaxNLocator(5, integer=True))
    cci_splate_ax.set_title("a", x=0.05, y=1.0, pad=-12.).set_fontweight("bold")
    fc_oplate_ax.set_xlabel("Output plate FC measurement (cells/μL)")
    fc_oplate_ax.set_ylabel("# of samples")
    fc_oplate_ax.yaxis.set_major_locator(
        mpl.ticker.MaxNLocator(5, integer=True))
    fc_oplate_ax.legend(
        title="Src. well rel. seeding density",
        ncol=2,
        handlelength=1.,
        labelspacing=0.5,
        columnspacing=0.5,
        handletextpad=0.15,
        frameon=False
        )
    fc_oplate_ax.set_title("b", x=0.05, y=1.0, pad=-12.).set_fontweight("bold")

    text_output = TextOutput(fig_id=f"{fig_id}_stats", text_lines=text_lines)

    #sim_results = cci_norm_fc_sim.get_sim_results()
    #mc_counts = [
    #    y for result in sim_results for y in result.fc_counts.values()]
    #mc_concs = np.asarray([
    #    y * 1e-3 / cci_norm_fc_sim.FC_SAMPLE_VOL_ML for y in mc_counts])
    #mc_conc_mean = np.mean(mc_concs)
    #mc_conc_std = np.std(mc_concs)
    #mc_conc_cv = mc_conc_std / mc_conc_mean
    #fc_concs = plot_data[0] # TODO: select by dataset name or subfig name
    #fc_conc_mean = np.mean(fc_concs)
    #fc_conc_std = np.std(fc_concs)
    #fc_conc_cv = fc_conc_std / fc_conc_mean

    # mc_fig = PubFigure(
    #     fig_id=f"{fig_id}_mc",
    #     width_in=MC_FIG_WIDTH_IN,
    #     height_in=MC_FIG_HEIGHT_IN
    #     )
    # ax = mc_fig[0][0]
    # mc_hist_xmax = math.ceil(mc_conc_mean * 1.5 / 10.) * 10.
    # ax.hist(
    #     fc_concs,
    #     bins=HIST_N_BINS,
    #     range=(0., mc_hist_xmax),
    #     alpha=0.75,
    #     density=True,
    #     label=(
    #         f"FC measurement\n(mean={fc_conc_mean:.1f}, "
    #         f"CV={100. * fc_conc_cv:.1f}%, n={len(fc_concs)})"
    #         ),
    #     )
    # ax.hist(
    #     mc_concs,
    #     bins=HIST_N_BINS,
    #     range=(0., mc_hist_xmax),
    #     histtype='step',
    #     ec='red',
    #     lw=1.5,
    #     density=True,
    #     label=(
    #         f"MC predicted distribution\n(mean={mc_conc_mean:.1f}, "
    #         f"CV={100. * mc_conc_cv:.1f}%, n={len(mc_concs)/(10**3):.0f}k)"
    #         )
    #     )
    # ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.05))
    # ax.set_xlim(0., mc_hist_xmax)
    # ax.set_xlabel("Cell concentration (cells/μL)")
    # ax.set_ylabel("PD")
    # ax.set_ylim(0, .15)
    # ax.legend(loc="upper center")

    return (
        [
            overview_fig,
            groups_fig,
            # mc_fig
            ],
        [text_output]
        )
