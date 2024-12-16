import itertools
import pickle
import logging
import numpy as np
from pathlib import Path
from typing import (
    Any, Callable, Collection, Dict, Iterable, List, Optional, Union)

import flowkit
import matplotlib

from ...figure import PubFigure, TextOutput
from ...analysis_common import (
    build_gating_strategy, count_sample, process_fcs_files)
#from ...mc_cell_suspension import CellSuspension
from ... import global_config
from .cci_offline import CciOfflineResult


logger = logging.getLogger(__name__)


DESCRIPTION = "CCI ground truth"
FIG_WIDTH_IN = 2.5
FIG_HEIGHT_IN = 2.25

DATA_DIR = "input_data"
USE_DATASET = "2023-10-26 cci gt"

SAMPLE_GROUP_WELLS = {
    'NC': [f'F{col}' for col in [1, 2, 3, 4]],
    '100k': [f'E{col}' for col in [1, 2, 3, 4]],
    '200k': [f'D{col}' for col in [1, 2, 3, 4]],
    '300k': [f'C{col}' for col in [1, 2, 3, 4]],
    '400k': [f'B{col}' for col in [1, 2, 3, 4]],
    '500k': [f'A{col}' for col in [1, 2, 3, 4]],
    }
CCI_SAMPLE_ORDER = ['500k', '400k', '300k', '200k', '100k', 'NC']
CCI_MIN_VALID_PIXELS = 0.7e6
USE_WELLS = [f"{row}{col}" for row in "ABCDEF" for col in [1]]
TUBE_VOL_CELLS_UL = 1000.
TUBE_VOL_BEADS_UL = 200.
BEAD_STOCK_BEADS_PER_UNIT = 0.52e5
BEAD_STOCK_UL_PER_UNIT = 50.
FC_SAMPLE_LOADING_TIME_S = 50.
FC_TIME_TOL_S = 1.
SIM_N_RUNS = 1000
SIM_STOCK_VOL_UL = 10e3 # FIXME
SIM_CCI_SAMPLE_VOL_UL = 30.
SIM_CCI_IMAGING_VOL_UL = 0.8 # FIXME
SIM_CCI_N_VALID_PIXELS = 0.7e6 # FIXME

WHISKER_COLOR = "#065b86"
BODY_COLOR = "#00a0dd"
FIT_LINE_COLOR = "#ff7d53"

# Feb 23 run
#USE_DATASET = "2023-02-23 cci gt"
#TUBE_VOL_CELLS_UL = 2000.
#TUBE_VOL_BEADS_UL = 100.


def simulate_cci_reads(n_samples: int, cell_conc: float, n_runs: int,
                       rand_gen: Optional[np.random.Generator] = None):
    cpps = []
    for i in range(n_runs):
        run_cpps = []
        stock = CellSuspension.from_conc(
            vol=SIM_STOCK_VOL_UL, conc=cell_conc, generator=rand_gen)
        n_samples_rem = n_samples
        while n_samples_rem > 0:
            ganged_sample = stock.take_sample(SIM_CCI_SAMPLE_VOL_UL * 8)
            for j in range(8):
                ind_sample = ganged_sample.take_sample(SIM_CCI_SAMPLE_VOL_UL)
                img_sample = ind_sample.take_sample(SIM_CCI_IMAGING_VOL_UL)
                run_cpps.append(img_sample.n_cells / SIM_CCI_N_VALID_PIXELS)
            n_samples_rem -= 8
        cpps.extend(run_cpps[:n_samples])
    return cpps


def get_sim_cci_cpps(group_name: str, n_samples: int, cell_conc: float,
                     rand_gen: np.random.Generator, n_runs: int = SIM_N_RUNS):
    call_args = locals()
    cache_path = (
        Path(__file__).parent.joinpath("cache")
        .joinpath(f"{group_name}.sim_cache.pickle")
        )
    match_info = {
        key: call_args[key]
        for key in ['group_name', 'n_samples', 'cell_conc', 'n_runs']
        }
    if cache_path.is_file():
        with cache_path.open("rb") as f:
            cache_data = pickle.load(f)
        if all(match_info[key] == cache_data[key] for key in match_info):
            logger.info(
                f"Using cached sim data for group {group_name!r} "
                f"({match_info})"
                )
            return cache_data['cpps']
    logger.info(
        f"Simulating CCI reads for group {group_name!r} ({match_info})")
    cpps = simulate_cci_reads(
            n_samples=n_samples,
            cell_conc=cell_conc,
            rand_gen=rand_gen,
            n_runs=n_runs
            )
    with cache_path.open("wb") as f:
        cache_data = {
            'group_name': group_name,
            'n_samples': n_samples,
            'cell_conc': cell_conc,
            'n_runs': n_runs,
            'cpps': cpps
            }
        pickle.dump(cache_data, f)
    return cpps


def mkoutputs_cci_gt(fig_id: str):
    dataset_dir = (
        Path(__file__).parent.joinpath(DATA_DIR).joinpath(USE_DATASET)
        )
    if not dataset_dir.is_dir():
        raise Exception(f"Bad dataset dir path: {str(dataset_dir)!r}")

    # Process CCI results
    cci_results_dir = dataset_dir.joinpath("cci_results")
    cci_results = [
        CciOfflineResult.read_from_file(path)
        for path in sorted(cci_results_dir.glob("*.json"))
        ]
    cci_sample_cycle = itertools.cycle(CCI_SAMPLE_ORDER)
    group_name_for_cci_scan_id = {
        scan_id: next(cci_sample_cycle)
        for scan_id in sorted({r.scan_id() for r in cci_results})
        }
    cci_rejects_lowva = [
        r.meas_id() for r in cci_results
        if r.n_valid_pixels < CCI_MIN_VALID_PIXELS]
    logger.info(
        f"Rejecting {len(cci_rejects_lowva)} CCI result(s) "
        "for image obstruction: "
        + ", ".join(f"{x!r}" for x in cci_rejects_lowva))
    cci_rejects_lowcount = [
        r.meas_id() for r in cci_results
        if group_name_for_cci_scan_id[r.scan_id()] != 'NC'
        and r.cells_per_ml < 10e3
        ]
    logger.info(
        f"Rejecting {len(cci_rejects_lowcount)} CCI result(s) "
        "for blurry/bad image: "
        + ", ".join(f"{x!r}" for x in cci_rejects_lowcount))
    cci_rejects = cci_rejects_lowva + cci_rejects_lowcount
    cci_results = [
        result for result in cci_results
        if result.meas_id() not in cci_rejects
        ]
    cci_results_by_group: Dict[str, List[CciOfflineResult]] = {
        group_name: [] for group_name in CCI_SAMPLE_ORDER
        }
    for result in cci_results:
        cci_results_by_group[
            group_name_for_cci_scan_id[result.scan_id()]].append(result)
    cpp_by_sample_name = {
        sample_name: [r.n_cells / r.n_valid_pixels for r in sub_results]
        for (sample_name, sub_results) in cci_results_by_group.items()
        }
    logger.info(f"Loaded {len(cci_results)} CCI results "
                f"from {str(cci_results_dir)!r}")

    # Process flow cytometry data
    fcs_data_dir = dataset_dir.joinpath("fcs_data")
    fcs_paths = fcs_data_dir.glob("*_*_*.[fF][cC][sS]")
    gating_strategy = build_gating_strategy()
    fcs_results = process_fcs_files(
        fcs_paths,
        cb_fn=lambda s: count_sample(
            sample=s,
            gating_strategy=gating_strategy,
            sample_loading_time_s=FC_SAMPLE_LOADING_TIME_S,
            time_tol_s=FC_TIME_TOL_S
            ),
        only_wells=USE_WELLS
        )
    logger.info(f"Analyzed {len(fcs_results)} FCS files "
                f"from {str(fcs_data_dir)!r}")

    fcs_results_for_group = {
        group_name: [
            result for result in fcs_results.values()
            if result['well_name'] in group_wells
            ]
        for (group_name, group_wells) in SAMPLE_GROUP_WELLS.items()
        }
    totals_for_group = {
        group_name: {
            key: sum(
                result[key] for result in fcs_results_for_group[group_name])
            for key in ['n_cells', 'n_beads']
            }
        for group_name in SAMPLE_GROUP_WELLS
        }
    bead_corrected_cell_concs = {}
    fc_sample_vols = {}
    for group_name, count_totals in totals_for_group.items():
        total_cell_count = count_totals['n_cells']
        total_bead_count = count_totals['n_beads']
        sample_vol = (
            (
                total_bead_count * BEAD_STOCK_UL_PER_UNIT
                / BEAD_STOCK_BEADS_PER_UNIT
                )
            * (1 + TUBE_VOL_CELLS_UL/TUBE_VOL_BEADS_UL)
            )
        seeded_cell_density = (
            (total_cell_count / sample_vol)
            * (TUBE_VOL_CELLS_UL + TUBE_VOL_BEADS_UL) / TUBE_VOL_CELLS_UL
            )
        bead_corrected_cell_concs[group_name] = seeded_cell_density
        fc_sample_vols[group_name] = sample_vol
        logger.info(
            f"FC sample group {group_name!r}: {total_cell_count} cells, "
            f"{total_bead_count} beads in {sample_vol:.1f} uL sample "
            f"=> {seeded_cell_density:.1f}k cells/mL"
            )

    fit_xs = []
    fit_ys = []
    for sample_name in SAMPLE_GROUP_WELLS.keys():
        if sample_name == 'NC':
            continue
        cpp = cpp_by_sample_name[sample_name]
        fit_ys.extend(cpp)
        fit_xs.extend([bead_corrected_cell_concs[sample_name]]*len(cpp))

    x_arr = np.array(fit_xs, dtype='float64')
    y_arr = np.array(fit_ys, dtype='float64')
    fit_slope = np.mean(x_arr.dot(y_arr)) / np.mean(x_arr.dot(x_arr))
    r_squared = 1. - np.mean((y_arr - x_arr * fit_slope)**2.) / np.var(y_arr)

    combo_fig = PubFigure(
        fig_id="fig_2_panel_a",
        width_in=FIG_WIDTH_IN,
        height_in=FIG_HEIGHT_IN,
        #n_rows=2,
        n_rows=1,
        gridspec_kw={
            'height_ratios': [
                1.0,
#                0.6,
                ],
            },
        share_x=True,
        )

    x_poss = [bead_corrected_cell_concs[k] for k in CCI_SAMPLE_ORDER]
    y_datas = [cpp_by_sample_name[k] for k in CCI_SAMPLE_ORDER]
    cvs = [np.std(y_data) / np.mean(y_data) for y_data in y_datas]

    violinplot_colls = combo_fig[0][0].violinplot(
        y_datas,
        positions=x_poss,
        showmeans=True,
        widths=50.,
        )
    for body in violinplot_colls['bodies']:
        body.set_facecolor(BODY_COLOR)
        body.set_alpha(0.55)
    for coll in [
            coll for (coll_name, coll) in violinplot_colls.items()
            if coll_name in ['cmeans', 'cmins', 'cmaxes', 'cbars', 'cmedians']
            ]:
        coll.set_linewidth(1.)
        coll.set_edgecolor(WHISKER_COLOR)

    label_params = [
        ((x_pos, np.mean(y_datas[i])), cvs[i])
        for (i, x_pos) in enumerate(x_poss)
        ]
    for xy, cv in label_params:
        if xy[1] < 1e-4:
            continue
        combo_fig[0][0].annotate(
            f"<keepsize>CV={cv * 100:.1f}%",
            fontsize=global_config.LEGEND_FONT_SIZE,
            xy=xy, xytext=(xy[0] + 15., xy[1] - 7.0e-5)
            )

    x_max = max(fit_xs) * 1.1
    combo_fig[0][0].plot(
        [0., x_max],
        [0., fit_slope * x_max],
        color=FIT_LINE_COLOR, lw=1.5, alpha=0.8,
        label=f"$y = ({1e-3*fit_slope/1e-10:.2f} \\times 10^{{-10}} "
              f"\\mathrm{{mL/pixel}}) \\cdot x$\n$(R^2 = {r_squared:.2f})$"
        )
    combo_fig[0][0].xaxis.set_tick_params(rotation=90)
    combo_fig[0][0].set_xticks(x_poss)
    combo_fig[0][0].ticklabel_format(axis='y', style='sci', scilimits=(-4, -4))
    combo_fig[0][0].set_xlabel(
        "Bead-corrected FC meas. (cells/μL)")
    combo_fig[0][0].set_xlim(-50., x_max + 125.)
    combo_fig[0][0].set_ylim(combo_fig[0][0].get_ylim()[0], 7e-4)
    combo_fig[0][0].set_ylabel("CCI reading (cells/valid px)")
    combo_fig[0][0].legend(
        loc="upper left",
        frameon=False,
        )
#    combo_fig[0][0].set_title("(A)", loc='left')

#    rand_gen = CellSuspension.new_rand_gen()
#    sim_cpp_data_by_sample_name = {
#        group_name: get_sim_cci_cpps(
#            group_name=group_name,
#            n_samples=len(cpp_by_sample_name[group_name]),
#            cell_conc=bead_corrected_cell_concs[group_name],
#            rand_gen=rand_gen
#            )
#        for group_name in SAMPLE_GROUP_WELLS
#        }
    cpp_cv_by_sample_name = {}
#    sim_cpp_cv_by_sample_name = {}
    text_lines = ["group, cci_n, cci_mean, cci_cv, "
                  "fc_n_cells, fc_n_beads, fc_vol, fc_conc, "
#                  "sim_cci_mean, sim_cci_cv"]
                  ]
    for group_name in SAMPLE_GROUP_WELLS:
        cpp_data = cpp_by_sample_name[group_name]
        cpp_mean = np.mean(cpp_data)
        cpp_std = np.std(cpp_data)
        cpp_cv = cpp_std / cpp_mean
#        sim_cpp_data = sim_cpp_data_by_sample_name[group_name]
#        sim_cpp_mean = np.mean(sim_cpp_data)
#        sim_cpp_std = np.std(sim_cpp_data)
#        sim_cpp_cv = sim_cpp_std / sim_cpp_mean
        text_lines.append(f"{group_name:4s}, "
                          f"{len(cpp_by_sample_name[group_name]):3d}, "
                          f"{cpp_mean:.3e}, "
                          f"{100. * cpp_cv:5.1f}%, "
                          f"{totals_for_group[group_name]['n_cells']:6d}, "
                          f"{totals_for_group[group_name]['n_beads']:6d}, "
                          f"{fc_sample_vols[group_name]:.1f}, "
                          f"{bead_corrected_cell_concs[group_name]:6.1f}, "
#                          f"{sim_cpp_mean:.3e}, "
#                          f"{100. * sim_cpp_cv:5.1f}%"
                          )
        cpp_cv_by_sample_name[group_name] = cpp_cv
#        sim_cpp_cv_by_sample_name[group_name] = sim_cpp_cv

#    cvs_group_order = [x for x in SAMPLE_GROUP_WELLS if x != 'NC']
#    combo_fig[1][0].plot(
#        [
#            bead_corrected_cell_concs[group_name] 
#            for group_name in cvs_group_order
#            ],
#        [
#            100. * cpp_cv_by_sample_name[group_name]
#            for group_name in cvs_group_order
#            ],
#        marker="o",
#        label="CCI data"
#        )
#    combo_fig[1][0].plot(
#        [
#            bead_corrected_cell_concs[group_name]
#            for group_name in cvs_group_order
#            ],
#        [
#            100. * sim_cpp_cv_by_sample_name[group_name]
#            for group_name in cvs_group_order
#            ],
#        marker="o",
#        label="MC predicted sampling error"
#        )
#    combo_fig[1][0].set_xlabel("Cell concentration (cells/μL)")
#    combo_fig[1][0].set_ylabel("CV (%)")
#    combo_fig[1][0].set_ylim(0., 15.)
#    combo_fig[1][0].legend(frameon=False)
#    combo_fig[1][0].set_title("(B)", loc='left')

    text_output = TextOutput(fig_id=f"{fig_id}_stats", text_lines=text_lines)

    return [combo_fig], [text_output]
