import functools
import itertools
import logging
from pathlib import Path
import random
import re
from typing import Dict, List, Optional, Set, Tuple

import aicsimageio
import matplotlib as mpl
import numpy as np

from ...analysis_common import PathOrStr, MissingInputFile
from ...data_loaders import FovScoreLogCsv
from ...figure import PubFigure, TextOutput
from ... import global_config


logger = logging.getLogger(__name__)


DESCRIPTION = "Example Z-projections from OpenCell acquisitions from human- and ACCS-seeded wells"
FIG_WIDTH_IN = 7.0
FIG_HEIGHT_IN = 2.5
N_ROWS = 2
N_COLS = 6

DATA_DIR = "input_data"
IMG_FNAME_GLOB_FMT = "*MMStack_*-{well_name}-{site_id}.ome.tif"
SELECTED_DSETS = ["PML0435", "PML0442"]
DSET_LABELS = {
    'PML0435': "Manually-\nseeded plate",
    'PML0442': "ACCS-seeded\nplate"
    }
N_TOP_WELLS = 10
WELL_RANKING_N_TOP_SITES = 4
MISSING_SCORE_VAL = -1
CH_DAPI = 0
CH_GFP = 1
MPL_CMAP = 'gray'

FRACTILES = [0.10, 0.5, 0.90]
FRACTILE_LABELS = [
    "10th percentile FOV score",
    "Median FOV score",
    "90th percentile FOV score"
    ]
RANDOM_SEED = 499


class ImageProcessing:
    def __init__(self, data_dir: PathOrStr,
                 cache_dir: Optional[PathOrStr] = None):
        self.cache_dir_path = (
            Path(cache_dir) if cache_dir is not None else None)
        self.data_dir_path = Path(data_dir)
        self._last_img_path = None
        self._last_img_data = None

    def get_all_run_ids(self):
        return [x.name for x in self.data_dir_path.iterdir() if x.is_dir()]

    def get_raw_image_glob(self, run_id: str, well_name: str, site_id: int
                           ) -> Tuple[Path, str]:
        img_dir_path = self.data_dir_path / run_id / "raw_data"
        fname_glob = IMG_FNAME_GLOB_FMT.format(
            well_name=well_name, site_id=site_id)
        return img_dir_path, fname_glob

    def get_raw_image_path(self, run_id: str, well_name: str, site_id: int
                           ) -> Path:
        img_dir, fname_glob = self.get_raw_image_glob(
            run_id=run_id, well_name=well_name, site_id=site_id)
        if not img_dir.is_dir():
            raise Exception(f"Missing/invalid image data dir for {run_id!r}: "
                            f"{str(img_dir)!r}")
        try:
            img_path = next(img_dir.glob(fname_glob))
        except StopIteration:
            raise MissingInputFile(
                f"Can't find image file for run {run_id!r} well {well_name!r} "
                f"site {site_id}: {str(img_dir / fname_glob)!r}"
                )
        return img_path

    def get_zstack_cache_path(
            self, run_id: str, well_name: str, site_id: int, ch_id: int
            ) -> Optional[Path]:
        if self.cache_dir_path is None:
            return None
        return (
            self.cache_dir_path
            / f"zproj_{run_id}_{well_name}_{site_id}_{ch_id}.tif"
            )

    def generate_zproj(
            self, run_id: str, well_name: str, site_id: int, ch_id: int):
        raw_img_path = self.get_raw_image_path(
            run_id=run_id, well_name=well_name, site_id=site_id)
        if (self._last_img_path is not None
                and raw_img_path == self._last_img_path):
            inp_img = self._last_img
        else:
            logger.info(f"Loading raw image file: {str(raw_img_path)!r}")
            inp_img = aicsimageio.AICSImage(raw_img_path, chunk_dims="CZYX")
            self._last_img_path = raw_img_path
            self._last_img = inp_img
        inp_img_data = inp_img.get_image_data("ZYX", C=ch_id)
        logger.info(
            f"Running z-projection for {raw_img_path.name!r} channel {ch_id}")
        zproj_img = np.max(inp_img_data, axis=0)
        return zproj_img

    def get_zproj(self, run_id: str, well_name: str, site_id: int, ch_id: int,
                  force_generate: bool = False, save_to_cache: bool = True):
        cache_path = self.get_zstack_cache_path(
            run_id=run_id,
            well_name=well_name,
            site_id=site_id,
            ch_id=ch_id
            )
        if (cache_path is not None) and not force_generate:
            if cache_path.is_file():
                logger.info(
                    f"Loading cached z-projection image for run {run_id!r} "
                    f"well {well_name!r} site {site_id} channel {ch_id}")
                return (aicsimageio.AICSImage(cache_path)
                    .get_image_data("YX", C=0, Z=0))
        logger.info(
            f"Generating new z-projection image for run {run_id!r} "
            f"well {well_name!r} site {site_id} channel {ch_id}")
        zproj_img = self.generate_zproj(
            run_id=run_id, well_name=well_name, site_id=site_id, ch_id=ch_id)
        if save_to_cache and (cache_path is not None):
            aicsimageio.AICSImage(zproj_img).save(cache_path)
        return zproj_img

    def convert_img_to_np(self, img):
        return np.array(self.ij.py.from_java(img), dtype='uint16')


class Statisticking:
    def __init__(self, data_dir: PathOrStr):
        self.data_dir_path = Path(data_dir)
        fov_score_csv_paths = {
            subdir.name: subdir / "logs" / "fov_scoring" / "fov-score-log.csv"
            for subdir in self.data_dir_path.iterdir()
            }
        self.fov_score_data = {
            key: FovScoreLogCsv.from_path(path)
            for (key, path) in fov_score_csv_paths.items()
            if path.is_file()
            }
        logger.info(
            "Loaded FOV scores for datasets: "
            + ", ".join(f"{key!r}" for key in self.fov_score_data)
            )
        self._sites_with_images: Dict[str, List[Tuple[str, int]]] = {}
        self._image_paths: Dict[str, List[Path]] = {}
        self._sorted_scores: Dict[str, List[float]] = {}
        self.random = random.Random()
        self.random.seed(RANDOM_SEED, version=2)

    def get_all_image_paths_for_run(self, run_id: str) -> List[Path]:
        if run_id in self._image_paths:
            return list(self._image_paths[run_id])
        img_fname_glob = IMG_FNAME_GLOB_FMT.format(well_name="*", site_id="*")
        self._image_paths[run_id] = list(
            (self.data_dir_path / run_id / "raw_data").glob(img_fname_glob)
            )
        if not self._image_paths[run_id]:
            logger.warning(f"There don't seem to be images for run {run_id!r}")
        return list(self._image_paths[run_id])

    def get_sites_with_images(self, run_id: str) -> List[Tuple[str, int]]:
        if run_id in self._sites_with_images:
            return list(self._sites_with_images[run_id])
        all_image_paths = self.get_all_image_paths_for_run(run_id)
        pat = re.compile(
            r'.*-(?P<well_name>[a-zA-Z0-9]+)-(?P<site_id>\d+)\.ome\.tif')
        results = []
        for img_path in all_image_paths:
            match = pat.match(img_path.name)
            if match is None:
                raise Exception(
                    "Couldn't parse image filename for some reason: "
                    + repr(img_path.name)
                    )
            well_name = match.group('well_name')
            site_id = int(match.group('site_id'))
            results.append((well_name, site_id))
        self._sites_with_images[run_id] = list(results)
        return results

    def get_score_for_site(self, run_id: str, well_name: str, site_id: int
                           ) -> Optional[float]:
        return self.fov_score_data[run_id].records[well_name][site_id].score

    def get_scores_for_selected_sites(
            self, run_id: str, from_available_images: bool = False
            ) -> dict[tuple[str, int], Optional[float]]:
        if from_available_images:
            return self.get_scores_for_sites_with_images(run_id)
        else:
            return {
                (well_name, site_id): score
                for well_name in self.fov_score_data[run_id].records
                for (site_id, score) in self.get_top_n_sites_for_well(
                    run_id, well_name)
                }

    def get_scores_for_sites_with_images(
            self, run_id: str
            ) -> dict[tuple[str, int], Optional[float]]:
        sites_with_images = self.get_sites_with_images(run_id)
        return {
            (well_name, site_id):
                self.get_score_for_site(run_id, well_name, site_id)
            for (well_name, site_id) in sites_with_images
            }

    def get_sites_closest_to_score(
            self, run_id: str, score: float, n_sites: int = 1
            ) -> List[Tuple[Tuple[str, int], Optional[float]]]:
        score_list = list(
            self.get_scores_for_selected_sites(run_id).items())

        def key_fn(item):
            (well_name, site_id), item_score = item
            return abs(score - item_score)
        score_list.sort(key=key_fn)
        score_0 = score_list[0][1]
        identical_score_list = [
            item for item in score_list if item[1] == score_0]
        if len(identical_score_list) >= n_sites:
            self.random.shuffle(identical_score_list)
            return identical_score_list[:n_sites]
        return score_list[:n_sites]

    def get_score_for_fractile(self, run_id: str, fractile: float):
        if run_id not in self._sorted_scores:
            self._sorted_scores[run_id] = sorted(
                score for (site_addr, score)
                in self.get_scores_for_selected_sites(run_id).items()
                )
        scores = self._sorted_scores[run_id]
        if not scores:
            raise Exception(f"No scores to pick from for run {run_id!r}??")
        return np.quantile(scores, fractile)

    def get_sites_closest_to_fractile(
            self, run_id: str, fractile: float, n_sites: int = 1
            ) -> List[Tuple[Tuple[str, int], float]]:
        fractile_score = self.get_score_for_fractile(run_id, fractile)
        return self.get_sites_closest_to_score(
            run_id=run_id, score=fractile_score, n_sites=n_sites)

    def get_top_n_sites_for_well(
            self,
            run_id: str,
            well_name: str,
            n_sites: int = WELL_RANKING_N_TOP_SITES,
            ) -> List[Tuple[int, float]]:
        site_scores = [
            (rec.site_id, rec.score)
            for rec in self.fov_score_data[run_id].records[well_name].values()
            if rec.score is not None
            ]
        site_scores.sort(key=lambda x: x[1], reverse=True)
        return site_scores[:n_sites]

    def get_top_n_wells(
            self,
            run_id: str,
            n_wells: int = N_TOP_WELLS,
            n_top_sites_per_well: int = WELL_RANKING_N_TOP_SITES,
            missing_score_val: float = MISSING_SCORE_VAL,
            ) -> List[Tuple[str, float]]:
        top_n_sites_for_well = {
            well_name: self.get_top_n_sites_for_well(
                run_id, well_name, n_sites=n_top_sites_per_well)
            for well_name in self.fov_score_data[run_id].records
            }
        well_avg_scores: List[Tuple[str, float]] = []
        for well_name in top_n_sites_for_well:
            scores = [
                score for (site_id, score) in top_n_sites_for_well[well_name]]
            while len(scores) < n_top_sites_per_well:
                scores.append(missing_score_val)
            well_avg_scores.append((well_name, np.mean(scores)))
        well_avg_scores.sort(key=lambda x: x[1], reverse=True)
        return well_avg_scores[:n_wells]


def mkoutputs_zstack_gallery(fig_id: str
                             ) -> Tuple[List[PubFigure], List[TextOutput]]:
    data_dir_path = Path(__file__).parent / DATA_DIR
    cache_dir_path = Path(__file__).parent / "cache"
    proc = ImageProcessing(data_dir=data_dir_path, cache_dir=cache_dir_path)
    stats = Statisticking(data_dir=data_dir_path)

    for run_id in SELECTED_DSETS:
        logger.info(
            f"Top-scoring {N_TOP_WELLS} wells (based on sum of top "
            f"{WELL_RANKING_N_TOP_SITES} site scores) for {run_id!r}: "
            + ", ".join(
                f"{well_name!r}({score:.2f})"
                for (well_name, score)
                in stats.get_top_n_wells(run_id)
                )
            )

    fig = PubFigure(
        fig_id=fig_id,
        width_in=FIG_WIDTH_IN,
        height_in=FIG_HEIGHT_IN,
        n_rows=N_ROWS,
        n_cols=N_COLS
        )

    reqd_img_relpaths = []
    cols_per_dset = N_COLS // len(FRACTILES)
    for row_idx, run_id in enumerate(SELECTED_DSETS):
        col_idx = 0
        for fractile_label, fractile in zip(FRACTILE_LABELS, FRACTILES):
            sites = stats.get_sites_closest_to_fractile(
                run_id, fractile, n_sites=cols_per_dset)
            logger.info(
                f"Selected sites for run {run_id!r} fractile {fractile:0.2f}: "
                + ", ".join(
                    f"{well_name}-{site_id}(score={score:0.2f})"
                    for ((well_name, site_id), score) in sites
                    )
                )
            for i, ((well_name, site_id), score) in enumerate(sites):
                if not i and not row_idx:
                    fig.fig.text(
                        0.0, 1.08,
                        f"<keepsize>{fractile_label}",
                        transform=fig[0][col_idx].transAxes,
                        ha="left", va="bottom"
                        )
                    fig.fig.add_artist(
                        mpl.patches.ConnectionPatch(
                            (0., 1.05),
                            (1., 1.05),
                            coordsA="axes fraction",
                            coordsB="axes fraction",
                            axesA=fig[0][col_idx],
                            axesB=fig[0][col_idx + len(sites) - 1],
                            lw=1.0, ec="black",
                            )
                        )

                dapi_img = proc.get_zproj(
                    run_id=run_id,
                    well_name=well_name,
                    site_id=site_id,
                    ch_id=CH_DAPI
                    )
                #gfp_img = proc.get_zproj(
                #    run_id=run_id,
                #    well_name=well_name,
                #    site_id=site_id,
                #    ch_id=CH_GFP
                #    )
                plot_img = dapi_img
                ax = fig[row_idx][col_idx]
                ax.imshow(
                    plot_img,
                    vmin=np.quantile(plot_img, 0.2),
                    vmax=np.quantile(plot_img, 0.95),
                    interpolation='antialiased',
                    interpolation_stage='rgba',
                    cmap=MPL_CMAP
                    )
                ax.text(
                    0.05, 0.05,
                    f"<keepsize>score:\n{score:.2f}",
                    color="white",
                    transform=ax.transAxes,
                    ha="left",
                    fontsize=global_config.TICKLABEL_FONT_SIZE,
                    bbox={
                        'facecolor': "#303030",
                        'edgecolor': 'none',
                        'pad': 1.
                        }
                    )
                ax.set_xticks([])
                ax.set_yticks([])
                if not col_idx:
                    ax.annotate(
                        f"<keepsize>{DSET_LABELS[run_id]}",
                        xy=(0, 0.5),
                        xytext=(-1.2, 0),
                        xycoords="axes fraction",
                        textcoords="offset fontsize",
                        ha="center",
                        va="center",
                        rotation=90,
                        #fontsize=global_config.XYLABEL_FONT_SIZE,
                        )
                raw_img_dir, raw_img_glob = proc.get_raw_image_glob(
                    run_id=run_id, well_name=well_name, site_id=site_id)
                reqd_img_relpaths.append(str(
                    raw_img_dir.relative_to(data_dir_path) / raw_img_glob
                    ))

                col_idx += 1
            
    # TODO maybe do something less dumb
    layout_enginge = fig.fig.get_layout_engine()
    rect = layout_enginge.get()['rect']
    rect = (rect[0] + 0.005, rect[1], rect[2] - 0.01, rect[3] - 0.06)  # lol
    layout_enginge.set(rect=rect)

    text_output = TextOutput(
        f"{fig_id}_reqd_imgs",
        ["Required image files: "] + reqd_img_relpaths
        )

    return [fig], [text_output]
