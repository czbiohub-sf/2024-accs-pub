import itertools
import logging
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union

import matplotlib as mpl
import matplotlib.figure
import PIL.Image
import PIL.ImageDraw
import PIL.ImageFont

from . import global_config

PathOrStr = Union[Path, str]


logger = logging.getLogger(__name__)

for font_path in mpl.font_manager.findSystemFonts(
        fontpaths=[
            Path(__file__).parent.joinpath(
                global_config.TEXT_FONT_PATH).parent
            ]  # FIXME
        ):
    mpl.font_manager.fontManager.addfont(font_path)


class Outputtable:
    fig_id: str

    def _get_out_path(self, out_dir: PathOrStr,
                      ext: str, create_dir: bool) -> Path:
        out_dir_path = Path(out_dir)
        if create_dir and not out_dir_path.is_dir():
            out_dir_path.mkdir(parents=True)
        return out_dir_path.joinpath(f"{self.fig_id}.{ext}")


class TextOutput(Outputtable):
    def __init__(self, fig_id: str, text_lines: Iterable[str]):
        self.fig_id = fig_id
        self.text_lines = list(text_lines)

    def get_textfile_path(
            self,
            out_dir: PathOrStr = Path(global_config.OUT_DIR_TEXT),
            create_dir: bool = True
            ) -> Path:
        return self._get_out_path(
            out_dir=out_dir,
            ext="txt",
            create_dir=create_dir
            )

    def save_text(self, path: Optional[PathOrStr] = None):
        path = (
            Path(path) if path is not None
            else self.get_textfile_path()
            ).resolve()
        with path.open("w") as f:
            f.write("\n".join(self.text_lines) + "\n")
        logger.info(
            f"Saved text output for {self.fig_id!r} to {str(path)!r}")


class PubFigureBase(Outputtable):
    def __init__(self, fig_id: str, width_in: float, height_in: float):
        self.fig_id = fig_id
        self.width_in = width_in
        self.height_in = height_in
        self.fig = mpl.figure.Figure(
            figsize=(width_in, height_in),
            layout=global_config.FIG_LAYOUT,
            )
        self._formatting_update_done = False

    @staticmethod
    def _update_axs_formatting(axs):
        # TODO: set legend border style
        for ax in axs:
            legend = ax.get_legend()
            if legend is not None:
                for text in [legend.get_title()] + legend.get_texts():
                    text.set_fontsize(global_config.LEGEND_FONT_SIZE)
                legend.get_frame().set_alpha(1.)
            for x in [ax._left_title, ax._right_title, ax.title]:  # gross :(
                x.set_fontsize(global_config.TITLE_FONT_SIZE)
            for x in [ax.xaxis, ax.yaxis]:
                x.label.set_fontsize(global_config.XYLABEL_FONT_SIZE)
                x.get_offset_text().set_fontsize(
                    global_config.TICKLABEL_FONT_SIZE)
                x.set_tick_params(labelsize=global_config.TICKLABEL_FONT_SIZE)
            for x in ax.spines.values():
                x.set_linewidth(global_config.FRAME_LINEWIDTH_PT)
            for x in itertools.chain(
                    ax.xaxis.get_ticklines(), ax.yaxis.get_ticklines()):
                x.set_linewidth(global_config.TICKS_LINEWIDTH_PT)
                x.set_markeredgewidth(global_config.TICKS_LINEWIDTH_PT)

    def _update_fig_formatting(self, fig):
        if self._formatting_update_done:
            return
        # FIXME clean this section up after the hax
        for item in fig.findobj(mpl.text.Text):
            old_fontsize = item.get_fontsize()
            item.set_family(["HK Grotesk", "Arial"]) # FIXME
            if item.get_text().startswith("<keepsize>"):
                item.set_text(item.get_text().rsplit(">", 1)[-1])
                item.set_fontsize(old_fontsize)
            else:
                item.set_fontsize(global_config.DEFAULT_FONT_SIZE)
        self._update_axs_formatting(fig.get_axes())
        self._formatting_update_done = True

    def _add_purple_stamp(self, img_path: PathOrStr, stamp_text: str):
        font = PIL.ImageFont.truetype(
            str(Path(__file__).parent.joinpath(
                global_config.PURPLESTAMP_FONT_PATH)),
            global_config.PURPLESTAMP_FONT_SIZE
            )
        with PIL.Image.open(img_path) as img:
            text_bbox = PIL.ImageDraw.Draw(img).textbbox(
                xy=(0, 0), text=stamp_text, font=font)
            text_layer = PIL.Image.new(
                mode="RGBA",
                size=(text_bbox[2] + 20, text_bbox[3] + 10),
                color=(255, 0, 255, 127)
                )
            draw = PIL.ImageDraw.Draw(text_layer)
            draw.text(xy=(10, 5), text=stamp_text,
                      color=(255, 255, 255, 255), font=font)
        img.alpha_composite(text_layer, dest=(20, 20))
        img.convert("RGB").save(img_path)

    def get_display_img_path(
            self,
            out_dir: PathOrStr = Path(global_config.OUT_DIR_DISPLAY),
            create_dir: bool = True
            ):
        return self._get_out_path(
            out_dir=out_dir,
            ext=global_config.IMG_FMT_DISPLAY,
            create_dir=create_dir
            )

    def get_print_img_path(
            self,
            out_dir: PathOrStr = Path(global_config.OUT_DIR_PRINT),
            create_dir: bool = True
            ):
        return self._get_out_path(
            out_dir=out_dir,
            ext=global_config.IMG_FMT_PRINT,
            create_dir=create_dir
            )

    def _before_saving(self):
        pass

    def save_img(self, path: PathOrStr, dpi: float, format: str,
                 pil_kwargs: Optional[Dict[str, Any]] = None):
        self._before_saving()
        out_path = str(Path(path).resolve())
        savefig_kwargs = {'pil_kwargs': pil_kwargs} if pil_kwargs else {}
        self.fig.savefig(
            out_path, dpi=dpi, format=format, **savefig_kwargs)

    def save_display_img(self, path: Optional[PathOrStr] = None,
                         purple_stamp: bool = True):
        path_resolved = (
            Path(path) if path is not None
            else self.get_display_img_path()
            ).resolve()
        self.save_img(
            str(path_resolved),
            dpi=global_config.IMG_DPI_DISPLAY,
            format=global_config.IMG_FMT_DISPLAY
            )
        if purple_stamp:
            print_fname = self.get_print_img_path(create_dir=False).name
            self._add_purple_stamp(
                path_resolved,
                (
                    f"print: {print_fname}\n"
                    f"{self.width_in:.2f}\" x {self.height_in:.2f}\""),
                )
        logger.info(
            f"Saved preview image for {self.fig_id!r} to {str(path_resolved)!r}")
        return path_resolved

    def save_print_img(self, path: Optional[PathOrStr] = None):
        path_resolved = (
            Path(path) if path is not None
            else self.get_print_img_path()
            ).resolve()
        pil_kwargs = None
        if global_config.IMG_COMPRESSION_PRINT:
            pil_kwargs = {'compression': global_config.IMG_COMPRESSION_PRINT}
        self.save_img(
            path_resolved,
            dpi=global_config.IMG_DPI_PRINT,
            format=global_config.IMG_FMT_PRINT,
            pil_kwargs=pil_kwargs
            )
        logger.info(
            f"Saved print image for {self.fig_id!r} to {str(path_resolved)!r}")
        return path_resolved

    def save_imgs(self):
        disp_path = self.save_display_img()
        print_path = self.save_print_img()
        return (disp_path, print_path)


class PubSuperFigure(PubFigureBase):
    def __init__(self, fig_id: str, width_in: float, height_in: float,
                 n_rows: int = 1, n_cols: int = 1):
        super().__init__(fig_id=fig_id, width_in=width_in, height_in=height_in)
        self.subfigs = self.fig.subfigures(
            nrows=n_rows, ncols=n_cols, squeeze=False)

    def __getitem__(self, idx):
        return self.subfigs[idx]

    def _before_saving(self):
        for subfig in itertools.chain(*self.subfigs):
            self._update_fig_formatting(subfig)


class PubFigure(PubFigureBase):
    def __init__(self, fig_id: str, width_in: float, height_in: float,
                 n_rows: int = 1, n_cols: int = 1,
                 share_x: bool = False, share_y: bool = False,
                 gridspec_kw: Optional[Dict[str, Any]] = None):
        super().__init__(fig_id=fig_id, width_in=width_in, height_in=height_in)
        self.axs = self.fig.subplots(
            nrows=n_rows, ncols=n_cols,
            sharex=share_x, sharey=share_y,
            squeeze=False,
            gridspec_kw=gridspec_kw or {},
            )

    def __getitem__(self, idx):
        return self.axs[idx]

    def _before_saving(self):
        self._update_fig_formatting(self.fig)
