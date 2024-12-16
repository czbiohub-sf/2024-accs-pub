from typing import Final


IMG_FMT_PRINT = "svg"
IMG_DPI_PRINT = 600.  # Doesn't matter for eps/pdf unless there's raster stuff in the figure
IMG_COMPRESSION_PRINT = None
IMG_FMT_DISPLAY = "png"
IMG_DPI_DISPLAY = 150.
WIDTH_MAX_IN = 2250/300
HEIGHT_MAX_IN = 2625/300
TEXT_FONT_PATH = "fonts/HKGrotesk-Regular.ttf" # FIXME
TITLE_FONT_SIZE = 11.
XYLABEL_FONT_SIZE = 8.
TICKLABEL_FONT_SIZE = 8.
LEGEND_FONT_SIZE = 7.
DEFAULT_FONT_SIZE = 9.
PURPLESTAMP_FONT_PATH = "fonts/RobotoMono-Regular.ttf"
PURPLESTAMP_FONT_SIZE = 22
FRAME_LINEWIDTH_PT = 0.5
TICKS_LINEWIDTH_PT = 0.5
FIG_LAYOUT: Final = 'constrained'
OUT_DIR_PRINT = "out_print"
OUT_DIR_DISPLAY = "out_disp"
OUT_DIR_TEXT = "out_text"
BRAND_COLORS = {
    'lavendish': "#c395d7",
    'rustish': "#ac4b3a",
    'orangish': "#ff7d53",
    'purplish': "#5b26f1",
    'sf_cyan_bright': "#00a0dd",
    'sf_cyan': "#0d7cb5",
    'sf_cyan_med': "#065b86",
    'net_gray_700': "#bdc3c6",
    }
