import logging

from . import figs as figs_pkg


logger = logging.getLogger(__name__)


def generate_figs(fig_mod_name: str):
    for fig_name, fig_maker in figs_pkg.fig_mods[fig_mod_name].fig_makers.items():
        logger.info(f"From {fig_mod_name!r} running {fig_name!r}...")
        figs, texts = fig_maker(fig_name)
        logger.info(f"Generated {len(figs)} figure(s) from {fig_name!r}: " + ", ".join(f"{fig.fig_id!r}" for fig in figs))
        logger.info(f"Generated {len(texts)} text report(s) from {fig_name!r}: " + ", ".join(f"{txt.fig_id!r}" for txt in texts))
        for fig in figs:
            fig.save_imgs()
        for text in texts:
            text.save_text()


def generate_all_figs():
    for mod_name in figs_pkg.fig_mods.keys():
        generate_figs(mod_name)
