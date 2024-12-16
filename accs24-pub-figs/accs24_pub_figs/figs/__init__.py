import os
import importlib
import logging
import pkgutil


logger = logging.getLogger(__name__)


__all__ = []
fig_mods = {}
for mod_info in pkgutil.iter_modules([os.path.dirname(__file__)]):
    module = importlib.import_module(
        "." + mod_info.name, "accs24_pub_figs.figs")
    __all__.append(mod_info.name)
    fig_mods[mod_info.name] = module
    logger.info(f"Loaded figure generator module {mod_info.name!r}")
