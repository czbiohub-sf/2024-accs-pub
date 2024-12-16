from . import analysis as _analysis
DESCRIPTION = _analysis.DESCRIPTION

fig_makers = {
    'fov_scores': _analysis.mkoutputs_fov_scores,
    }
