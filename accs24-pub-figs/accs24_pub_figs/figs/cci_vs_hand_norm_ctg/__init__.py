from . import analysis as _analysis
DESCRIPTION = _analysis.DESCRIPTION

fig_makers = {
    'cci_vs_hand_norm_ctg': _analysis.mkoutputs_cci_vs_hand_norm,
    }
