import logging
from pathlib import Path
from typing import Any, Callable, Collection, Dict, Iterable, Optional, Union

import flowkit


CELLS_GATE_DIMENSIONS = ('FSC-A', 'SSC-A')
CELLS_GATE_VERTICES = [
    (47e3, 15e3),
    (79e3, 9.4e3),
    (229e3, 29e3),
    (451e3, 95e3),
    (459e3, 194e3),
    (348e3, 228e3),
    (169e3, 140e3),
    (44e3, 37e3),
    ]
BEADS_GATE_RANGES = {
    'BB700-H': (2e3, 2e4),
    'BB515-H': (2e3, 4e4),
    }


logger = logging.getLogger(__name__)


PathOrStr = Union[Path, str]


class MissingInputFile(Exception):
    pass


def standardize_well_name(x: str) -> str:
    row = x[0]
    col = int(x[1:])
    return f"{row}{col}"


def process_fcs_files(paths: Iterable[PathOrStr], cb_fn: Callable[[Path], Any],
                      only_wells: Optional[Collection[str]] = None
                      ) -> Dict[str, Any]:
    return {
        path.name: cb_fn(sample)
        for path, sample in (
            (Path(fcs_path), flowkit.Sample(fcs_path)) for fcs_path in paths
            )
        if only_wells is None or (
            standardize_well_name(sample.get_metadata()['well id'])
            in only_wells
            )
        }


def build_gating_strategy() -> flowkit.GatingStrategy:
    beads_gate = flowkit.gates.RectangleGate(
        'beads',
        dimensions=[
            flowkit.Dimension(ch_name, range_min=ch_rmin, range_max=ch_rmax)
            for ch_name, (ch_rmin, ch_rmax) in BEADS_GATE_RANGES.items()
            ]
        )
    cells_base_gate = flowkit.gates.PolygonGate(
        'cellsBase',
        dimensions=[flowkit.Dimension(x) for x in CELLS_GATE_DIMENSIONS],
        vertices=list(CELLS_GATE_VERTICES)
        )
    cells_gate = flowkit.gates.BooleanGate(
        'cells',
        bool_type='AND',
        gate_refs=[
            {'ref': "cellsBase", 'path': ("root", ), 'complement': False},
            {'ref': "beads", 'path': ("root", ), 'complement': True}
            ]
        )
    gating_strategy = flowkit.GatingStrategy()
    gating_strategy.add_gate(beads_gate, ('root', ))
    gating_strategy.add_gate(cells_base_gate, ('root', ))
    gating_strategy.add_gate(cells_gate, ('root', 'cellsBase'))
    return gating_strategy


def check_fcs_time_data(sample: flowkit.Sample, sample_loading_time_s: float,
                        time_tol_s: float) -> bool:
    times = sample.get_channel_events(
        sample.get_channel_index('Time'), source='raw')
    if (
            max(times) > sample_loading_time_s + time_tol_s
            or min(times) > time_tol_s
            ):
        logger.warning(
            f"Possible bad FC sample (suspicious time data): {sample.id!r}")
        return True
    return False


def count_sample(sample: flowkit.Sample,
                 gating_strategy: flowkit.GatingStrategy,
                 sample_loading_time_s: float,
                 time_tol_s: float,
                 skip_wells: Collection[str] = ()) -> Dict[str, Any]:
    well_name = standardize_well_name(sample.get_metadata()['well id'])
    info = {
        'well_name': well_name,
        'analyzed': False,
        'suspicious': False
        }
    if well_name in skip_wells:
        return info
    info.update({
        'suspicious': check_fcs_time_data(
            sample,
            sample_loading_time_s=sample_loading_time_s,
            time_tol_s=time_tol_s
            ),
        'analyzed': True,
        })
    results = gating_strategy.gate_sample(sample)
    info.update({
        f'n_{gate_name}': results.get_gate_count(gate_name)
        for gate_name in ['cells', 'beads']
        })
    return info
