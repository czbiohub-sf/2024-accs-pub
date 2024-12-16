from dataclasses import dataclass
import numpy as np
from pathlib import Path
from typing import Dict, List
import xml.etree.ElementTree

from ...data_loaders import NpfRunLog


DATA_DIR = "input_data"


@dataclass(frozen=True)
class LumDataPt:
    ps_name: str
    well_name: str
    sample_idx: int
    time: float
    lum: int


def load_softmax_xml(f) -> Dict[str, Dict[str, List[LumDataPt]]]:
    root = xml.etree.ElementTree.parse(f)
    data_by_section = {}
    for pss_el in root.findall("PlateSections"):
        for ps_el in pss_el.findall("PlateSection"):
            ps_name = ps_el.attrib['Name']
            wl_els = ps_el.find("Wavelengths").findall("Wavelength")
            if len(wl_els) != 1:
                raise Exception(
                    "This implementation doesn't support "
                    "plate sections with multiple wavelengths")
            wl_el = wl_els[0]
            section_data = {}
            for well_el in wl_el.find("Wells").findall("Well"):
                well_name = well_el.attrib['Name']
                lum_data = [
                    int(x) for x in well_el.find("RawData").text.split()]
                time_data = [
                    float(x) for x in well_el.find("TimeData").text.split()]
                section_data[well_name] = [
                    LumDataPt(
                        ps_name=ps_name,
                        well_name=well_name,
                        sample_idx=i,
                        time=t,
                        lum=y
                        )
                    for (i, (t, y)) in enumerate(zip(time_data, lum_data))
                    ]
            data_by_section[ps_name] = section_data
    return data_by_section
