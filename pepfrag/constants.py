#! /usr/bin/env python3
"""
A module to define some constants used throughout peptide processing.

"""
import dataclasses
import enum
from typing import Dict


class MassType(enum.Enum):
    """
    An enumeration representing the possible mass types.

    Note that the values of these enumerations correspond to their index in
    `Mass`, and similarly in the C++ code underneath methods such as
    calculate_mass.

    """
    mono = 0  #: Monoisotopic mass
    avg = 1  #: Average mass


@dataclasses.dataclass
class Mass:
    """
    Represents a mass pair of monoisotopic and average masses.

    Args:
        mono: Monoisotopic mass.
        avg: Average mass.

    """
    mono: float
    avg: float


AA_MASSES: Dict[str, Mass] = {
    'G': Mass(57.02146372069, 57.051402191402),
    'A': Mass(71.03711378515, 71.078019596249),
    'S': Mass(87.03202840472, 87.077424520567),
    'P': Mass(97.05276384961, 97.115372897831),
    'V': Mass(99.06841391407, 99.131254405943),
    'T': Mass(101.04767846918, 101.104041925414),
    'C': Mass(103.00918495955, 103.142807002376),
    'I': Mass(113.08406397853, 113.157871810790),
    'L': Mass(113.08406397853, 113.157871810790),
    'N': Mass(114.04292744138, 114.102804382804),
    'D': Mass(115.02694302429, 115.087565341620),
    'Q': Mass(128.05857750584, 128.129421787651),
    'K': Mass(128.09496301519, 128.172515776292),
    'E': Mass(129.04259308875, 129.114182746467),
    'M': Mass(131.04048508847, 131.19604181207),
    'H': Mass(137.05891185847, 137.139515217458),
    'F': Mass(147.06841391407, 147.174197992883),
    'R': Mass(156.10111102405, 156.185922199184),
    'Y': Mass(163.06332853364, 163.173602917201),
    'W': Mass(186.07931295073, 186.210313751855)
}

FIXED_MASSES: Dict[str, float] = {
    # iTRAQ tag
    "tag": 304.20536,
    "H2O": 18.01056468403,
    "H": 1.007276466879,
    "CO": 27.99491461957,
    "CO2": 43.989830,
    "NH3": 17.02654910112,
    "N": 14.003074,
    # Carbamidomethylation
    "cys_c": 57.021464
}
