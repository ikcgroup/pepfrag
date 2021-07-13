#! /usr/bin/env python3
"""
This module is used to fragment a peptide, including its modifications, to
generate theoretical ions useful for annotating mass spectra (MSn).

"""
from __future__ import annotations

import dataclasses
import enum
from typing import Dict, List, Optional, Sequence, Tuple, Union

from cpepfrag import calculate_mass, generate_ions

from .constants import AA_MASSES, FIXED_MASSES, MassType


Ion = Tuple[float, str, int]


@dataclasses.dataclass(frozen=True)
class ModSite:
    """
    Class representing an instance of `mod_name` at position `site`.

    Args:
        mass: Mass of the modification.
        site: Position of the modification. Integer for sequence position,
              'nterm' for N-terminus or 'cterm' for C-terminus.
        mod: Name of the modification.

    """
    mass: float
    site: Union[int, str]
    mod: str


class IonType(enum.Enum):
    """
    Enumeration of possible fragment ion types.

    """
    precursor = 1  #: Precursor ions
    imm = 2  #: Immonium ions
    b = 3  #: b-type ions
    y = 4  #: y-type ions
    a = 5  #: a-type ions
    c = 6  #: c-type ions
    z = 7  #: z-type ions
    x = 8  #: x-type ions


IonTypesDict = Dict[IonType, List[Union[str, Tuple[str, float]]]]
# The dictionary format to be passed to the C++ extension
CIonTypesDict = Dict[int, List[Tuple[str, float]]]


DEFAULT_IONS: IonTypesDict = {
    IonType.precursor: ["H2O", "NH3", "CO2"],
    IonType.imm: [],
    IonType.b: ["H2O", "NH3", "CO"],
    IonType.y: ["NH3", "H2O"],
    IonType.a: [],
    IonType.c: [],
    IonType.z: [],
    IonType.x: []
}


AA_TYPE_MASSES = {
    (mass_type, aa): getattr(masses, mass_type.name)
    for aa, masses in AA_MASSES.items()
    for mass_type in [MassType.mono, MassType.avg]
}


def _reformat_ion_types(ion_types: IonTypesDict) -> CIonTypesDict:
    """
    Reformats the `ion_types` dictionary to convert string neutral losses to
    tuples and use the integer value of the IonType enumeration.

    Args:
        ion_types: The selected ion type dictionary.

    Returns:
        Reformatted ion type dictionary.

    """
    new_ion_types = {}
    for ion_type, losses in ion_types.items():
        for i, loss in enumerate(losses):
            if isinstance(loss, str):
                try:
                    losses[i] = (loss, FIXED_MASSES[loss])
                except KeyError:
                    raise KeyError(
                        f'Unknown neutral loss: {loss}. Consider using the mass'
                        'directly.'
                    )
        new_ion_types[ion_type.value] = losses
    return new_ion_types


class UnknownModificationSite(Exception):
    """
    An exception to represent the detection of an unknown/uninterpretable
    modification site.

    """


class Peptide:
    """
    A class to represent a peptide, including its charge state and any
    modifications, including PTMs and quantitative tags. The class should be
    used to fragment the peptides for mass spectrum annotation.

    Attributes:
        mass_type: Type of masses used in calculations (see :class:`MassType`).
        radical: Flag indicating whether the peptide is a radical peptide.

    """

    __slots__ = ("seq", "charge", "mods", "mass_type", "radical",)

    def __init__(
            self,
            sequence: str,
            charge: int,
            modifications: Sequence[ModSite],
            mass_type: MassType = MassType.mono,
            radical: bool = False
    ):
        """
        Initializes the Peptide object.

        Args:
            sequence: The peptide sequence (single character format).
            charge: The charge state of the peptide.
            modifications: The modifications applied to the peptide.
            mass_type: The type of masses used in calculations
                       (see :class:`MassType`).
            radical: Flag indicating whether the peptide is a radical peptide.
                     This flag influences the ion candidates generated during
                     fragmentation.

        """
        self.seq = sequence
        self.charge = charge
        self.mods = modifications
        self.mass_type: MassType = mass_type
        self.radical: bool = radical

    @property
    def peptide_mass(self) -> List[float]:
        """
        The mass of the peptide along the sequence, with each position
        calculated separately.

        Note:
            In the returned list, index 0 is the N-terminus mass, while index -1
            is the C-terminus mass.

        """
        return self.calculate_mass()

    @property
    def mass(self) -> float:
        """
        Total mass of the peptide, including modifications.

        """
        pep_mass = self.peptide_mass
        mass = sum(pep_mass) + FIXED_MASSES["H2O"]
        return mass

    @property
    def mz(self) -> float:
        """
        Calculates the mass-to-charge ratio of the peptide.

        """
        return (self.mass / self.charge) + FIXED_MASSES["H"]

    def __repr__(self) -> str:
        """
        Constructs the official representation of the Peptide object.

        Returns:
            Official representation of the Peptide object.

        """
        out = {s: getattr(self, s) for s in self.__class__.__slots__}
        return f"<{self.__class__.__name__} {out}>"

    def __str__(self) -> str:
        """
        Constructs the string representation of the Peptide object.

        Returns:
            String representation.

        """
        out = {
            "seq": self.seq,
            "charge": self.charge,
            "mods": self.mods,
            "mass_type": self.mass_type,
            "radical": self.radical
        }
        return f"<{self.__class__.__name__} {out}>"

    def __hash__(self):
        """
        Implements the hash function for the Peptide object.

        """
        return hash((self.seq, self.charge, tuple(self.mods)))

    def __eq__(self, other: object) -> bool:
        """
        Implements the equality test for the Peptide object.

        """
        if not isinstance(other, Peptide):
            return NotImplemented
        return (self.seq, self.charge, self.mods) == \
               (other.seq, other.charge, other.mods)

    def calculate_mass(self) -> List[float]:
        """
        Calculates the theoretical mass of the peptide along the sequence,
        including any modifications.

        Returns:
            Masses along the peptide sequence. Index 0 is the N-terminus mass,
            while index -1 is the C-terminus mass.

        """
        return calculate_mass(
            self.seq,
            self.mods,
            self.mass_type.value
        )

    def fragment(
            self,
            ion_types: Optional[IonTypesDict] = None
    ) -> List[Ion]:
        """
        Fragments the peptide to generate the ion types specified.

        Args:
            ion_types: Dictionary of :class:`IonType` s to list of configured
                       neutral losses. Only fragments for :class:`IonType` s
                       specified here will be generated.

        Returns:
            List of generated ions, as tuples of `(fragment mass, ion label,
            sequence position)`.

        """
        if ion_types is None:
            ion_types = DEFAULT_IONS

        ion_types = _reformat_ion_types(ion_types)

        return self._fragment(ion_types)

    def _fragment(
            self,
            ion_types: CIonTypesDict
    ) -> List[Ion]:
        """
        Fragments the peptide to generate the ion types specified.

        Args:
            ion_types (dict):

        """
        b_masses, y_masses = self._ion_masses()
        return generate_ions(
            ion_types,
            self.mass,
            self.peptide_mass[1:-1],
            b_masses,
            y_masses,
            self.charge,
            self.radical,
            self.seq
        )

    def _ion_masses(self) -> Tuple[List[float], List[float]]:
        """
        Generates the theoretical ion masses of the peptide based on its
        sequence and modification masses.

        Returns:
            Tuple of two lists: the b and y ion masses.

        """
        seq_len = len(self.seq)

        pep_mass = self.peptide_mass

        y_base = FIXED_MASSES["H2O"] + pep_mass[-1]
        b_base = pep_mass[0]
        # Reverse the sequence portion of the list (i.e. pep_mass[1:-1])
        rev_seq_masses = pep_mass[-2:0:-1]
        y_ions = []
        b_ions = []
        for ii in range(seq_len):
            y_base += rev_seq_masses[ii]
            y_ions.append(y_base)
            b_base += pep_mass[ii + 1]
            b_ions.append(b_base)

        return b_ions, y_ions
