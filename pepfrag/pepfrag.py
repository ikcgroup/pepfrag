#! /usr/bin/env python3
"""
This module is used to fragment a peptide, including its modifications, to
generate theoretical ions useful for annotating mass spectra (MSn).

"""
from __future__ import annotations

import collections
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .constants import AA_MASSES, MassType
from .ion_generators import FIXED_MASSES, Ion, IonType, IonGenerator


PeptideMass = collections.namedtuple("PeptideMass", ["nterm", "seq", "cterm"])

ModSite = collections.namedtuple("ModSite", ["mass", "site", "mod"])


DEFAULT_IONS: Dict[IonType, Dict[str, Any]] = {
    IonType.precursor: {},
    IonType.imm: {},
    IonType.b: {},
    IonType.y: {},
    IonType.a: {},
    IonType.c: {},
    IonType.z: {}
}


class UnknownModificationSite(Exception):
    """
    An exception to represent the detection of an unknown/uninterpretable
    modification site.

    """


class Peptide():
    """
    A class to represent a peptide, including its charge state and any
    modifications, including PTMs and quantitative tags. The class should be
    used to fragment the peptides for mass spectrum annotation.

    """

    __slots__ = ("_seq", "_charge", "_mods", "mass_type", "radical",
                 "fragment_ions",)

    def __init__(self, sequence: str, charge: int,
                 modifications: Sequence[ModSite],
                 mass_type: MassType = MassType.mono,
                 radical: bool = False) -> None:
        """
        Initializes the Peptide object.

        Args:
            sequence (str): The peptide sequence (single character format).
            charge (int): The charge state of the peptide.
            modifications (list): TODO: ModSites?

        """
        self.seq = sequence
        self.charge = charge
        self.mods = modifications
        self.mass_type = mass_type
        self.radical = radical

        self.fragment_ions: Optional[List[Ion]] = None

    @property
    def seq(self) -> str:
        """
        Returns the peptide sequence.

        """
        return self._seq

    @seq.setter
    def seq(self, seq: str):
        """
        Sets the sequence and clears the fragment_ions.

        """
        self.clean_fragment_ions()
        self._seq = seq

    @property
    def charge(self) -> int:
        """
        Returns the peptide charge.

        """
        return self._charge

    @charge.setter
    def charge(self, charge: int):
        """
        Sets the charge and clears the fragment_ions.

        """
        self.clean_fragment_ions()
        self._charge = charge

    @property
    def mods(self) -> List:
        """
        Returns the list of ModSites.

        """
        return self._mods

    @mods.setter
    def mods(self, mods: List):
        """
        Sets the modifications and clears the fragment_ions.

        """
        self.clean_fragment_ions()
        self._mods = mods

    @property
    def peptide_mass(self) -> PeptideMass:
        """
        Returns the mass of the peptide, including modifications.

        """
        return self.calculate_mass()

    @property
    def mass(self) -> float:
        """
        Returns the mass of the peptide, indcluding modifications, as a float.

        """
        pep_mass = self.peptide_mass
        mass = sum(pep_mass.seq) + FIXED_MASSES["H2O"]
        if pep_mass.nterm is not None:
            mass += pep_mass.nterm
        return mass

    def __repr__(self) -> str:
        """
        Constructs the official representation of the Peptide object.

        Returns:
            string: Official representation of the Peptide object.

        """
        out = {s: getattr(self, s) for s in self.__class__.__slots__}
        return f"<{self.__class__.__name__} {out}>"

    def __str__(self) -> str:
        """
        Constructs the string representation of the Peptide object.

        Returns:
            string.

        """
        out = {
            "seq": self.seq,
            "charge": self.charge,
            "mods": self.mods,
            "mass_type": self.mass_type,
            "radical": self.radical,
            "fragment_ions": f"{len(self.fragment_ions) if self.fragment_ions is not None else 0} ions"
        }
        return f"<{self.__class__.__name__} {out}>"

    def __hash__(self):
        """
        Implements the hash function for the Peptide object.

        """
        return hash((self.seq, self.charge))

    def __eq__(self, other: object) -> bool:
        """
        Implements the equality test for the Peptide object.

        """
        if not isinstance(other, Peptide):
            raise NotImplementedError()
        return (self.seq, self.charge, self.mods) == \
               (other.seq, other.charge, other.mods)

    def clean_fragment_ions(self):
        """
        Cleans the cached fragment ions.

        """
        self.fragment_ions = None

    def calculate_mass(self) -> PeptideMass:
        """
        Calculates the theoretical mass of the peptide, including
        any modifications.

        Returns:
            PeptideMass

        Raises:
            UnknownModificationSite

        """
        nterm_mass, cterm_mass = None, None
        # Initialize added modification masses for each sequence position
        site_mod_masses = [0.] * len(self.seq)

        for mod in self.mods:
            if isinstance(mod.site, int):
                site_mod_masses[mod.site - 1] += mod.mass
            else:
                site = mod.site.lower().replace("-", "")
                if site == "cterm":
                    cterm_mass = mod.mass
                elif site == "nterm":
                    nterm_mass = mod.mass
                else:
                    raise UnknownModificationSite()

        seq_masses = [getattr(AA_MASSES[aa], self.mass_type.name) +
                      site_mod_masses[ii] for ii, aa in enumerate(self.seq)]

        return PeptideMass(nterm_mass, seq_masses, cterm_mass)

    def fragment(self,
                 ion_types: Dict[IonType, Dict[str, Any]] = DEFAULT_IONS,
                 force: bool = False) -> List[Ion]:
        """
        Fragments the peptide to generate the ion types specified.

        Args:
            ion_types (dict):

        """
        # If fragment_ions already exists or force=False, use the cached ions
        if self.fragment_ions is None or force:
            # Cache the new ions
            self.fragment_ions = self._fragment(ion_types)

        return self.fragment_ions

    def _fragment(self,
                  ion_types: Dict[IonType, Dict[str, Any]]) -> List[Ion]:
        """
        Fragments the peptide to generate the ion types specified.

        Args:
            ion_types (dict):

        """
        mass = self.mass
        pep_mass = self.peptide_mass

        b_masses, y_masses = self._ion_masses()

        ions: List[Ion] = []

        for ion_type in ion_types:
            generator: IonGenerator = IonGenerator.factory(ion_type)
            if ion_type == IonType.precursor:
                ions.extend(generator(mass, self.charge, len(self.seq),
                                      self.mods, radical=self.radical,
                                      **ion_types[ion_type]))
            elif ion_type == IonType.imm:
                ions.extend(generator(pep_mass.seq, self.charge, self.seq,
                                      self.mods, **ion_types[ion_type]))
            else:
                masses = None
                if ion_type in [IonType.b, IonType.a, IonType.c]:
                    masses = b_masses
                elif ion_type in [IonType.y, IonType.z]:
                    masses = y_masses

                if masses is None:
                    raise RuntimeError(
                        f"Invalid IonType {ion_type} specified")

                ions.extend(
                    generator(masses, self.charge, radical=self.radical,
                              **ion_types[ion_type]))

        return ions

    def _ion_masses(self) -> Tuple[List[float], List[float]]:
        """
        Generates the theoretical ion masses of the peptide based on its
        sequence and modification masses.

        Returns:
            Tuple of two lists: the b and y ion masses.

        """
        seq_len = len(self.seq)

        pep_mass = self.peptide_mass

        # The base mass for y-type ions
        y_base = (FIXED_MASSES["H2O"] if pep_mass.cterm is None
                  else pep_mass.cterm)
        rev_seq_masses = pep_mass.seq[::-1]
        y_ions = [y_base + sum(rev_seq_masses[:ii])
                  for ii in range(1, seq_len + 1)]

        # The base mass for b-type ions
        b_base = 0. if pep_mass.nterm is None else pep_mass.nterm
        b_ions = [b_base + sum(pep_mass.seq[:ii])
                  for ii in range(1, seq_len + 1)]

        return b_ions, y_ions
