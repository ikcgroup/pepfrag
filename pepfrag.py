#! /usr/bin/env python3
"""
This module is used to fragment a peptide, including its modifications, to
generate theoretical ions useful for annotating mass spectra (MSn).

"""
import collections
import enum

from ion_generators import FIXED_MASSES, IonType, IonGenerator


PeptideMass = collections.namedtuple("PeptideMass", ["nterm", "seq", "cterm"])

    
# TODO: deduplicate with constants.py
class MassType(enum.Enum):
    mono = enum.auto()
    avg = enum.auto()
    
Mass = collections.namedtuple('Mass', ['mono', 'avg'])

AA_MASSES = {
    'G': Mass(57.02146, 57.052),
    'A': Mass(71.03711, 71.078),
    'S': Mass(87.03203, 87.078),
    'P': Mass(97.05276, 97.117),
    'V': Mass(99.06841, 99.133),
    'T': Mass(101.04768, 101.105),
    'C': Mass(103.00918, 103.144),
    'I': Mass(113.08406, 113.160),
    'L': Mass(113.08406, 113.160),
    'N': Mass(114.04292, 114.104),
    'D': Mass(115.02693, 115.089),
    'Q': Mass(128.05857, 128.131),
    'K': Mass(128.09495, 128.174),
    'E': Mass(129.04258, 129.116),
    'M': Mass(131.04048, 131.198),
    'H': Mass(137.05891, 137.142),
    'F': Mass(147.06841, 147.177),
    'R': Mass(156.10110, 156.188),
    'Y': Mass(163.06332, 163.170),
    'W': Mass(186.07931, 186.213)
}
    
    
DEFAULT_IONS = {
    IonType.precursor: {},
    IonType.imm: {},
    IonType.b: {},
    IonType.y: {},
    IonType.a: {},
    IonType.c: {},
    IonType.z: {}
}


class Peptide():
    """
    A class to represent a peptide, including its charge state and any
    modifications, including PTMs and quantitative tags. The class should be
    used to fragment the peptides for mass spectrum annotation.

    """
    def __init__(self, sequence, charge, modifications, mass_type=MassType.mono,
                 radical=False):
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
        
    @property
    def peptide_mass(self):
        return self.calculate_mass()
        
    @property
    def mass(self):
        pep_mass = self.peptide_mass
        mass = sum(pep_mass.seq) + FIXED_MASSES["H2O"]
        if pep_mass.nterm is not None:
            mass += pep_mass.nterm
        return mass
        
    def __repr__(self):
        """
        Constructs the official representation of the Peptide object.
        
        Returns:
            string: Official representation of the Peptide object.
            
        """
        return f"<{self.__class__.__name__} {self.__dict__}>"
        
    def calculate_mass(self):
        """
        Calculates the theoretical mass of the peptide, including
        any modifications.
            
        Returns:
            PeptideMass

        """
        nterm_mass, cterm_mass = None, None
        # Initialize added modification masses for each sequence position
        site_mod_masses = [0.] * len(self.seq)
        
        for mod in self.mods:
            if isinstance(mod.site, int):
                site_mod_masses[mod.site - 1] += mod.mass
            elif mod.site == "cterm":
                cterm_mass = mod.mass
            else:
                nterm_mass = mod.mass
        
        seq_masses = [getattr(AA_MASSES[aa], self.mass_type.name) +
                      site_mod_masses[ii] for ii, aa in enumerate(self.seq)]
                      
        return PeptideMass(nterm_mass, seq_masses, cterm_mass)

    def fragment(self, ion_types=DEFAULT_IONS):
        """
        Fragments the peptide to generate the ion types specified.
        
        Args:
            ion_types (dict):
            
        """
        mass = self.mass
            
        bm, ym = self._ion_masses()
        
        ions = []
        
        for ion_type in ion_types:
            generator = IonGenerator.factory(ion_type)
            if ion_type == IonType.precursor:
                ions.extend(generator(mass, self.charge, len(self.seq),
                                      self.mods, radical=self.radical,
                                      **ion_types[ion_type]))
            elif ion_type == IonType.imm:
                ions.extend(generator(self.mass.seq, self.charge, self.seq,
                                      self.mods, **ion_types[ion_type]))
            else:
                masses = None
                if ion_type in [IonType.b, IonType.a, IonType.c]:
                    masses = bm
                elif ion_type in [IonType.y, IonType.z]:
                    masses = ym
                    
                if masses is None:
                    raise RuntimeError(
                        f"Invalid IonType {ion_type} specified")
                        
                ions.extend(
                    generator(masses, self.charge, radical=self.radical,
                              **ion_types[ion_type]))
            
        return ions
        
    def _ion_masses(self):
        """
        Generates the theoretical ion masses of the peptide based on its
        sequence and modification masses.
            
        Returns:
            Tuple of two lists: the b and y ion masses.

        """
        seq_len = len(self.seq)
        
        pep_mass = self.peptide_mass

        # The base mass for y-type ions
        y_base = (FIXED_MASSES["H2O"] if self.pep_mass.cterm is None
                  else self.pep_mass.cterm)
        rev_seq_masses = self.pep_mass.seq[::-1]
        y_ions = [y_base + sum(rev_seq_masses[:ii])
                  for ii in range(1, seq_len + 1)]
                  
        # The base mass for b-type ions
        b_base = 0. if self.pep_mass.nterm is None else self.pep_mass.nterm
        b_ions = [b_base + sum(self.pep_mass.seq[:ii])
                  for ii in range(1, seq_len + 1)]
        
        return b_ions, y_ions