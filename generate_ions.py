#!/usr/bin/env python

import abc
import collections
import enum

Ion = collections.namedtuple("Ion", ["mass", "label", "pos"])

FIXED_MASSES = {
    "tag": 304.20536,
    "H2O": 18.006067,
    "H": 1.0073,
    "CO": 28.0101,
    "NH3": 16.998767,
    "cys": 57.021464,
    "N": 14.003074,
    "CO2": 43.989830
}

class IonType(enum.Enum):
    precursor = enum.auto()
    b = enum.auto()
    y = enum.auto()
    a = enum.auto()
    c = enum.auto()
    z = enum.auto()
    
    
class IonGenerator(metaclass=abc.ABCMeta):
    @staticmethod
    def factory(ion_type):
        if ion_type == IonType.precursor:
            return PrecursorIonGenerator()
        if ion_type == IonType.b:
            return BIonGenerator()
        if ion_type == IonType.y:
            return YIonGenerator()
        if ion_type == IonType.a:
            return AIonGenerator()
        if ion_type == IonType.c:
            return CIonGenerator()
        if ion_type == IonType.z:
            return ZIonGenerator()
        raise RuntimeError(f"Invalid IonType {ion_type} specified")
        
    def __call__(self, *args, **kwargs):
        return self.generate(*args, **kwargs)
        
    @abc.abstractmethod
    def generate(self, *args, **kwargs):
        pass
        
# BIG TODO: deduplicate ion generation functions        
class PrecursorIonGenerator(IonGenerator):
    def generate(self, mass, seq_len, charge,
                 neutral_losses=["H2O", "NH3", "CO2"], radical=False):
        """
        Generates the precursor ions for the peptide.
        
        Args:
            mass (float): The mass of the peptide.
            seq_len (int): The length of the peptide sequence.
            charge (int): The charge state of the peptide.
            
        Returns:
            The list of precursor Ions.

        """
        ions = []
        for cs in range(charge):
            charge_symbol = f"{'•' if radical else ''}{cs + 1 if cs > 0 else ''}+"
            ions.append(Ion(mass + FIXED_MASSES["H"],
                         f"[M+H][{charge_symbol}]", seq_len))
                         
            if radical:
                ions.append(Ion(mass, f"M[{charge_symbol}]", seq_len))
                
            ions.extend([
                Ion(mass - FIXED_MASSES[nl], f"[M-{nl}][{charge_symbol}]", seq_len)
                for nl in neutral_losses])
                
        return ions
        
        
class BIonGenerator(IonGenerator):
    def generate(self, b_masses, charge, neutral_losses=["H2O", "NH3", "CO"],
                 radical=False):
        """
        Generates the b-type ions for each of the given peptide sequence masses.
        
        Args:
            b_masses (list): The b-ion masses ordered by b-ion number.
            charge (int): The charge state of the peptide.
            
        Returns:
            The list of b-type Ions.

        """
        ions = []
        for idx, bm in enumerate(b_masses):
            bmh = bm + FIXED_MASSES["H"]

            ions.append(Ion(bmh, f"b{idx + 1}[+]", idx + 1))
            
            if radical:
                ions.append(Ion(bm, f"[b{idx + 1}-H][•+]", idx + 1))
                
            ions.extend([
                Ion(bmh - FIXED_MASSES[nl], f"[b{idx + 1}-{nl}][+]", idx + 1)
                for nl in neutral_losses])
                
        all_ions = ions
        for cs in range(1, charge):
            all_ions.extend(_charge_ions(ions, cs + 1))
        
        return all_ions
        
        
class YIonGenerator(IonGenerator):
    def generate(self, y_masses, charge, neutral_losses=["NH3", "H2O"], radical=False):
        """
        Generates the y-type ions for each of the given peptide sequence masses.
        
        Args:
            y_masses (list): The y-ion masses ordered by y-ion number.
            charge (int): The charge state of the peptide.
            
        Returns:
            The list of y Ion namedtuples.

        """
        ions = []
        for idx, ym in enumerate(y_masses[:-1]):
            ymh = ym + FIXED_MASSES["H"]
            
            ions.append(Ion(ymh, f"y{idx + 1}[+]", idx + 1))
            
            if radical:
                # TODO
                pass
                
            ions.extend([
                Ion(ymh - FIXED_MASSES[nl], f"[y{idx + 1}-{nl}][+]", idx + 1)
                for nl in neutral_losses])
        
        all_ions = ions
        for cs in range(1, charge):
            all_ions.extend(_charge_ions(ions, cs + 1))
        
        return all_ions
        
        
class AIonGenerator(IonGenerator):
    def generate(self, b_masses, charge, neutral_losses=[], radical=False):
        """
        """
        ions = []
        for idx, bm in enumerate(b_masses):
            am = bm + FIXED_MASSES["H"] - FIXED_MASSES["CO"]
            
            ions.append(Ion(am, f"a{idx + 1}[+]", idx + 1))
            
            if radical:
                ions.extend([
                    Ion(am - FIXED_MASSES["H"], f"[a{idx + 1}-H][•+]", idx + 1),
                    Ion(am + FIXED_MASSES["H"], f"[a{idx + 1}+H][•+]", idx + 1)])
                    
            ions.extend([
                Ion(am - FIXED_MASSES[nl], f"[a{idx+1}-{nl}][+]", idx + 1)
                for nl in neutral_losses])
        
        all_ions = ions
        for cs in range(1, charge):
            all_ions.extend(_charge_ions(ions, cs + 1))
        
        return all_ions
        
        
class CIonGenerator(IonGenerator):
    def generate(self, b_masses, charge, neutral_losses=[], radical=False):
        """
        """
        ions = []
        for idx, bm in enumerate(b_masses):
            cm = bm + FIXED_MASSES["N"] + 3 * FIXED_MASSES["H"]
            
            ions.append(Ion(cm, f"c{idx + 1}[+]", idx + 1))
            
            if radical:
                ions.append(
                    Ion(cm + 2 * FIXED_MASSES["H"], f"[c{idx + 1}+2H][•+]", idx + 1))
                    
            ions.extend([
                Ion(am - FIXED_MASSES[nl], f"[c{idx+1}-{nl}][+]", idx + 1)
                for nl in neutral_losses])
        
        all_ions = ions
        for cs in range(1, charge):
            all_ions.extend(_charge_ions(ions, cs + 1))
        
        return all_ions
        
        
class ZIonGenerator(IonGenerator):
    def generate(self, y_masses, charge, neutral_losses=[], radical=False):
        """
        """
        ions = []
        for idx, ym in enumerate(y_masses[:-1]):
            zm = ym - FIXED_MASSES["N"] - 3 * FIXED_MASSES["H"]
            
            ions.append(Ion(zm, f"z{idx + 1}[+]", idx + 1))
            
            if radical:
                ions.append(
                    Ion(zm - FIXED_MASSES["H"], f"[z{idx + 1}-H][•+]", idx + 1))
                
            ions.extend([
                Ion(zm - FIXED_MASSES[nl], f"[z{idx + 1}-{nl}][+]", idx + 1)
                for nl in neutral_losses])
        
        all_ions = ions
        for cs in range(1, charge):
            all_ions.extend(_charge_ions(ions, cs + 1))
        
        return all_ions
   

def _charge_ions(ions, charge):
    """
    Generates the multiply charged ions from the singly charged ions provided.
    
    Args:
        ions (list): A list of Ions to charge.
        charge (int): The target charge state of the Ions.
        
    Returns:
        A list of multiply charged Ions.

    """
    return [Ion((ion.mass + (charge - 1) * FIXED_MASSES["H"]) / float(charge),
                ion.label.replace("+", f"{charge}+"), ion.pos)
            for ion in ions]