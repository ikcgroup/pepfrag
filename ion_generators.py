#!/usr/bin/env python
"""
A module providing classes for generating fragment ions of a peptide.

"""

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

# IonGenerator is the base class for all ion generators. It provides a factory
# method to construct the relevant generator based on the given IonType enum.
# Subclasses of IonGenerator must implement a generate method, which is called
# when the generator is called. To extend the base class' generate method,
# at least the abstract base_ions method must be overridden by subclasses.
# Further customization can be achieved by providing overrides for fix_mass,
# radical and neutral_losses.
class IonGenerator(metaclass=abc.ABCMeta):
    '''
    This class is the base class for all ion generators and provides a
    construction factory based on the IonType.

    '''
    @staticmethod
    def factory(ion_type):
        '''
        Constructs the relevant subclass of IonGenerator based on the
        ion_type.

        Args:
            ion_type (IonType): An enumeration value for one of the supported
                                types.

        Returns:
            IonGenerator (subclass)

        '''
        return ion_type.value()

    def __call__(self, *args, **kwargs):
        '''
        Makes the IonGenerator class callable. When called, the generate
        method, which must be overridden, will be called using all of the
        given arguments.

        Returns:
            List of generated Ions.

        '''
        return self.generate(*args, **kwargs)

    @abc.abstractmethod
    def generate(self, masses, charge, neutral_losses, radical):
        '''
        Generates the fragment ions for the peptide, with the type according
        to the IonGenerator used.

        This method is abstract and must be overridden by any subclasses of
        IonGenerator. However, an implementation is provided to reduce code
        duplication. In most cases, the subclass should refer its arguments
        back to this method.

        Note that this is designed in this manner to allow each IonGenerator
        subclass to define its own default arguments, e.g., for
        neutral_losses.

        Args:
            masses (list): A list of peptide ion masses (float).
            charge (int): The charge state of the peptide.
            neutral_losses (list): A list of strings indicated for which
                                   neutral losses fragment ions should be
                                   created. These must be found in
                                   FIXED_MASSES.
            radical (bool): A boolean flag indicating whether radical ions
                            should be generated.

        Returns:
            List of generated Ions.

        '''
        ions = []
        for pos, mass in enumerate(masses):
            ion_mass = self.fix_mass(mass)

            ions.extend(self.base_ions(ion_mass, pos))

            if radical:
                ions.extend(self.radical(ion_mass, pos))

            ions.extend(self.neutral_losses(ion_mass, pos, neutral_losses))

        all_ions = ions
        for _charge in range(1, charge):
            all_ions.extend(_charge_ions(ions, _charge + 1))

        return all_ions

    @abc.abstractmethod
    def base_ions(self, mass, pos):
        '''
        Generates the basic fragment ion(s) for the type. For example, b1[+]
        and y[1]+.

        This method is abstract and must be overridden by any subclasses.

        Args:
            mass (float): The mass of the ion.
            pos (int): The index in the peptide sequence of the ion.

        Returns:
            List of generated Ions.

        '''

    def fix_mass(self, mass):
        '''
        Given the b/y type ion mass, modifies the mass according to the ion
        type under consideration.

        This method should be overridden as necessary.

        Args:
            mass (float): The b/y-type ion mass.

        Returns:
            The modified mass as a float.

        '''
        return mass

    def radical(self, mass, pos):
        '''
        Generates fragment ions specific to radical peptides. This method is
        only called (using the base class generate method) if the radical flag
        is set to True.

        This method should be overridden as necessary.

        Args:
            mass (float): The mass of the ion.
            pos (int): The index in the peptide sequence of the ion.

        Returns:
            List of generated Ions.

        '''
        return []

    def neutral_losses(self, mass, pos, neutral_losses):
        '''
        Generates fragment ions deriving from neutral losses.

        This method should be overridden as necessary.

        Args:
            mass (float): The mass of the ion.
            pos (int): The index in the peptide sequence of the ion.
            neutral_losses (list): A list of strings indicated for which
                                   neutral losses fragment ions should be
                                   created. These must be found in
                                   FIXED_MASSES.

        Returns:
            List of generated Ions.

        '''
        return []


# TODO: refactor this to be similar to the other IonGenerators - will require
# the interface to be modified for extra parameters
class PrecursorIonGenerator(IonGenerator):
    """
    Generator for precursor ions.

    """
    def generate(self, mass, charge, seq_len, mods,
                 neutral_losses=["H2O", "NH3", "CO2"], itraq=False,
                 radical=False):
        """
        Generates the precursor ions for the peptide.

        Args:
            mass (float): The mass of the peptide.
            charge (int): The charge state of the peptide.
            seq_len (int): The length of the peptide sequence.
            mods (list of ModSites): The modifications applied to the peptide.

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

            if any(ms.mod == 'iTRAQ8plex' and ms.site in ["cterm", "nterm"]
                   for ms in mods):
                ions.append(Ion(mass + FIXED_MASSES["H"] - FIXED_MASSES["tag"],
                                f"M-iT8[{charge_symbol}]", seq_len))

        return ions

    def base_ions(self, mass, pos):
        return []


class ImmoniumIonGenerator(IonGenerator):
    """
    Generator for immonium ions.

    """
    def generate(self, seq_masses, charge, sequence, mods, neutral_losses=None,
                 radical=False):
        """
        Pass specific defaults back to base class.

        """
        self.sequence = sequence
        self.mod_sites = [mod.site - 1 for mod in mods
                          if isinstance(mod.site, int)]

        return super().generate(seq_masses, charge, neutral_losses, radical)

    def base_ions(self, mass, pos):
        return [Ion(mass,
                    f"imm({self.sequence[pos]}"
                    f"{'*' if pos in self.mod_sites else ''})", 0)]

    def fix_mass(self, mass):
        return mass - FIXED_MASSES["CO"] + FIXED_MASSES["H"]


class BIonGenerator(IonGenerator):
    """
    Generator for b type ions.

    """
    def generate(self, b_masses, charge, neutral_losses=None, radical=False):
        """
        Pass specific defaults back to base class.

        """
        if neutral_losses is None:
            neutral_losses = ["H2O", "NH3", "CO"]
        return super().generate(b_masses, charge, neutral_losses, radical)

    def base_ions(self, mass, pos):
        return [Ion(mass, f"b{pos + 1}[+]", pos + 1)]

    def fix_mass(self, mass):
        return mass + FIXED_MASSES["H"]

    def radical(self, mass, pos):
        return [Ion(mass, f"[b{pos + 1}-H][•+]", pos + 1)]

    def neutral_losses(self, mass, pos, neutral_losses):
        return [Ion(mass - FIXED_MASSES[nl], f"[b{pos + 1}-{nl}][+]", pos + 1)
                for nl in neutral_losses]


class YIonGenerator(IonGenerator):
    """
    Generator for y type ions.

    """
    def generate(self, y_masses, charge, neutral_losses=None, radical=False):
        """
        Pass specific defaults back to base class.

        """
        if neutral_losses is None:
            neutral_losses = ["NH3", "H2O"]
        return super().generate(y_masses, charge, neutral_losses, radical)

    def base_ions(self, mass, pos):
        return [Ion(mass, f"y{pos + 1}[+]", pos + 1)]

    def fix_mass(self, mass):
        return mass + FIXED_MASSES["H"]

    def radical(self, mass, pos):
        # TODO
        return []

    def neutral_losses(self, mass, pos, neutral_losses):
        return [Ion(mass - FIXED_MASSES[nl], f"[y{pos + 1}-{nl}][+]", pos + 1)
                for nl in neutral_losses]


class AIonGenerator(IonGenerator):
    """
    Generator for a type ions.

    """
    def generate(self, b_masses, charge, neutral_losses=None, radical=False):
        """
        Pass specific defaults back to base class.

        """
        return super().generate(b_masses, charge, neutral_losses, radical)

    def base_ions(self, mass, pos):
        return [Ion(mass, f"a{pos + 1}[+]", pos + 1)]

    def fix_mass(self, mass):
        return mass + FIXED_MASSES["H"] - FIXED_MASSES["CO"]

    def radical(self, mass, pos):
        return [Ion(mass - FIXED_MASSES["H"], f"[a{pos + 1}-H][•+]", pos + 1),
                Ion(mass + FIXED_MASSES["H"], f"[a{pos + 1}+H][•+]", pos + 1)]

    def neutral_losses(self, mass, pos, neutral_losses):
        return [Ion(mass - FIXED_MASSES[nl], f"[a{pos + 1}-{nl}][+]", pos + 1)
                for nl in neutral_losses]


class CIonGenerator(IonGenerator):
    """
    Generator for c type ions.

    """
    def generate(self, b_masses, charge, neutral_losses=None, radical=False):
        """
        Pass specific defaults back to base class.

        """
        # Before passing to the base class method, the last mass is removed
        # since the c{seq_len} ion is not a sensible annotation
        return super().generate(b_masses[:-1], charge, neutral_losses, radical)

    def base_ions(self, mass, pos):
        return [Ion(mass, f"c{pos + 1}[+]", pos + 1)]

    def fix_mass(self, mass):
        return mass + FIXED_MASSES["N"] + 3 * FIXED_MASSES["H"]

    def radical(self, mass, pos):
        return [Ion(mass + 2 * FIXED_MASSES["H"], f"[c{pos + 1}+2H][•+]",
                    pos + 1)]

    def neutral_losses(self, mass, pos, neutral_losses):
        return [Ion(mass - FIXED_MASSES[nl], f"[c{pos + 1}-{nl}][+]", pos + 1)
                for nl in neutral_losses]


class ZIonGenerator(IonGenerator):
    """
    Generator for z type ions.

    """
    def generate(self, y_masses, charge, neutral_losses=None, radical=False):
        """
        Pass specific defaults back to base class.

        """
        return super().generate(y_masses, charge, neutral_losses, radical)

    def base_ions(self, mass, pos):
        return [Ion(mass, f"z{pos + 1}[+]", pos + 1)]

    def fix_mass(self, mass):
        return mass - FIXED_MASSES["N"] - 3 * FIXED_MASSES["H"]

    def radical(self, mass, pos):
        return [Ion(mass - FIXED_MASSES["H"], f"[z{pos + 1}-H][•+]", pos + 1)]

    def neutral_losses(self, mass, pos, neutral_losses):
        return [Ion(mass - FIXED_MASSES[nl], f"[z{pos + 1}-{nl}][+]", pos + 1)
                for nl in neutral_losses]


class IonType(enum.Enum):
    """
    An enumeration mapping the ion types to their generator classes.

    """
    precursor = PrecursorIonGenerator
    imm = ImmoniumIonGenerator
    b = BIonGenerator
    y = YIonGenerator
    a = AIonGenerator
    c = CIonGenerator
    z = ZIonGenerator


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
