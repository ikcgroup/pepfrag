import unittest
from typing import Dict, List, Tuple

import numpy as np

from pepfrag.pepfrag import (
    IonType, MassType, ModSite, Peptide, _reformat_ion_types
)


def ions_to_dict(ions: List[Tuple[float, str, int]]) -> Dict[str, float]:
    return {ion[1]: ion[0] for ion in ions}


class TestPeptideMass(unittest.TestCase):
    """
    Tests for the pepfrag.Peptide class, focusing on mass calculation.

    """
    def test_peptide_mass_mono_no_mods(self):
        """Tests that the calculated mono mass is correct."""
        peptide = Peptide('AAA', 2, [])
        self.assertAlmostEqual(231.12, peptide.mass, 2)
        self.assertAlmostEqual(116.57, peptide.mz, 2)

    def test_peptide_mass_avg_no_mods(self):
        """Tests that the calculated average mass is correct."""
        peptide = Peptide(
            'AAA',
            2,
            [],
            mass_type=MassType.avg
        )
        self.assertAlmostEqual(231.24, peptide.mass, 2)
        self.assertAlmostEqual(116.63, peptide.mz, 2)

    def test_peptide_mass_mono_nterm_mod(self):
        """
        Tests that the calculated mono mass (with an N-terminal modification)
        is correct.
        """
        peptide = Peptide(
            'AAA',
            2,
            [ModSite(304.20536, 'nterm', 'iTRAQ8plex')]
        )
        self.assertAlmostEqual(535.33, peptide.mass, 2)
        self.assertAlmostEqual(268.67, peptide.mz, 2)

    def test_peptide_mass_mono_cterm_mod(self):
        """
        Tests that the calculated mono mass (with a C-terminal modification)
        is correct.
        """
        peptide = Peptide('AAA', 2, [ModSite(21.981943, 'cterm', 'Cation:Na')])
        self.assertAlmostEqual(253.10, peptide.mass, 2)
        self.assertAlmostEqual(127.56, peptide.mz, 2)

    def test_peptide_multiple_mods(self):
        """
        Tests that the calculated mono mass is correct for a peptide with many
        modifications.

        """
        peptide = Peptide(
            'AYHGMLPWK',
            3,
            [
                ModSite(304.20536, 'nterm', 'iTRAQ8plex'),
                ModSite(44.985078, 2, 'Nitro'),
                ModSite(15.994915, 5, 'Oxidation'),
                ModSite(15.994915, 7, 'Oxidation'),
                ModSite(31.989829, 8, 'Dioxidation'),
                ModSite(21.981943, 'cterm', 'Cation:Na')
            ]
        )
        self.assertAlmostEqual(1536.70, peptide.mass, 2)
        self.assertAlmostEqual(513.24, peptide.mz, 2)

    def test_invalid_residue(self):
        """
        Tests that invalid residues result in a KeyError.

        """
        peptide = Peptide('AUA', 2, [])
        with self.assertRaisesRegex(KeyError, r'Invalid residue detected: U'):
            mass = peptide.mass

    def test_mod_site_out_of_range(self):
        """
        Tests behaviour when a modification has a site beyond the peptide's
        length.

        """
        peptide = Peptide('ALPK', 2, [ModSite(1000., 100, 'TestMod')])
        self.assertAlmostEqual(427.28, peptide.mass, 2)


class TestReformatIonTypeDictionary(unittest.TestCase):
    """
    Tests for the reformat_ion_types function.

    """
    def test(self):
        ion_types = {
            IonType.b: ['NH3', ('testmod', 19.)],
            IonType.precursor: [('testmod2', 101.), 'H2O']
        }

        expected_ion_types = {
            IonType.b.value: [('NH3', 17.02654910112), ('testmod', 19.)],
            IonType.precursor.value: [('testmod2', 101.),
                                      ('H2O', 18.01056468403)]
        }

        self.assertEqual(expected_ion_types, _reformat_ion_types(ion_types))


class TestPeptideFragmentation(unittest.TestCase):
    """
    Tests for the pepfrag.Peptide class, focusing on fragmentation.

    """
    def _test_ions(self, expected: Dict[str, float], actual: Dict[str, float]):
        for ion, mass in expected.items():
            self.assertAlmostEqual(mass, actual[ion], places=2)
        self.assertEqual(set(), set(actual.keys()) - set(expected.keys()))

    def test_basic(self):
        peptide = Peptide('AAA', 2, [])
        ions = peptide.fragment()
        self.assertIsNotNone(ions)

    def test_basic_numpy(self):
        peptide = Peptide('AAA', np.int32(2), [])
        ions = peptide.fragment()
        self.assertIsNotNone(ions)

    def test_basic_mod(self):
        peptide = Peptide('AAYK', 2, [ModSite(1., 2, 'testmod')])
        ions = peptide.fragment()
        self.assertIsNotNone(ions)

    def test_basic_mod_numpy(self):
        peptide = Peptide('AAYK', 2, [ModSite(1., np.int32(2), 'testmod')])
        with self.assertRaisesRegex(RuntimeError, r'Modification site was not an integer or a string'):
            peptide.fragment()

    def test_custom_ion_types(self):
        peptide = Peptide('AAA', 2, [])
        ions = peptide.fragment(ion_types={
            IonType.b: []
        })

        self.assertIsNotNone(ions)
        self.assertTrue(
            all(ion.startswith('b') for _, ion, _ in ions)
        )

    def test_all_ion_types(self):
        peptide = Peptide('AAAMLPK', 2, [])
        ions = peptide.fragment(ion_types={
            IonType.b: [],
            IonType.y: [],
            IonType.a: [],
            IonType.c: [],
            IonType.z: [],
            IonType.precursor: [],
            IonType.imm: []
        })

        self.assertIsNotNone(ions)

    def test_a_ions(self):
        peptide = Peptide('AFCWK', 1, [])
        ions = peptide.fragment(ion_types={IonType.a: []})
        expected = {
            'a1[+]': 44.0495,
            'a2[+]': 191.1179,
            'a3[+]': 294.1271,
            'a4[+]': 480.2064,
        }
        self._test_ions(expected, ions_to_dict(ions))

    def test_b_ions(self):
        peptide = Peptide('AFCWK', 1, [])
        ions = peptide.fragment(ion_types={IonType.b: []})
        expected = {
            'b1[+]': 72.0444,
            'b2[+]': 219.1128,
            'b3[+]': 322.1220,
            'b4[+]': 508.2013,
        }
        self._test_ions(expected, ions_to_dict(ions))

    def test_c_ions(self):
        peptide = Peptide('AFCWK', 1, [])
        ions = peptide.fragment(ion_types={IonType.c: []})
        expected = {
            'c1[+]': 89.0709,
            'c2[+]': 236.1394,
            'c3[+]': 339.1485,
            'c4[+]': 525.2279,
        }
        self._test_ions(expected, ions_to_dict(ions))

    def test_x_ions(self):
        peptide = Peptide('AFCWK', 1, [])
        ions = peptide.fragment(ion_types={IonType.x: []})
        expected = {
            'x1[+]': 173.0921,
            'x2[+]': 359.1714,
            'x3[+]': 462.1806,
            'x4[+]': 609.2490,
        }
        self._test_ions(expected, ions_to_dict(ions))

    def test_y_ions(self):
        peptide = Peptide('AFCWK', 1, [])
        ions = peptide.fragment(ion_types={IonType.y: []})
        expected = {
            'y1[+]': 147.1128,
            'y2[+]': 333.1921,
            'y3[+]': 436.2013,
            'y4[+]': 583.2697,
        }
        self._test_ions(expected, ions_to_dict(ions))

    def test_z_ions(self):
        peptide = Peptide('AFCWK', 1, [])
        ions = peptide.fragment(ion_types={IonType.z: []})
        expected = {
            'z1[+]': 131.0941,
            'z2[+]': 317.1734,
            'z3[+]': 420.1826,
            'z4[+]': 567.2510,
        }
        self._test_ions(expected, ions_to_dict(ions))

    def test_invalid_residue(self):
        peptide = Peptide('AUA', 2, [])
        with self.assertRaisesRegex(KeyError, r'Invalid residue detected: U'):
            peptide.fragment()

    def test_neutral_losses(self):
        peptide = Peptide('AAAK', 2, [])
        ions = peptide.fragment(ion_types={
            IonType.b: ['NH3', ('testLoss', 9.)],
            IonType.precursor: ['CO2', ('testLoss2', 13.)],
            IonType.imm: [('testLoss3', 2.04), 'H2O']
        })

        expected_ions = {
            (44.049475632459, 'imm(A)', 0),
            (42.009475632459, '[imm1-testLoss3][+]', 1),
            (26.038910948429, '[imm1-H2O][+]', 1),
            (42.009475632459, '[imm2-testLoss3][+]', 2),
            (26.038910948429, '[imm2-H2O][+]', 2),
            (42.009475632459, '[imm3-testLoss3][+]', 3),
            (26.038910948429, '[imm3-H2O][+]', 3),
            (101.107324862499, 'imm(K)', 0),
            (72.044390252029, 'b1[+]', 1),
            (55.01784115090901, '[b1-NH3][+]', 1),
            (63.044390252029004, '[b1-testLoss][+]', 1),
            (143.081504037179, 'b2[+]', 2),
            (126.054954936059, '[b2-NH3][+]', 2),
            (134.081504037179, '[b2-testLoss][+]', 2),
            (214.118617822329, 'b3[+]', 3),
            (197.09206872120902, '[b3-NH3][+]', 3),
            (205.118617822329, '[b3-testLoss][+]', 3),
            (107.56294714460401, 'b3[2+]', 3),
            (99.04967259404401, '[b3-NH3][2+]', 3),
            (103.06294714460401, '[b3-testLoss][2+]', 3),
            (99.06732486249899, '[imm4-testLoss3][+]', 4),
            (83.096760178469, '[imm4-H2O][+]', 4),
            (21.508376049669, '[imm3-testLoss3][2+]', 3),
            (13.523093707653999, '[imm3-H2O][2+]', 3),
            (360.22414552154896, '[M+H][+]', 4),
            (316.234315521549, '[M-CO2][+]', 4),
            (347.22414552154896, '[M-testLoss2][+]', 4),
            (180.615710994214, '[M+H][2+]', 4),
            (158.620795994214, '[M-CO2][2+]', 4),
            (174.115710994214, '[M-testLoss2][2+]', 4),
            (50.037300664689, '[imm4-testLoss3][2+]', 4),
            (42.052018322674, '[imm4-H2O][2+]', 4)
        }

        self.assertEqual(expected_ions, set(ions))

    def test_unknown_string_neutral_loss(self):
        peptide = Peptide('AAAK', 2, [])
        with self.assertRaises(KeyError):
            peptide.fragment(ion_types={
                IonType.b: ['test']
            })

    def test_radical(self):
        peptide = Peptide('AAA', 2, [], radical=True)
        ions = peptide.fragment(ion_types={
            IonType.b: [],
            IonType.y: [],
            IonType.a: [],
            IonType.c: [],
            IonType.z: [],
            IonType.precursor: [],
            IonType.imm: []
        })

        self.assertIsNotNone(ions)


class TestPeptides(unittest.TestCase):
    def test_peptides_equal(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        peptide2 = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        self.assertEqual(peptide, peptide2)

    def test_peptides_not_equal_seq(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        peptide2 = Peptide(
            'ATSYPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        self.assertNotEqual(peptide, peptide2)

    def test_peptides_not_equal_charge(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        peptide2 = Peptide(
            'ATSMPLK', 3,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        self.assertNotEqual(peptide, peptide2)

    def test_peptides_not_equal_mods(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        peptide2 = Peptide(
            'ATSMPLK', 3,
            [ModSite(23.01, 'nterm', 'testmod')]
        )
        self.assertNotEqual(peptide, peptide2)

    def test_peptides_not_equal_non_instance(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        self.assertNotEqual(peptide, ('ATSMPLK', 2))

    def test_peptide_str(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        expected_str = \
            "<Peptide {'seq': 'ATSMPLK', 'charge': 2, 'mods': " \
            "[ModSite(mass=23.01, site='nterm', mod='testmod'), " \
            "ModSite(mass=19.24, site=2, mod='testmod2')], 'mass_type': " \
            "<MassType.mono: 0>, 'radical': False}>"
        self.assertEqual(expected_str, str(peptide))

    def test_peptide_repr(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        expected_repr = \
            "<Peptide {'seq': 'ATSMPLK', 'charge': 2, 'mods': " \
            "[ModSite(mass=23.01, site='nterm', mod='testmod'), "\
            "ModSite(mass=19.24, site=2, mod='testmod2')], 'mass_type': "\
            "<MassType.mono: 0>, 'radical': False}>"
        self.assertEqual(expected_repr, repr(peptide))

    def test_peptides_equal_hash(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        peptide2 = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        self.assertEqual(hash(peptide), hash(peptide2))
