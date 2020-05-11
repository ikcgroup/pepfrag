import unittest

from pepfrag import IonType, MassType, ModSite, Peptide


class TestPeptideMass(unittest.TestCase):
    """
    Tests for the pepfrag.Peptide class, focusing on mass calculation.

    """
    def test_peptide_mass_mono_no_mods(self):
        """Tests that the calculated mono mass is correct."""
        peptide = Peptide('AAA', 2, [])
        self.assertAlmostEqual(231.12, peptide.mass, 2)

    def test_peptide_mass_avg_no_mods(self):
        """Tests that the calculated average mass is correct."""
        peptide = Peptide(
            'AAA',
            2,
            [],
            mass_type=MassType.avg
        )
        self.assertEqual(231.25, peptide.mass)

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

    def test_peptide_mass_mono_cterm_mod(self):
        """
        Tests that the calculated mono mass (with a C-terminal modification)
        is correct.
        """
        peptide = Peptide('AAA', 2, [ModSite(21.981943, 'cterm', 'Cation:Na')])
        self.assertAlmostEqual(253.10, peptide.mass, 2)

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


class TestPeptideFragmentation(unittest.TestCase):
    """
    Tests for the pepfrag.Peptide class, focusing on fragmentation.

    """
    def test_basic(self):
        peptide = Peptide('AAA', 2, [])
        peptide.fragment()
        self.assertIsNotNone(peptide.fragment_ions)
        peptide.clean_fragment_ions()
        self.assertIsNone(peptide.fragment_ions)

    def test_fragment_invalidation(self):
        peptide = Peptide('AAA', 2, [])
        peptide.fragment()
        self.assertIsNotNone(peptide.fragment_ions)
        peptide.seq = 'AA'
        self.assertIsNone(peptide.fragment_ions)

    def test_custom_ion_types(self):
        peptide = Peptide('AAA', 2, [])
        peptide.fragment(ion_types={
            IonType.b.value: []
        })

        self.assertIsNotNone(peptide.fragment_ions)
        self.assertTrue(
            all(ion.startswith('b') for _, ion, _ in peptide.fragment_ions)
        )

    def test_invalid_ion_types(self):
        peptide = Peptide('AAA', 2, [])
        with self.assertRaises(RuntimeError):
            # noinspection PyTypeChecker
            peptide.fragment(ion_types={
                IonType.b: []
            })
