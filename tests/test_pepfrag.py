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
        self.assertAlmostEqual(231.24, peptide.mass, 2)

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

    def test_invalid_residue(self):
        """
        Tests that invalid residues result in a KeyError.

        """
        peptide = Peptide('AUA', 2, [])
        with self.assertRaisesRegex(KeyError, r'Invalid residue detected: U'):
            mass = peptide.mass


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

    def test_all_ion_types(self):
        peptide = Peptide('AAAMLPK', 2, [])
        peptide.fragment(ion_types={
            IonType.b.value: [],
            IonType.y.value: [],
            IonType.a.value: [],
            IonType.c.value: [],
            IonType.z.value: [],
            IonType.precursor.value: [],
            IonType.imm.value: []
        })

        self.assertIsNotNone(peptide.fragment_ions)

    def test_invalid_ion_types(self):
        peptide = Peptide('AAA', 2, [])
        with self.assertRaises(RuntimeError):
            # noinspection PyTypeChecker
            peptide.fragment(ion_types={
                IonType.b: []
            })

    def test_invalid_residue(self):
        peptide = Peptide('AUA', 2, [])
        with self.assertRaisesRegex(KeyError, r'Invalid residue detected: U'):
            peptide.fragment()

    def test_radical(self):
        peptide = Peptide('AAA', 2, [], radical=True)
        peptide.fragment(ion_types={
            IonType.b.value: [],
            IonType.y.value: [],
            IonType.a.value: [],
            IonType.c.value: [],
            IonType.z.value: [],
            IonType.precursor.value: [],
            IonType.imm.value: []
        })

        self.assertIsNotNone(peptide.fragment_ions)


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
            "<MassType.mono: 0>, 'radical': False, 'fragment_ions': '0 ions'}>"
        self.assertEqual(expected_str, str(peptide))

    def test_peptide_repr(self):
        peptide = Peptide(
            'ATSMPLK', 2,
            [ModSite(23.01, 'nterm', 'testmod'), ModSite(19.24, 2, 'testmod2')]
        )
        expected_repr = \
            "<Peptide {'_seq': 'ATSMPLK', '_charge': 2, '_mods': " \
            "[ModSite(mass=23.01, site='nterm', mod='testmod'), "\
            "ModSite(mass=19.24, site=2, mod='testmod2')], 'mass_type': "\
            "<MassType.mono: 0>, 'radical': False, 'fragment_ions': None}>"
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
