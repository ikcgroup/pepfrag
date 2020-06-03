Usage
=====

Peptide Construction
--------------------

:mod:`~pepfrag` provides one key public class: :class:`~pepfrag.Peptide`.
This class includes public methods for calculating the mass of the peptide,
including any configured modifications (:class:`~pepfrag.ModSite` s), and the
peptide fragment ions, with configurable neutral losses.

A :class:`~pepfrag.Peptide` can be constructed from its amino acid sequence,
charge state and modifications, for example:

.. code-block:: python

    from pepfrag import ModSite, Peptide

    peptide = Peptide(
        "ABCMPK",
        2,
        (ModSite(15.994915, 4, "Oxidation"), ModSite(304.20536, "nterm", "iTRAQ8plex"))
    )


:class:`~pepfrag.Peptide` modifications are defined using a sequence of
:class:`~pepfrag.ModSite` instances.

Additional keyword arguments are available, allowing the use of average masses instead
of monoisotopic masses and introducing radical peptide fragment generation.

Fragment Generation
-------------------

Fragment ions can be generated using the :func:`~pepfrag.Peptide.fragment` method;
for efficiency when the same :class:`~pepfrag.Peptide` instance is used repeatedly,
the resulting fragments are cached in the :attr:`~pepfrag.Peptide.fragment_ions` attribute.
This cache is invalidated if the instance attributes are changed.

The generated fragment ions can be customized using the ``ion_types`` argument to
:func:`~pepfrag.Peptide.fragment`, which takes a dictionary mapping the desired
:class:`~pepfrag.IonType` s to their planned neutral losses. The default is:

.. code-block:: python

    from pepfrag import IonType

    DEFAULT_IONS = {
        IonType.precursor: ['H2O', 'NH3', 'CO2'],
        IonType.imm: [],
        IonType.b: ['H2O', 'NH3', 'CO'],
        IonType.y: ['NH3', 'H2O'],
        IonType.a: [],
        IonType.c: [],
        IonType.z: []
    }

The generated ions can be changed by providing a custom ``ion_types`` dictionary
when calling :func:`~pepfrag.Peptide.fragment`, for example:

.. code-block:: python

    from pepfrag import IonType, Peptide

    peptide = Peptide('AMYK', 2, [])
    peptide.fragment(ion_types={
        IonType.precursor: [],
        IonType.b: ['NH3'],
        IonType.y: ['H2O']
    })

outputs the following fragment ions, including precursor ions, `b` ions with `NH3`
losses and `y` ions with `H2O` losses:

.. code-block:: python

    [
        (72.044390252029, 'b1[+]', 1),
        (55.01784115090901, '[b1-NH3][+]', 1),
        (147.11280416609898, 'y1[+]', 1),
        (129.10223948206897, '[y1-H2O][+]', 1),
        (203.084875340499, 'b2[+]', 2),
        (186.058326239379, '[b2-NH3][+]', 2),
        (310.17613269973896, 'y2[+]', 2),
        (292.16556801570897, '[y2-H2O][+]', 2),
        (366.14820387413897, 'b3[+]', 3),
        (349.121654773019, '[b3-NH3][+]', 3),
        (183.57774017050897, 'b3[2+]', 3),
        (175.06446561994898, '[b3-NH3][2+]', 3),
        (441.21661778820896, 'y3[+]', 3),
        (423.206053104179, '[y3-H2O][+]', 3),
        (221.11194712754397, 'y3[2+]', 3),
        (212.10666478552898, '[y3-H2O][2+]', 3),
        (512.253731573359, '[M+H][+]', 4),
        (256.63050402011896, '[M+H][2+]', 4)
    ]

Customizing Neutral Losses
^^^^^^^^^^^^^^^^^^^^^^^^^^

:mod:`~pepfrag` includes a number of common neutral losses available using only their
string names. These are: `NH3`, `H2O`, `CO2` and `CO`.

Additional neutral losses can be specified using a tuple of `(label, mass)`.
For example:

.. code-block:: python

    from pepfrag import IonType

    ion_types = {
        IonType.b: [('testLoss1', 17.04), 'NH3]
    }

This would generate `b` ions, along with `b-testLoss1` and `b-NH3` fragment ions.
