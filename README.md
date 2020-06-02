# pepfrag

![Python package](https://github.com/ikcgroup/pepfrag/workflows/Python%20package/badge.svg)
[![codecov](https://codecov.io/gh/ikcgroup/pepfrag/branch/master/graph/badge.svg)](https://codecov.io/gh/ikcgroup/pepfrag)


pepfrag is a library for generating possible dissociation fragment ions of peptides
in tandem mass spectrometry experiments.

# Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [License](#license)

# Installation

`pepfrag` can be installed from [PyPI](https://pypi.org/project/pepfrag/).

### Compatibility

`pepfrag` is written using Python 3 and should be compatible with most 
operating systems. The package has been tested on
- Windows 10
- MacOS 10.15
- Ubuntu 18.04 LTS

Because `pepfrag` includes C/C++ extensions, installation requires the 
presence of a C++ 11 compatible compiler on your machine.

### Instructions

1. Install Python 3 (>= version 3.6).
2. Install the latest version of the `pepfrag` library using `pip`:
```shell script
pip install pepfrag
```

# Usage

### `Peptide` Construction

`pepfrag` provides one key public class: `Peptide`. This class includes public methods
for computing the mono/average mass of the peptide, including any configured modifications
(`ModSite`s), and the peptide fragment ions, with configurable neutral losses.

A `Peptide` can be constructed from its amino acid sequence, charge state and modifications,
for example:
```python
from pepfrag import MassType, ModSite, Peptide

peptide = Peptide(
    "ABCMPK", 
    2, 
    (ModSite(15.994915, 4, "Oxidation"), ModSite(304.20536, "nterm", "iTRAQ8plex")),
    mass_type=MassType.mono
)
```

`Peptide` modifications are defined using a sequence of `ModSite` instances, which
are `namedtuples` defined by the mass of the modification (float), the site of the
modification (string for terminal modifications or 1-indexed integer otherwise) and
the name of the modification.

The `Peptide` constructor has two keyword parameters:
- `mass_type`:
    - Description: The type of mass to calculate.
    - Type: `MassType` enumeration.
    - Default: `MassType.mono`.
- `radical`:
    - Description: Flag indicating whether radical cation fragments should be 
    generated.
    - Type: bool.
    - Default: `False`.
    
### Fragment Generation
    
Fragment ions can be generated using the `fragment` method; for efficiency when the
same `Peptide` instance is used repeatedly, the resulting fragments are cached on the
instance. This cache is invalidated if the instance attributes are changed.

The `fragment` method has two keyword parameters:
- `ion_types`:
    - Description: The types of fragment ion species to generate. The default
    setup would generate precursor, immonium, b, y, a, c and z ions.
    - Type: dictionary mapping `IonType` enumeration instances to possible 
    neutral loss species, represented as strings for the most common neutral losses
    (i.e. those configured in `pepfrag`, namely `'H2O'`, `'NH3'`, `'CO2'` and `'CO'`).
    For additional losses, these can be input as tuples of 
    `(name (str), mass (float))`.
    - Default: 
    ```python
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
    ```
  - Example:
  ```python
  from pepfrag import IonType, Peptide
  
  peptide = Peptide('APYSMLK', 2, [])
  
  peptide.fragment(ion_types={
      IonType.b: ['NH3', ('customLoss', 17.04)]
  })
  ```
- `force`:
    - Description: Flag indicating whether fragments should be forcibly regenerated,
    *i.e.* bypassing the cached ions.
    - Type: bool.
    - Default: `False`.

# License

`pepfrag` is released under the terms of the [MIT license](https://github.com/ikcgroup/pepfrag/blob/master/LICENSE).