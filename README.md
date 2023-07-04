# biodatatypes
Pure-Python package for handling biological sequence datatypes as Enum objects.

## Installation

```bash
pip install biodatatypes
```

## Basic usage

```python
from biodatatypes import Nucleotide, AminoAcid, Codon
from biodatatypes import NucleotideSequence, AminoAcidSequence, CodonSequence

# Nucleotide
nucleotide_a = Nucleotide['A']
nucleotide_c = Nucleotide.C
nucleotide_g = Nucleotide.from_str('G')
nucleotide_t = Nucleotide(4)
gap = Nucleotide.from_str('-')
also_gap = Nucleotide.Gap

# AminoAcid
amino_acid_ala = AminoAcid['Ala']
amino_acid_arg = AminoAcid.Arg
amino_acid_asn = AminoAcid.from_str('N')
amino_acid_gly = AminoAcid.from_str('Gly')
amino_acid_asp = AminoAcid(4)
stop = AminoAcid['Stop']
also_stop = AminoAcid.from_str('Ter')
also_also_stop = AminoAcid.from_str('*')

# Codon
codon_gca = Codon['GCA']
codon_gcg = Codon.GCG
codon_gct = Codon.from_str('GCT')
codon_aat = Codon(4)
codon_atg = Codon.start_codon()

# NucleotideSequence
nucleotide_sequence = NucleotideSequence.from_str('ATGAAACGATAG')
gapped_nucleotide_sequence = NucleotideSequence.from_str('ATG-AACGA--AG')
masked_nucleotide_sequence = NucleotideSequence.from_str('ATG#AA##ATAG')
print(nucleotide_sequence)  # ATGAAACGATAG
print(repr(nucleotide_sequence))  # ATGAAACGATAG
print(gapped_nucleotide_sequence)  # ATG-AACGA--AG
print(masked_nucleotide_sequence)  # ATG#AA##ATAG

# AminoAcidSequence
amino_acid_sequence = AminoAcidSequence.from_str('ACDE')
gapped_amino_acid_sequence = AminoAcidSequence.from_str('A-C-E')
masked_amino_acid_sequence = AminoAcidSequence.from_str('A#DE')
print(amino_acid_sequence)  # ACDE
print(gapped_amino_acid_sequence)  # A-C-E
print(masked_amino_acid_sequence)  # A#DE

# CodonSequence
codon_sequence = CodonSequence.from_str('ATGAAACGATAG')
gapped_codon_sequence = CodonSequence.from_str('ATGAAA---CGATAG')
masked_codon_sequence = CodonSequence.from_str('ATG###CGATAG')
print(codon_sequence)  # ATGAAACGATAG
print(repr(codon_sequence))  # ATG AAA CGA TAG
print(gapped_codon_sequence)  # ATGAAA---CGATAG
print(masked_codon_sequence)  # ATG###CGATAG
```

## Making custom datatypes

While default datatypes `Nucleotide`, `AminoAcid`, and `Codon` are provided, it is possible to create custom nucleotide, amino acid, and codon datatypes by subclassing `NucleotideEnum`, `AminoAcidEnum`, and `CodonEnum` respectively.

### Custom nucleotide enum

```python
from biodatatypes.unit.base import NucleotideEnum, AminoAcidEnum, CodonEnum
from biodatatypes.unit.mixins import GapTokenMixin, MaskTokenMixin

# Create a custom nucleotide enum
class MyNucleotide(GapTokenMixin, NucleotideEnum):
    A = 1
    C = 2
    G = 3
    T = 4
    Gap = 5

my_a = MyNucleotide['A']
my_c = MyNucleotide.C
my_g = MyNucleotide.from_str('G')
my_t = MyNucleotide(4)
my_gap = MyNucleotide.from_str('-')

# Use methods inherited from NucleotideEnum
print(my_a.is_purine())  # True
print(my_c.is_purine())  # False
print(my_a.is_gap())  # False
print(my_gap.is_gap())  # True
print(my_a.is_standard())  # True
print(my_t.to_onehot())  # [0, 0, 0, 0, 1]
print(my_c.to_complement())  # G
```

`NucleotideEnum` contains methods for handling standard nucleotides.

To extend functionality to handle gaps, use `GapTokenMixin` together with the corresponding enum class (e.g. `NucleotideEnum` for nucleotides). 
`GapTokenMixin` adds `is_gap()` method to check if a token is a gap. 
It expects that the gap token enum is named `Gap` (case-sensitive).

To extend functionality to handle masks, use `MaskTokenMixin` together with the corresponding enum class. 
`MaskTokenMixin` adds `is_mask()` method to check if a token is a mask.
It expects that the mask token enum is named `Mask` (case-sensitive).

To extend functionality to handle unspecified "other" nucleotides, use `OtherTokenMixin` together with the corresponding enum class. 
`OtherTokenMixin` adds `is_other()` method to check if a token is an unspecified "other" nucleotide.
It expects that the "other" token enum is named `Other` (case-sensitive).

To extend functionality to handle gaps, masks, and unspecified "other" nucleotides, use `SpecialTokenMixin` together with the corresponding enum class. 
`SpecialTokenMixin` adds `is_special()` method to check if a token is a gap, mask, or other.
The same enum name requirements apply as above.

### Custom amino acid enum

Similarly, custom amino acid enums can be created by subclassing `AminoAcidEnum`. 
Amino acid enums are expected to use the the IUPAC three-letter amino acid codes as enum names.

If the termination signal token is included, it is expected to be named `Stop` (case-sensitive).
The termination token has a one-letter code of `*` and a three-letter code of `Ter`.

```python
from biodatatypes.unit.base import AminoAcidEnum
from biodatatypes.unit.mixins import GapTokenMixin, MaskTokenMixin


# Create a custom amino acid enum
# Add mixins to extend functionality when non-standard amino acid tokens (gap, mask) are expected
class MyAminoAcid(MaskTokenMixin, GapTokenMixin, AminoAcidEnum):
    Ala = 1
    Arg = 2
    Asn = 3
    Asp = 4
    Gap = 5
    Mask = 6
    Stop = 7

my_ala = MyAminoAcid['Ala']
my_arg = MyAminoAcid.Arg
my_asn = MyAminoAcid.from_str('N')
my_gap = MyAminoAcid.from_str('-')
my_mask = MyAminoAcid.Mask
my_stop = MyAminoAcid['Stop']

print(my_asn.is_polar())  # True
print(my_ala.is_polar())  # False
print(my_gap.is_gap())  # True
print(my_ala.is_gap())  # False
print(my_mask.is_mask())  # True
print(my_asn.has_amide())  # True
print(my_asn.to_one_letter())  # N
```

The same mixins used for `NucleotideEnum` can be used for `AminoAcidEnum` to extend functionality to handle gaps, masks, and unspecified "other" amino acids.
The same enum name requirements apply as previously mentioned for making custom nucleotide enums.


### Custom codon enum

When creating a custom codon enum, aside from specifying the enumerations using the IUPAC three-letter codon codes, it is also necessary to specify associated nucleotide and amino acid enums by overriding the `nucleotide_class` and `aminoacid_class` getter methods.

The `nucleotide_class` property is used when calling `from_nucleotides` and `to_nucleotides()` methods to convert the codon to and from a triplet nucleotide sequence.
The `aminoacid_class` property is used when calling the `translate()` method to translate the codon to an amino acid based on the standard genetic code.
It is not necessary to create custom nucleotide and amino acid enums for this purpose as the included `Nucleotide` and `AminoAcid` can be used, but it is necessary to specify the corresponding enum classes.

```python
from biodatatypes import Nucleotide, AminoAcid
from biodatatypes.unit.base import CodonEnum
from biodatatypes.unit.mixins import GapTokenMixin, MaskTokenMixin


# Create a custom codon enum
# Add mixins to extend functionality when non-standard codon tokens (gap, mask) are expected
class MyCodon(MaskTokenMixin, GapTokenMixin, CodonEnum):
    GCA = 1
    GCG = 2
    GCT = 3
    TAG = 4
    Gap = 5
    Mask = 6
    ATG = 7

    @property
    def nucleotide_class(self):
        return Nucleotide

    @property
    def aminoacid_class(self):
        return MyAminoAcid

my_gca = MyCodon['GCA']
my_gcg = MyCodon.GCG
my_gct = MyCodon.from_str('GCT')
my_gap = MyCodon.from_str('---')
my_mask = MyCodon.Mask
my_atg = MyCodon.start_codon()

print(my_gca.is_fourfold_degenerate())  # True, GCA, GCC, GCG, GCT encode Ala
print(my_atg.is_start_codon())  # True
print(my_atg.is_stop_codon())  # False
print(my_gap.is_gap())  # True
print(my_mask.is_mask())  # True
print(my_gca.translate())  # A

```

The same mixins can be used for `CodonEnum` to extend functionality to handle gaps, masks, and unspecified "other" codons. The same enum name requirements apply as previously mentioned for making custom nucleotide and amino acid enums.

## License
MIT License

## Author
Kent Kawashima
