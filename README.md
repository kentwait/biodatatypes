# biodatatypes
Biological sequence datatypes

## Installation

```bash
pip install biodatatypes
```

## Usage

```python
from biodatatypes import Nucleotide, AminoAcid, Codon
from biodatatypes import NucleotideSequence, AminoAcidSequence, CodonSequence

# Nucleotide
nucleotide_a = Nucleotide['A']
nucleotide_c = Nucleotide.c
nucleotide_g = Nucleotide.from_str('G')
nucleotide_t = Nucleotide(4)
gap = Nucleotide['-']
also_gap = Nucleotide.Gap

# AminoAcid
amino_acid_ala = AminoAcid['Ala']
amino_acid_arg = AminoAcid.Arg
amino_acid_asn = AminoAcid.from_str('Asn')
amino_acid_asp = AminoAcid(4)
stop = AminoAcid['Stop']
also_stop = AminoAcid.Stop

# Codon
codon_gca = Codon['GCA']
codon_gcg = Codon.GCG
codon_gct = Codon.from_str('GCT')
codon_aat = Codon(4)
stop = Codon['Stop']
also_stop = Codon.Stop

# NucleotideSequence
nucleotide_sequence = NucleotideSequence.from_str('ACGT')
gapped_nucleotide_sequence = NucleotideSequence.from_str('A-CG-T')
masked_nucleotide_sequence = NucleotideSequence.from_str('A#GT')

# AminoAcidSequence
amino_acid_sequence = AminoAcidSequence.from_str('ACDE')
gapped_amino_acid_sequence = AminoAcidSequence.from_str('A-C-E')
masked_amino_acid_sequence = AminoAcidSequence.from_str('A#DE')

# CodonSequence
codon_sequence = CodonSequence.from_str('ATGAAACGATAG')
gapped_codon_sequence = CodonSequence.from_str('ATGAAA---CGATAG')
masked_codon_sequence = CodonSequence.from_str('ATG###CGATAG')
```

## License
MIT License

## Author
Kent Kawashima
