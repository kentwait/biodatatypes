from enum import Enum


class AminoAcid(Enum):
    Ala = 1
    Arg = 2
    Asn = 3
    Asp = 4
    Cys = 5
    Glu = 6
    Gln = 7
    Gly = 8
    His = 9
    Ile = 10
    Leu = 11
    Lys = 12
    Met = 13
    Phe = 14
    Pro = 15
    Ser = 16
    Thr = 17
    Trp = 18
    Tyr = 19
    Val = 20
    # Non-standard amino acids
    Sec = 21  # Selenocysteine, U
    Pyl = 22  # Pyrrolysine, O
    # Ambiguous amino acids
    Asx = 23  # Asn or Asp, B
    Glx = 24  # Gln or Glu, Z
    Xle = 25  # Leu or Ile, J
    Xaa = 26  # Unknown or Ambiguous, X
    # Special tokens
    Mask = 27  # Masked, ?
    Ter = 28  # Terminator, *
    Gap = 0  # Single gap, -
    
    def __str__(self) -> str:
        return self.to_one_letter_code()
    
    @classmethod
    def from_one_letter_code(cls, one_letter_code: str) -> 'AminoAcid':
        return cls[ONE_TO_THREE_LETTER_AA[one_letter_code]]
    
    @classmethod
    def from_three_letter_code(cls, three_letter_code: str) -> 'AminoAcid':
        return cls[three_letter_code]
    
    def to_one_letter_code(self) -> str:
        return THREE_TO_ONE_LETTER_AA[self.name]
    
    def to_three_letter_code(self) -> str:
        return self.name
    
    def is_standard(self) -> bool:
        return self.name in STANDARD_AA
    
    def is_nonstandard(self) -> bool:
        return (
            self.name == 'Sec' or 
            self.name == 'Pyl')
    
    def is_ambiguous(self) -> bool:
        return (
            self.name == 'Asx' or 
            self.name == 'Glx' or 
            self.name == 'Xle' or
            self.name == 'Mask')
    
    def is_unknown(self) -> bool:
        return self.name == 'Xaa'
    
    def is_other(self) -> bool:
        return self.name == 'Xaa'
    
    def is_terminator(self) -> bool:
        return self.name == 'Ter'
    
    def is_masked(self) -> bool:
        return self.name == 'Mask'
    
    def is_gap(self) -> bool:
        return self.name == 'Gap'
    
    def is_special(self) -> bool:
        return (
            self.name == 'Mask' or
            self.name == 'Ter' or
            self.name == 'Gap'
        )
    
    # Properties
    
    def is_polar(self) -> bool:
        return self.name in POLAR_AA

    def is_nonpolar(self) -> bool:
        return not self.is_polar()
    
    def is_acidic(self) -> bool:
        return self.name in ACIDIC_AA
    
    def is_basic(self) -> bool:
        return self.name in BASIC_AA
    
    def is_aromatic(self) -> bool:
        return self.name in AROMATIC_AA
    
    def is_hydrophobic(self) -> bool:
        return self.name in HYDROPHOBIC_AA

    # Contains particular elements/functional groups
    
    def has_sulfur(self) -> bool:
        return self.name in SULFUR_AA
    
    def has_amide(self) -> bool:
        return self.name in AMIDE_AA
    
    def has_hydroxyl(self) -> bool:
        return self.name in HYDROXYL_AA

# Constants

# Naming
STANDARD_AA = [
    'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile',
    'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
]
ONE_TO_THREE_LETTER_AA = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'B': 'Asx',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'Z': 'Glx',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'W': 'Trp',
    'Y': 'Tyr',
    'V': 'Val',
    '*': 'Ter'
}
THREE_TO_ONE_LETTER_AA = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'Ter': '*',
}

# Properties
# Polar and nonpolar are disjoint
POLAR_AA = ['Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'His', 'Lys', 'Ser', 'Thr', 'Tyr']
NONPOLAR_AA = ['Ala', 'Gly', 'Ile', 'Leu', 'Met', 'Phe', 'Pro', 'Trp', 'Val']
# Acidic and basic are subsets of polar
ACIDIC_AA = ['Asp', 'Glu']
NEGATIVE_AA = ['Asp', 'Glu']
BASIC_AA = ['Arg', 'His', 'Lys']
POSITIVE_AA = ['Arg', 'His', 'Lys']
# Aromatic is a subset of nonpolar
AROMATIC_AA = ['His', 'Phe', 'Trp', 'Tyr']
# Aliphatic is a subset of nonpolar, but not aromatic
ALIPHATIC_AA = ['Ala', 'Gly', 'Ile', 'Leu', 'Pro', 'Val']
# Hydrophobic and hydrophilic are disjoint
HYDROPHOBIC_AA = ['Ala', 'Ile', 'Leu', 'Met', 'Phe', 'Pro', 'Trp', 'Val']
HYDROPHILIC_AA = ['Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'His', 'Lys', 'Ser', 'Thr', 'Tyr']
# Contains sulfur
SULFUR_AA = ['Cys', 'Met']
# Contains amide
AMIDE_AA = ['Asn', 'Gln']
# Contains hydroxyl
HYDROXYL_AA = ['Ser', 'Thr', 'Tyr']
