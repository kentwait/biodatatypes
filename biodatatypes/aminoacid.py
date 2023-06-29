from typing import List
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
    # Special tokens
    Ter = 21 # Terminator, *
    Mask = 22  # Masked, ?
    Special = 23  # Special, @
    Gap = 0  # Single gap, -
    # Non-standard amino acids
    Sec = 24  # Selenocysteine, U
    Pyl = 25  # Pyrrolysine, O
    # Ambiguous amino acids
    Asx = 26  # Asn or Asp, B
    Glx = 27  # Gln or Glu, Z
    Xle = 28  # Leu or Ile, J
    Xaa = 29  # Unknown or Ambiguous, X
    
    def __str__(self) -> str:
        return self.to_one_letter_code()
    
    def __repr__(self) -> str:
        return self.__str__()
    
    def __eq__(self, other):
        if isinstance(other, AminoAcid):
            return self.value == other.value
        elif isinstance(other, str):
            return str(self) == other
        else:
            return False
        
    def __hash__(self):
        return hash(str(self))
    
    @classmethod
    def from_str(cls, s: str) -> 'AminoAcid':
        if len(s) == 1:
            return cls.from_one_letter(s)
        elif len(s) == 3:
            return cls.from_three_letter(s)
        else:
            raise ValueError(f'Invalid string: {s}')
    
    @classmethod
    def from_one_letter(cls, one_letter_code: str) -> 'AminoAcid':
        try:
            return cls[FROM_ONE_LETTER_TOKEN[one_letter_code.upper()]]
        except KeyError:
            raise ValueError(f'Invalid one-letter code: {one_letter_code}')
    
    @classmethod
    def from_three_letter(cls, three_letter_code: str) -> 'AminoAcid':
        try:
            return cls[FROM_THREE_LETTER_TOKEN[three_letter_code.capitalize()]]
        except KeyError:
            raise ValueError(f'Invalid three-letter code: {three_letter_code}')
    
    @classmethod
    def from_onehot(cls, onehot: List[int]) -> 'AminoAcid':
        try:
            return cls(onehot.index(1))
        except ValueError:
            raise ValueError(f'Invalid onehot value: {onehot}')
    
    # Formatting methods
    
    def to_one_letter_code(self) -> str:
        return TO_ONE_LETTER_TOKEN[self.name]
    
    def to_three_letter_code(self) -> str:
        return self.name
    
    def to_onehot(self) -> List[int]:
        onehot = [0] * len(self.__class__)
        onehot[self.value] = 1
        return onehot
    
    # Symbol type
    
    def is_standard(self) -> bool:
        return self.name in STANDARD_AA
    
    def is_nonstandard(self) -> bool:
        return (
            self.name == 'Sec' or 
            self.name == 'Pyl')
        
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
            self.name == 'Gap' or
            self.name == 'Special'
        )
    
    def is_any(self):
        return self.name == 'Xaa'

    def is_unknown(self):
        return self.is_any()
    
    def is_other(self):
        return self.is_any()
   
    def is_ambiguous(self) -> bool:
        return (
            self.name == 'Asx' or 
            self.name == 'Glx' or 
            self.name == 'Xle' or 
            self.name == 'Xaa')

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
# One-letter lookup tables
FROM_ONE_LETTER_TOKEN = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
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
    'U': 'Sec',
    'O': 'Pyl',
    '#' : 'Mask',
    '*': 'Ter',
    '@': 'Special',
    '-' : 'Gap',
    'B': 'Asx',
    'Z': 'Glx',
    'J': 'Xle',
    'X': 'Xaa',
}
TO_ONE_LETTER_TOKEN = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
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
    'Sec': 'U',
    'Pyl': 'O',
    'Mask': '#',
    'Special': '@',
    'Gap': '-',
    'Ter': '*',
    'Asx': 'B',
    'Glx': 'Z',
    'Xle': 'J',
    'Xaa': 'X',
}
# Three-letter lookup tables
FROM_THREE_LETTER_TOKEN = {
    'Ala': 'Ala',
    'Arg': 'Arg',
    'Asn': 'Asn',
    'Asp': 'Asp',
    'Cys': 'Cys',
    'Glu': 'Glu',
    'Gln': 'Gln',
    'Gly': 'Gly',
    'His': 'His',
    'Ile': 'Ile',
    'Leu': 'Leu',
    'Lys': 'Lys',
    'Met': 'Met',
    'Phe': 'Phe',
    'Pro': 'Pro',
    'Ser': 'Ser',
    'Thr': 'Thr',
    'Trp': 'Trp',
    'Tyr': 'Tyr',
    'Val': 'Val',
    'Sec': 'Sec',
    'Pyl': 'Pyl',
    'Msk': 'Mask',
    'Spl': 'Special',
    'Gap': 'Gap',
    'Ter': 'Ter',
    'Asx': 'Asx',
    'Glx': 'Glx',
    'Xle': 'Xle',
    'Xaa': 'Xaa',
}
TO_THREE_LETTER_TOKEN = {
    'Ala': 'Ala',
    'Arg': 'Arg',
    'Asn': 'Asn',
    'Asp': 'Asp',
    'Cys': 'Cys',
    'Glu': 'Glu',
    'Gln': 'Gln',
    'Gly': 'Gly',
    'His': 'His',
    'Ile': 'Ile',
    'Leu': 'Leu',
    'Lys': 'Lys',
    'Met': 'Met',
    'Phe': 'Phe',
    'Pro': 'Pro',
    'Ser': 'Ser',
    'Thr': 'Thr',
    'Trp': 'Trp',
    'Tyr': 'Tyr',
    'Val': 'Val',
    'Sec': 'Sec',
    'Pyl': 'Pyl',
    'Mask': 'Msk',
    'Special': 'Spl',
    'Gap': 'Gap',
    'Ter': 'Ter',
    'Asx': 'Asx',
    'Glx': 'Glx',
    'Xle': 'Xle',
    'Xaa': 'Xaa',
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
