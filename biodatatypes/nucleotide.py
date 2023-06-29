from typing import List
from enum import Enum


class Nucleotide(Enum):
    A = 1
    C = 2
    G = 3
    T = 4
    U = 5
    # Special tokens
    Mask = 6  # Masked, #
    Special = 7  # Special, @
    Gap = 0  # Gap, -

    def __str__(self) -> str:
        if self.name == 'Gap':
            return '-'
        return self.name

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other):
        if isinstance(other, Nucleotide):
            return self.value == other.value
        elif isinstance(other, str):
            return str(self) == other
        else:
            return False

    def __hash__(self):
        return hash(str(self))

    @classmethod
    def from_str(cls, s: str) -> 'Nucleotide':
        try:
            return cls[TO_ONE_LETTER_SYMBOL[s.upper()]]
        except KeyError:
            raise ValueError(f'Invalid nucleotide: {s}')

    @classmethod
    def from_onehot(cls, onehot:List[int]) -> 'Nucleotide':
        try:
            return cls(onehot.index(1))
        except ValueError:
            raise ValueError(f'Invalid onehot value: {onehot}')

    # Formatting methods
    
    def to_onehot(self) -> List[int]:
        onehot = [0] * len(self.__class__)
        onehot[self.value] = 1
        return onehot

    # Symbol type

    def is_any(self):
        return self == Nucleotide.N

    def is_unknown(self):
        return self.is_any()
    
    def is_other(self):
        return self.is_any()
   
    def is_masked(self):
        return self == Nucleotide.Mask

    def is_gap(self):
        return self == Nucleotide.Gap
    
    def is_special(self):
        return self.is_masked() or self.is_gap()
 
    # Properties
    
    def is_purine(self):
        return self == Nucleotide.A or self == Nucleotide.G

    def is_pyrimidine(self):
        return self == Nucleotide.C or self == Nucleotide.T

    def is_strong(self):
        return self == Nucleotide.G or self == Nucleotide.C

    def is_weak(self):
        return self == Nucleotide.A or self == Nucleotide.T

    def is_amino(self):
        return self == Nucleotide.A or self == Nucleotide.C

    def is_keto(self):
        return self == Nucleotide.G or self == Nucleotide.T

  
class DegenNucleotide(Nucleotide):
    # I = 6  # Inosine
    # X = 7  # Xanthosine
    # P = 8  # Pseudouridine
    # Degenerate symbols
    R = 8  # A or G, purine
    Y = 9 # C or T, pyrimidine
    S = 10  # G or C, strong
    W = 11  # A or T, weak
    K = 12  # G or T, keto
    M = 13  # A or C, amino
    B = 14  # C or G or T, not A
    D = 15  # A or G or T, not C
    H = 16  # A or C or T, not G
    V = 17  # A or C or G, not T
    N = 18  # A or C or G or T, any

    # Symbol type
    
    # def is_standard(self):
    #     return (
    #         self == Nucleotide.A or 
    #         self == Nucleotide.C or 
    #         self == Nucleotide.G or 
    #         self == Nucleotide.T or 
    #         self == Nucleotide.U)   

    # def is_nonstandard(self):
    #     return (
    #         self == Nucleotide.I or 
    #         self == Nucleotide.X or 
    #         self == Nucleotide.P)

    def is_ambiguous(self):
        return (
            self == Nucleotide.R or 
            self == Nucleotide.Y or 
            self == Nucleotide.S or 
            self == Nucleotide.W or 
            self == Nucleotide.K or 
            self == Nucleotide.M or 
            self == Nucleotide.B or 
            self == Nucleotide.D or 
            self == Nucleotide.H or 
            self == Nucleotide.V or 
            self == Nucleotide.N)

    def is_degenerate(self):
        return self.is_ambiguous

    # Properties
    
    def is_purine(self):
        return super().is_purine() or self == Nucleotide.R

    def is_pyrimidine(self):
        return super().is_pyrimidine() or self == Nucleotide.Y

    def is_strong(self):
        return super().is_strong() or self == Nucleotide.S

    def is_weak(self):
        return super().is_weak() or self == Nucleotide.W

    def is_amino(self):
        return super().is_amino() or self == Nucleotide.M

    def is_keto(self):
        return super().is_keto() or self == Nucleotide.K

# Constants
TO_ONE_LETTER_SYMBOL = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U',
    'GAP': '-',
    'MASK': '#',
    'SPECIAL': '@',
    'R': 'R',
    'Y': 'Y',
    'S': 'S',
    'W': 'W',
    'K': 'K',
    'M': 'M',
    'B': 'B',
    'D': 'D',
    'H': 'H',
    'V': 'V',
    'N': 'N',
}
