from enum import Enum


class Nucleotide(Enum):
    A = 1
    C = 2
    G = 3
    T = 4
    U = 5
    R = 6  # A or G, purine
    Y = 7  # C or T, pyrimidine
    S = 8  # G or C, strong
    W = 9  # A or T, weak
    K = 10  # G or T, keto
    M = 11  # A or C, amino
    B = 12  # C or G or T, not A
    D = 13  # A or G or T, not C
    H = 14  # A or C or T, not G
    V = 15  # A or C or G, not T
    N = 16  # A or C or G or T, any
    X = 17  # A or C or G or T, any
    Gap = 0  # Gap

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
            return cls[s.upper()]
        except KeyError:
            raise ValueError(f'Invalid nucleotide: {s}')

    def is_gap(self):
        return self == Nucleotide.Gap
    
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
    
    def is_any(self):
        return self == Nucleotide.N or self == Nucleotide.X
    
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
            self == Nucleotide.N or
            self == Nucleotide.X)
