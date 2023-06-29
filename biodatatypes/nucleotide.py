from typing import List
from collections.abc import Sequence
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

    def __str__(self) -> str:
        return self.to_one_letter()

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
        """Convert a string to a Nucleotide.
        
        Parameters
        ----------
        s : str
            A string representing a nucleotide.
            
        Returns
        -------
        Nucleotide
            A Nucleotide object.
            
        Raises
        ------
        ValueError
            If the string does not represent a nucleotide.
            
        Examples
        --------
        >>> Nucleotide.from_str('A')
        A
        >>> Nucleotide.from_str('a')
        A
        >>> Nucleotide.from_str('Gap')
        -
        >>> Nucleotide.from_str('N')
        N
        >>> Nucleotide.from_str('X')
        Traceback (most recent call last):
        ValueError: Invalid nucleotide: X
        """
        try:
            return cls[s.capitalize()]
        except KeyError:
            raise ValueError(f'Invalid nucleotide: {s}')

    @classmethod
    def from_onehot(cls, onehot:Sequence[int]) -> 'Nucleotide':
        """Convert a onehot encoding to a Nucleotide.
        
        Parameters
        ----------
        onehot : Sequence[int]
            A onehot encoding of a nucleotide.
            
        Returns
        -------
        Nucleotide
            A Nucleotide object.
            
        Examples
        --------
        >>> Nucleotide.from_onehot([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        A
        >>> Nucleotide.from_onehot([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        -
        >>> Nucleotide.from_onehot([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        N
        >>> Nucleotide.from_onehot([0, 0, 0, 0, 0, 0, 0, 1])
        Traceback (most recent call last):
        ValueError: Invalid onehot value: [0, 0, 0, 0, 0, 0, 0, 1]
        >>> Nucleotide.from_onehot([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        Traceback (most recent call last):
        ValueError: Invalid onehot value: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        >>> Nucleotide.from_onehot([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        Traceback (most recent call last):
        ValueError: Invalid onehot value: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        """
        if (len(onehot) != len(cls)) or (sum(onehot) != 1):
            raise ValueError(f'Invalid onehot value: {onehot}')
        try:
            return cls(onehot.index(1))
        except ValueError:
            raise ValueError(f'Invalid onehot value: {onehot}')

    # Formatting methods
    
    def to_onehot(self) -> Sequence[int]:
        """Convert a Nucleotide to a onehot encoding.
        
        Returns
        -------
        Sequence[int]
            A onehot encoding of the nucleotide.
            
        Examples
        --------
        >>> Nucleotide['A'].to_onehot()
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        >>> Nucleotide['Gap'].to_onehot()
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        >>> Nucleotide['N'].to_onehot()
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        """
        onehot = [0] * len(self.__class__)
        onehot[self.value] = 1
        return onehot
    
    def to_one_letter(self) -> str:
        """Convert a Nucleotide to a one letter code.
        
        Returns
        -------
        str
            A one letter code for the nucleotide.
            
        Examples
        --------
        >>> a = Nucleotide['A']
        >>> a.to_one_letter()
        'A'
        >>> gap = Nucleotide['Gap']
        >>> gap.to_one_letter()
        '-'
        >>> n = Nucleotide['N']
        >>> n.to_one_letter()
        'N'
        """
        return TO_ONE_LETTER_TOKEN[self.name]

    # Symbol type

    def is_masked(self):
        """Return True if the nucleotide is masked.
        
        Returns
        -------
        bool
            True if the nucleotide is masked, False otherwise.
        
        See Also
        --------
        is_special : Checks if it is a special token
        
        Examples
        --------
        >>> Nucleotide['Mask'].is_masked()
        True
        >>> Nucleotide['A'].is_masked()
        False
        >>> Nucleotide['Gap'].is_masked()
        False
        """
        return self == Nucleotide.Mask

    def is_gap(self):
        """Return True if the nucleotide is a gap.
        
        Returns
        -------
        bool
            True if the nucleotide is a gap, False otherwise.
        
        See Also
        --------
        is_special : Checks if it is a special token
        
        Examples
        --------
        >>> Nucleotide['Gap'].is_gap()
        True
        >>> Nucleotide['A'].is_gap()
        False
        >>> Nucleotide['Mask'].is_gap()
        False
        """
        return self == Nucleotide.Gap
    
    def is_special(self):
        """Return True if it represents a special sequence.
        
        Returns
        -------
        bool
            True if the nucleotide token represents a gap, mask or special character.
        
        See Also
        --------
        is_masked : Checks if it is a masked token
        is_gap : Checks if it is a gap token
        
        Examples
        --------
        >>> gap = Nucleotide['Gap']
        >>> gap.is_special()
        True
        >>> mask = Nucleotide['Mask']
        >>> mask.is_special()
        True
        >>> special = Nucleotide['Special']
        >>> special.is_special()
        True
        >>> a = Nucleotide['A']
        >>> a.is_special()
        False
        >>> n = Nucleotide['N']
        >>> n.is_special()
        False
        """
        return (
            self.name == 'Gap' or
            self.name == 'Mask' or
            self.name == 'Special')
    
    def is_any(self):
        """Return True if it represents any nucleotide.
        
        Returns
        -------
        bool
            True if the nucleotide token represents any nucleotide.
        
        See Also
        --------
        is_ambiguous : Checks if it is an ambiguous token
        is_degenerate : Checks if it is a degenerate token
        
        Examples
        --------
        >>> n = Nucleotide['N']
        >>> n.is_any()
        True
        >>> gap = Nucleotide['Gap']
        >>> gap.is_any()
        False
        >>> mask = Nucleotide['Mask']
        >>> mask.is_any()
        False
        >>> a = Nucleotide['A']
        >>> a.is_any()
        False
        """
        return self.name == 'N'
   
    def is_ambiguous(self):
        """Return True if it represents an ambiguous nucleotide.
        
        Returns
        -------
        bool
            True if the nucleotide token represents an ambiguous nucleotide.
        
        See Also
        --------
        is_degenerate : Checks if it is a degenerate token
        
        Examples
        --------
        >>> Nucleotide['N'].is_ambiguous()
        True
        >>> Nucleotide['R'].is_ambiguous()
        True
        >>> Nucleotide['Gap'].is_ambiguous()
        False
        >>> Nucleotide['Mask'].is_ambiguous()
        False
        >>> Nucleotide['A'].is_ambiguous()
        False
        """
        return self.name in DEGENERATE_NUCLEOTIDES

    def is_degenerate(self):
        """Return True if it represents a degenerate nucleotide.
        
        Returns
        -------
        bool
            True if the nucleotide token represents a degenerate nucleotide.
        
        See Also
        --------
        is_ambiguous : Checks if it is an ambiguous token
        
        Examples
        --------
        >>> Nucleotide['N'].is_degenerate()
        True
        >>> Nucleotide['R'].is_degenerate()
        True
        >>> Nucleotide['Gap'].is_degenerate()
        False
        """
        return self.is_ambiguous()
 
    # Properties
    
    def is_purine(self):
        """Return True if the nucleotide is a purine.
        
        Returns
        -------
        bool
            True if the nucleotide is a purine, False otherwise.
            
        Examples
        --------
        >>> Nucleotide['A'].is_purine()
        True
        >>> Nucleotide['G'].is_purine()
        True
        >>> Nucleotide['R'].is_purine()  # R is A or G
        True
        >>> Nucleotide['C'].is_purine()
        False
        >>> Nucleotide['N'].is_purine()
        False
        >>> Nucleotide['Gap'].is_purine()
        False
        """
        return self == Nucleotide.A or self == Nucleotide.G or self == Nucleotide.R

    def is_pyrimidine(self):
        """Return True if the nucleotide is a pyrimidine.
        
        Returns
        -------
        bool
            True if the nucleotide is a pyrimidine, False otherwise.
            
        Examples
        --------
        >>> Nucleotide['C'].is_pyrimidine()
        True
        >>> Nucleotide['T'].is_pyrimidine()
        True
        >>> Nucleotide['Y'].is_pyrimidine()  # Y is C or T
        True
        >>> Nucleotide['A'].is_pyrimidine()
        False
        >>> Nucleotide['N'].is_pyrimidine()
        False
        >>> Nucleotide['Gap'].is_pyrimidine()
        False
        """
        return self == Nucleotide.C or self == Nucleotide.T or self == Nucleotide.Y

    def is_strong(self):
        """Return True if the nucleotide is a nucleotide involved in strong triple H-bond interaction.
        
        Returns
        -------
        bool
            True if the nucleotide is C or G, False otherwise.
            
        Examples
        --------
        >>> Nucleotide['C'].is_strong()
        True
        >>> Nucleotide['G'].is_strong()
        True
        >>> Nucleotide['S'].is_strong()  # S is C or G
        True
        >>> Nucleotide['A'].is_strong()
        False
        >>> Nucleotide['N'].is_strong()
        False
        >>> Nucleotide['Gap'].is_strong()
        False
        """
        return self == Nucleotide.G or self == Nucleotide.C or self == Nucleotide.S

    def is_weak(self):
        """Return True if the nucleotide is a nucleotide involved in weak double H-bond interaction.
        
        Returns
        -------
        bool
            True if the nucleotide is A or T, False otherwise.
        
        Examples
        --------
        >>> Nucleotide['A'].is_weak()
        True
        >>> Nucleotide['T'].is_weak()
        True
        >>> Nucleotide['W'].is_weak()  # W is A or T
        True
        >>> Nucleotide['C'].is_weak()
        False
        >>> Nucleotide['N'].is_weak()
        False
        >>> Nucleotide['Gap'].is_weak()
        False
        """
        return self == Nucleotide.A or self == Nucleotide.T or self == Nucleotide.W

    def is_amino(self):
        return self == Nucleotide.A or self == Nucleotide.C or self == Nucleotide.M

    def is_keto(self):
        return self == Nucleotide.G or self == Nucleotide.T or self == Nucleotide.K


# Constants

TO_ONE_LETTER_TOKEN = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'U',
    'Gap': '-',
    'Mask': '#',
    'Special': '@',
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
DEGENERATE_NUCLEOTIDES = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
