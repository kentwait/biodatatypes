from typing import Any
from collections.abc import Sequence
from enum import Enum

from biodatatypes.constants.nucleotide import *
from biodatatypes.constants.aminoacid import *
from biodatatypes.constants.codon import *


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
        >>> Nucleotide.from_str('-')
        -
        >>> Nucleotide.from_str('N')
        N
        >>> Nucleotide.from_str('X')
        Traceback (most recent call last):
        ValueError: Invalid nucleotide: X
        """
        try:
            return cls[FROM_NUCLEOTIDE_ONE_LETTER_TOKEN[s.upper()]]
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
        return TO_NUCLEOTIDE_ONE_LETTER_TOKEN[self.name]

    # Symbol type

    def is_standard(self):
        return self.name in ['A', 'C', 'G', 'T']

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
        """Convert a string to an AminoAcid object.
        
        Parameters
        ----------
        s : str
            A string representing an amino acid. It can be a one-letter code or
            a three-letter code.

        Returns
        -------
        AminoAcid
            An AminoAcid object.
            
        Raises
        ------
        ValueError
            If the string is not a valid amino acid.
            
        Examples
        --------
        >>> AminoAcid.from_str('A')
        A
        >>> AminoAcid.from_str('Ala')
        A
        >>> AminoAcid.from_str('gap')
        -
        >>> AminoAcid.from_str('!')
        Traceback (most recent call last):
        ValueError: Invalid one-letter code: !
        """
        if len(s) == 1:
            return cls.from_one_letter(s)
        elif len(s) == 3:
            return cls.from_three_letter(s)
        else:
            raise ValueError(f'Invalid string: {s}')
    
    @classmethod
    def from_one_letter(cls, one_letter_code: str) -> 'AminoAcid':
        """Convert a one-letter code to an AminoAcid object.
        
        Parameters
        ----------
        one_letter_code : str
            A one-letter code representing an amino acid.
        
        Returns
        -------
        AminoAcid
            An AminoAcid object.
            
        Raises
        ------
        ValueError
            If the one-letter code is not a valid amino acid.
            
        Examples
        --------
        >>> AminoAcid.from_one_letter('A')
        A
        >>> AminoAcid.from_one_letter('a')
        A
        >>> AminoAcid.from_one_letter('-')
        -
        >>> AminoAcid.from_one_letter('Asp')
        Traceback (most recent call last):
        ValueError: Invalid one-letter code: Asp
        >>> AminoAcid.from_one_letter('!')
        Traceback (most recent call last):
        ValueError: Invalid one-letter code: !
        """
        try:
            return cls[FROM_AMINO_ACID_ONE_LETTER_TOKEN[one_letter_code.upper()]]
        except KeyError:
            raise ValueError(f'Invalid one-letter code: {one_letter_code}')
    
    @classmethod
    def from_three_letter(cls, three_letter_code: str) -> 'AminoAcid':
        """Convert a three-letter code to an AminoAcid object.
        
        Parameters
        ----------
        three_letter_code : str
            A three-letter code representing an amino acid.
            
        Returns
        -------
        AminoAcid
            An AminoAcid object.
            
        Raises
        ------
        ValueError
            If the three-letter code is not a valid amino acid.
            
        Examples
        --------
        >>> AminoAcid.from_three_letter('Ala')
        A
        >>> AminoAcid.from_three_letter('ala')
        A
        >>> AminoAcid.from_three_letter('ALA')
        A
        >>> AminoAcid.from_three_letter('Gap')
        -
        >>> AminoAcid.from_three_letter('XXX')
        Traceback (most recent call last):
        ValueError: Invalid three-letter code: XXX
        """
        try:
            return cls[FROM_AMINO_ACID_THREE_LETTER_TOKEN[three_letter_code.capitalize()]]
        except KeyError:
            raise ValueError(f'Invalid three-letter code: {three_letter_code}')
    
    @classmethod
    def from_onehot(cls, onehot: Sequence[int]) -> 'AminoAcid':
        """Convert a onehot vector to an AminoAcid object.
        
        Parameters
        ----------
        onehot : Sequence[int]
            A onehot vector representing an amino acid.
            
        Returns
        -------
        AminoAcid
            An AminoAcid object.
            
        Raises
        ------
        ValueError
            If the onehot vector is not a valid amino acid.
            
        Examples
        --------
        >>> AminoAcid.from_onehot([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        A
        >>> AminoAcid.from_onehot([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        -
        >>> AminoAcid.from_onehot([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        X
        >>> AminoAcid.from_onehot([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        Traceback (most recent call last):
        ValueError: Invalid onehot value: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        if (len(onehot) != len(cls)) or (sum(onehot) != 1):
            raise ValueError(f'Invalid onehot value: {onehot}')
        try:
            return cls(onehot.index(1))
        except ValueError:
            raise ValueError(f'Invalid onehot value: {onehot}')
    
    # Formatting methods
    
    def to_one_letter_code(self) -> str:
        """Convert an AminoAcid object to a one-letter code.
        
        Returns
        -------
        str
            A one-letter code representing the amino acid.
        
        Examples
        --------
        >>> aa = AminoAcid.from_str('A')
        >>> aa.to_one_letter_code()
        'A'
        >>> gap = AminoAcid.from_str('-')
        >>> gap.to_one_letter_code()
        '-'
        >>> unknown = AminoAcid.from_str('X')
        >>> unknown.to_one_letter_code()
        'X'
        """
        return TO_AMINO_ACID_ONE_LETTER_TOKEN[self.name]
    
    def to_three_letter_code(self) -> str:
        """Convert an AminoAcid object to a three-letter code.
        
        Returns
        -------
        str
            A three-letter code representing the amino acid.
            
        Examples
        --------
        >>> aa = AminoAcid.from_str('A')
        >>> aa.to_three_letter_code()
        'Ala'
        >>> gap = AminoAcid.from_str('-')
        >>> gap.to_three_letter_code()
        'Gap'
        >>> unknown = AminoAcid.from_str('X')
        >>> unknown.to_three_letter_code()
        'Xaa'
        """
        return self.name
    
    def to_onehot(self) -> Sequence[int]:
        """Convert an AminoAcid object to a onehot vector.
        
        Returns
        -------
        Sequence[int]
            A onehot vector representing the amino acid.
            
        Examples
        --------
        >>> aa = AminoAcid.from_str('A')
        >>> aa.to_onehot()
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        >>> gap = AminoAcid.from_str('-')
        >>> gap.to_onehot()
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        >>> unknown_aa = AminoAcid.from_str('Xaa')
        >>> unknown_aa.to_onehot()
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        """
        onehot = [0] * len(self.__class__)
        onehot[self.value] = 1
        return onehot
    
    # Symbol type
    
    def is_standard(self) -> bool:
        """Return True if it is a standard amino acid.
            
        Returns
        -------
        bool
            True if the amino acid token represents a standard amino acid.
            
        Examples
        --------
        >>> AminoAcid['Ala'].is_standard()
        True
        >>> AminoAcid['Sec'].is_standard()
        False
        >>> AminoAcid['Xaa'].is_standard()
        False
        """
        return self.name in STANDARD_AA
    
    def is_nonstandard(self) -> bool:
        """Return True if the amino acid is a nonstandard amino acid.
        
        Returns
        -------
        bool
            True if the amino acid is a nonstandard amino acid selenocysteine or pyrrolysine.
            
        Examples
        --------
        >>> AminoAcid['Ala'].is_nonstandard()
        False
        >>> AminoAcid['Sec'].is_nonstandard()
        True
        >>> AminoAcid['Xaa'].is_nonstandard()
        False
        """
        return (
            self.name == 'Sec' or 
            self.name == 'Pyl')
        
    def is_terminator(self) -> bool:
        """Return True if it represents a terminator sequence.
        
        Returns
        -------
        bool
            True if the amino acid token is a terminator sequence.
            
        Examples
        --------
        >>> AminoAcid['Ter'].is_terminator()
        True
        >>> AminoAcid['Ala'].is_terminator()
        False
        >>> AminoAcid['Gap'].is_terminator()
        False
        """
        return self.name == 'Ter'
    
    def is_masked(self) -> bool:
        """Return True if masked.
        
        Returns
        -------
        bool
            True if the amino acid object is masked.
        
        Examples
        --------
        >>> AminoAcid['Mask'].is_masked()
        True
        >>> AminoAcid['Ala'].is_masked()
        False
        >>> AminoAcid['Gap'].is_masked()
        False
        """
        return self.name == 'Mask'
    
    def is_gap(self) -> bool:
        """Return True if it represents a gap.
        
        Returns
        -------
        bool
            True if the amino acid token represents a gap.
        
        Examples
        --------
        >>> AminoAcid['Gap'].is_gap()
        True
        >>> AminoAcid['Ala'].is_gap()
        False
        >>> AminoAcid['Mask'].is_gap()
        False
        """
        return self.name == 'Gap'
    
    def is_special(self) -> bool:
        """Return True if it represents a special sequence.
        
        Returns
        -------
        bool
            True if the amino acid token represents a special sequence (Mask, Ter, Gap, Special).
            
        Examples
        --------
        >>> AminoAcid['Mask'].is_special()
        True
        >>> AminoAcid['Gap'].is_special()
        True
        >>> AminoAcid['Ter'].is_special()
        True
        >>> AminoAcid['Special'].is_special()
        True
        >>> AminoAcid['Ala'].is_special()
        False
        >>> AminoAcid['Xaa'].is_special()
        False
        >>> AminoAcid['Sec'].is_special()
        False
        """
        return (
            self.name == 'Mask' or
            self.name == 'Ter' or
            self.name == 'Gap' or
            self.name == 'Special'
        )
    
    def is_any(self):
        """Return True if it represents any amino acid.
        
        Returns
        -------
        bool
            True if the amino acid token represents any amino acid.
        
        Examples
        --------
        >>> AminoAcid['Xaa'].is_any()
        True
        >>> AminoAcid['Ala'].is_any()
        False
        >>> AminoAcid['Gap'].is_any()
        False
        >>> AminoAcid['Mask'].is_any()
        False
        """
        return self.name == 'Xaa'
   
    def is_ambiguous(self) -> bool:
        """Return True if it represents an ambiguous amino acid token.
        
        Returns
        -------
        bool
            True if the amino acid token represents an ambiguous amino acid token.
            
        Examples
        --------
        >>> AminoAcid['Asx'].is_ambiguous()
        True
        >>> AminoAcid['Glx'].is_ambiguous()
        True
        >>> AminoAcid['Xle'].is_ambiguous()
        True
        >>> AminoAcid['Xaa'].is_ambiguous()
        True
        >>> AminoAcid['Ala'].is_ambiguous()
        False
        >>> AminoAcid['Gap'].is_ambiguous()
        False
        """
        return (
            self.name == 'Asx' or 
            self.name == 'Glx' or 
            self.name == 'Xle' or 
            self.name == 'Xaa')

    # Properties
    
    def is_polar(self) -> bool:
        """Return True if the amino acid identity is known and is polar.
        
        Returns
        -------
        bool
            True if the amino acid token represents a polar amino acid.
            
        Examples
        --------
        >>> AminoAcid['Arg'].is_polar()
        True
        >>> AminoAcid['Ala'].is_polar()
        False
        >>> AminoAcid['Xaa'].is_polar()
        False
        """
        return self.name in POLAR_AA

    def is_nonpolar(self) -> bool:
        """Return True if the amino acid identity is known and is nonpolar.
        
        Returns
        -------
        bool
            True if the amino acid token represents a nonpolar amino acid.
            
        Examples
        --------
        >>> AminoAcid['Ala'].is_nonpolar()
        True
        >>> AminoAcid['Arg'].is_nonpolar()
        False
        >>> AminoAcid['Xaa'].is_nonpolar()
        False
        """
        return self.name in NONPOLAR_AA
    
    def is_acidic(self) -> bool:
        """Return True if the amino acid identity is known and is acidic.
        
        Returns
        -------
        bool
            True if the amino acid token represents an acidic amino acid.
            
        Examples
        --------
        >>> AminoAcid['Asp'].is_acidic()
        True
        >>> AminoAcid['Ala'].is_acidic()
        False
        >>> AminoAcid['Xaa'].is_acidic()
        False
        """
        return self.name in ACIDIC_AA
    
    def is_basic(self) -> bool:
        """Return True if the amino acid identity is known and is basic.
        
        Returns
        -------
        bool
            True if the amino acid token represents a basic amino acid.
        
        Examples
        --------
        >>> AminoAcid['Arg'].is_basic()
        True
        >>> AminoAcid['Ala'].is_basic()
        False
        >>> AminoAcid['Xaa'].is_basic()
        False
        """
        return self.name in BASIC_AA
    
    def is_aromatic(self) -> bool:
        """Return True if the amino acid identity is known and is aromatic.
        
        Returns
        -------
        bool
            True if the amino acid token represents an aromatic amino acid.
            
        Examples
        --------
        >>> AminoAcid['Phe'].is_aromatic()
        True
        >>> AminoAcid['Ala'].is_aromatic()
        False
        >>> AminoAcid['Xaa'].is_aromatic()
        False
        """
        return self.name in AROMATIC_AA
    
    def is_hydrophobic(self) -> bool:
        """Return True if the amino acid identity is known and is hydrophobic.
        
        Returns
        -------
        bool
            True if the amino acid token represents a hydrophobic amino acid.
            
        Examples
        --------
        >>> AminoAcid['Ala'].is_hydrophobic()
        True
        >>> AminoAcid['Gly'].is_hydrophobic()
        False
        >>> AminoAcid['Xaa'].is_hydrophobic()
        False
        """
        return self.name in HYDROPHOBIC_AA

    # Contains particular elements/functional groups
    
    def has_sulfur(self) -> bool:
        """Return True if the amino acid identity is known and contains sulfur.
        
        Returns
        -------
        bool
            True if the amino acid token represents an amino acid with sulfur.
            
        Examples
        --------
        >>> AminoAcid['Cys'].has_sulfur()
        True
        >>> AminoAcid['Met'].has_sulfur()
        True
        >>> AminoAcid['Ala'].has_sulfur()
        False
        >>> AminoAcid['Xaa'].has_sulfur()
        False
        """
        return self.name in SULFUR_AA
    
    def has_amide(self) -> bool:
        """Return True if the amino acid identity is known and contains an amide.
        
        Returns
        -------
        bool
            True if the amino acid token represents an amino acid with an amide.
            
        Examples
        --------
        >>> AminoAcid['Asn'].has_amide()
        True
        >>> AminoAcid['Gln'].has_amide()
        True
        >>> AminoAcid['Ala'].has_amide()
        False
        >>> AminoAcid['Xaa'].has_amide()
        False
        """
        return self.name in AMIDE_AA
    
    def has_hydroxyl(self) -> bool:
        """Return True if the amino acid identity is known and contains a hydroxyl.
        
        Returns
        -------
        bool
            True if the amino acid token represents an amino acid with a hydroxyl.
            
        Examples
        --------
        >>> AminoAcid['Ser'].has_hydroxyl()
        True
        >>> AminoAcid['Thr'].has_hydroxyl()
        True
        >>> AminoAcid['Tyr'].has_hydroxyl()
        True
        >>> AminoAcid['Ala'].has_hydroxyl()
        False
        >>> AminoAcid['Xaa'].has_hydroxyl()
        False
        """
        return self.name in HYDROXYL_AA


class Codon(Enum):
    AAA = 1
    AAC = 2
    AAG = 3
    AAT = 4
    ACA = 5
    ACC = 6
    ACG = 7
    ACT = 8
    AGA = 9
    AGC = 10
    AGG = 11
    AGT = 12
    ATA = 13
    ATC = 14
    ATG = 15
    ATT = 16
    CAA = 17
    CAC = 18
    CAG = 19
    CAT = 20
    CCA = 21
    CCC = 22
    CCG = 23
    CCT = 24
    CGA = 25
    CGC = 26
    CGG = 27
    CGT = 28
    CTA = 29
    CTC = 30
    CTG = 31
    CTT = 32
    GAA = 33
    GAC = 34
    GAG = 35
    GAT = 36
    GCA = 37
    GCC = 38
    GCG = 39
    GCT = 40
    GGA = 41
    GGC = 42
    GGG = 43
    GGT = 44
    GTA = 45
    GTC = 46
    GTG = 47
    GTT = 48
    TAA = 49  # Stop codon
    TAC = 50
    TAG = 51  # Stop codon
    TAT = 52
    TCA = 53
    TCC = 54
    TCG = 55
    TCT = 56
    TGA = 57  # Stop codon
    TGC = 58
    TGG = 59
    TGT = 60
    TTA = 61
    TTC = 62
    TTG = 63
    TTT = 64
    # Special tokens
    Mask = 65
    Special = 66
    Gap = 0
    # Ambiguous codons
    NNN = 67
    
    def __str__(self) -> str:
        if self.name == 'Gap':
            return '---'
        elif self.name == 'Mask':
            return '###'
        elif self.name == 'Special':
            return '@@@'
        return self.name
    
    def __repr__(self) -> str:
        return self.__str__()
    
    def __eq__(self, other: Any) -> bool:
        if isinstance(other, Codon):
            return self.value == other.value
        elif isinstance(other, str):
            return str(self) == other
        else:
            return False
        
    def __hash__(self):
        return hash(str(self))
    
    @classmethod
    def from_str(cls, s: str) -> 'Codon':
        """Convert a string to a Codon object.
        
        Parameters
        ----------
        s : str
            The string to convert.
        
        Returns
        -------
        Codon
            The codon corresponding to the given string.
            
        Raises
        ------
        ValueError
            If the string is not length 3 or not a valid codon sequence.
            
        Examples
        --------
        >>> Codon.from_str('AAA')
        AAA
        >>> Codon.from_str('TAG')
        TAG
        >>> Codon.from_str('---')
        ---
        >>> Codon.from_str('###')
        ###
        >>> Codon.from_str('AA')
        Traceback (most recent call last):
        ValueError: Invalid string length: 2 (AA)
        >>> Codon.from_str('A-A')
        Traceback (most recent call last):
        ValueError: Invalid codon: A-A
        """
        if len(s) != 3:
            raise ValueError(f'Invalid string length: {len(s)} ({s})')
        try:
            return cls[FROM_TRIPLET_TOKEN[s]]
        except KeyError:
            raise ValueError(f'Invalid codon: {s}')

    @classmethod
    def from_nucleotides(cls, seq: Sequence[Nucleotide]) -> 'Codon':
        """Convert a sequence of nucleotides to a codon.
        
        Parameters
        ----------
        seq : Sequence[Nucleotide]
            The sequence of nucleotides.
            
        Returns
        -------
        Codon
            The codon corresponding to the given sequence of nucleotides.
            
        Raises
        ------
        ValueError
            If the sequence of nucleotides is invalid.
            
        Examples
        --------
        >>> Codon.from_nucleotides([Nucleotide.A, Nucleotide.A, Nucleotide.A])
        AAA
        >>> Codon.from_nucleotides([Nucleotide.T, Nucleotide.A, Nucleotide.G])
        TAG
        >>> Codon.from_nucleotides([Nucleotide.Gap, Nucleotide.Gap, Nucleotide.Gap])
        ---
        >>> Codon.from_nucleotides([Nucleotide.Mask, Nucleotide.Mask, Nucleotide.Mask])
        ###
        >>> Codon.from_nucleotides([Nucleotide.A, Nucleotide.A])
        Traceback (most recent call last):
        ValueError: Invalid nucleotide sequence: [A, A]
        >>> Codon.from_nucleotides([Nucleotide.A, Nucleotide.Gap, Nucleotide.A])
        Traceback (most recent call last):
        ValueError: Invalid nucleotide sequence: [A, -, A]
        """
        try:
            str_triplet = ''.join([str(n) for n in seq])
            return cls.from_str(str_triplet)
        except ValueError:
            raise ValueError(f'Invalid nucleotide sequence: {seq}')

    @classmethod
    def from_onehot(cls, onehot: Sequence[int]) -> 'Codon':
        """Convert a onehot vector to a Codon object.
        
        Parameters
        ----------
        onehot : Sequence[int]
            The one-hot vector.
            
        Returns
        -------
        Codon
            The codon corresponding to the given one-hot vector.
            
        Raises
        ------
        ValueError
            If the one-hot vector is invalid.
        """
        if (len(onehot) != len(cls)) or (sum(onehot) != 1):
            raise ValueError(f'Invalid onehot value: {onehot}')
        try:
            cls(onehot.index(1))
        except ValueError:
            raise ValueError(f'Invalid one-hot value: {onehot}')

    @classmethod
    def from_nucleotide_onehot(cls, onehot: Sequence[Sequence[int]]) -> 'Codon':
        """Convert sequence of Nucleotide onehot vectors to a Codon object.
        
        Parameters
        ----------
        onehot : Sequence[Sequence[int]]
            The one-hot vector.
            
        Returns
        -------
        Codon
            The codon corresponding to the given one-hot vector.
            
        Raises
        ------
        ValueError
            If the one-hot vector is invalid.
        """
        if len(onehot) != 3:
            raise ValueError(f'Invalid number of onehot vectors: {onehot}')
        try:
            return cls.from_nucleotides([Nucleotide.from_onehot(onehot[i]) for i in range(3)])
        except ValueError:
            raise ValueError(f'Invalid one-hot value: {onehot}')

    @classmethod
    def start_codon(cls) -> 'Codon':
        """Return the start codon.
        
        Returns
        -------
        Codon
            The start codon ATG.
            
        Examples
        --------
        >>> Codon.start_codon()
        ATG
        """
        return cls.ATG

    # Fomatting methods
    
    def to_onehot(self) -> Sequence[int]:
        """Convert the codon to a one-hot vector.
        
        Returns
        -------
        Sequence[int]
            The one-hot vector.
        """
        onehot = [0] * len(self.__class__)
        onehot[self.value - 1] = 1
        return onehot
    
    def to_str(self) -> str:
        """Convert the codon to a string triplet.
        
        Returns
        -------
        str
            The string triplet.
            
        Examples
        --------
        >>> Codon.AAA.to_str()
        'AAA'
        >>> Codon.TAG.to_str()
        'TAG'
        >>> Codon.Gap.to_str()
        '---'
        >>> Codon.Mask.to_str()
        '###'
        """
        return str(self)
    
    def to_nucleotides(self) -> Sequence[Nucleotide]:
        """Convert the codon to a sequence of nucleotides.
        
        Returns
        -------
        Sequence[Nucleotide]
            The sequence of nucleotides.
            
        Examples
        --------
        >>> Codon.AAA.to_nucleotides()
        [A, A, A]
        >>> Codon.TAG.to_nucleotides()
        [T, A, G]
        >>> Codon.Gap.to_nucleotides()
        [-, -, -]
        >>> Codon.Mask.to_nucleotides()
        [#, #, #]
        """
        return [Nucleotide.from_str(n) for n in self.to_str()]
    
    # Symbol type
    
    def is_start_codon(self) -> bool:
        return self.name == 'ATG'
    
    def is_stop_codon(self) -> bool:
        return self.name in STOP_CODONS
    
    def is_gap(self) -> bool:
        return self.name == '---'
    
    def is_masked(self) -> bool:
        return self.name == '###'
    
    def is_special(self) -> bool:
        return (
            self.name == 'Mask' or
            self.name == 'Gap' or
            self.name == 'Special'
        )
        
    def is_any(self) -> bool:
        return self.name == 'NNN'
    
    def is_ambiguous(self) -> bool:
        return self.is_any()
    
    def is_degenerate(self) -> bool:
        return self.is_any()
    
    # Properties
    
    def is_twofold_degenerate(self) -> bool:
        return self.name in TWOFOLD_DEGENERATE_CODONS
    
    def is_threefold_degenerate(self) -> bool:
        return self.name in THREEFOLD_DEGENERATE_CODONS
    
    def is_fourfold_degenerate(self) -> bool:
        return self.name in FOURFOLD_DEGENERATE_CODONS

    def is_sixfold_degenerate(self) -> bool:
        return self.name in SIXFOLD_DEGENERATE_CODONS

    def translate(self) -> AminoAcid:
        """Translate the codon to an amino acid.
        
        Returns
        -------
        AminoAcid
            The amino acid corresponding to the codon.
            
        Raises
        ------
        ValueError
            If the codon is not a valid codon.
            
        Examples
        --------
        >>> Codon.ATG.translate()
        M
        >>> Codon.AAA.translate()
        K
        >>> Codon.TAG.translate()
        *
        >>> Codon.Gap.translate()
        -
        """
        if self.is_stop_codon():
            return AminoAcid.Ter
        try:
            return AminoAcid[TRANSLATION_TABLE[self.name]]
        except KeyError:
            raise ValueError(f'Cannot be translated: {self}')
