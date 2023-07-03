from typing import Any
from collections.abc import Sequence
from enum import Enum

from biodatatypes.constants.nucleotide import *
from biodatatypes.constants.aminoacid import *
from biodatatypes.constants.codon import *


class BioToken(Enum):
    """Base class for all biological sequence unit tokens as Enum members.

    This subclasses the Enum datatype but does not define any members.
    This class is intended to be subclassed for each type of unit token: 
    Nucleotide, AminoAcid, Codon. 
    Each subclass should define its own set of tokens as class attributes. 
    
    Attributes
    ----------
    name : str
        The name of the token. This is used to identify the token and may be also be its string representations.
        
        Standard tokens should be named after their string representation.
        Special tokens like those representing a gap, mask, or other are always named 'Gap', 'Mask', and 'Other' respectively.
        
    value : int
        The integer value of the token. This is used for indexing into
        onehot vectors.
        
        Gap tokens are always assigned the value 0.
        Standard tokens are always assigned the values from 1 to n.
        Special tokens are always assigned the the highest values, and the Mask token is always assigned the highest value.
    """
    def __repr__(self) -> str:
        return self.__str__()
    
    def __eq__(self, other: object) -> bool:
        if isinstance(other, BioToken):
            return self.value == other.value
        elif isinstance(other, str):
            return str(self) == other
        else:
            return False
    
    def __hash__(self) -> int:
        return hash(str(self))
    
    def __str__(self) -> str:
        ...
        
    # Factory methods
    
    @classmethod
    def from_str(cls, s: str) -> 'BioToken':
        ...
    
    @classmethod
    def from_int(cls, i: int) -> 'BioToken':
        try:
            return cls(i)
        except ValueError:
            raise ValueError(f'Invalid index value for {cls.__name__}: {i}')
    
    @classmethod
    def from_onehot(cls, onehot:Sequence[int]) -> 'BioToken':
        if (len(onehot) != len(cls)) or (sum(onehot) != 1):
            raise ValueError(f'Invalid onehot vector for {cls.__name__}: {onehot}')
        try:
            return cls(onehot.index(1))
        except ValueError:
            raise ValueError(f'Invalid onehot vector for {cls.__name__}: {onehot}')
        
    # Conversion methods
    
    def to_onehot(self) -> Sequence[int]:
        onehot = [0] * len(self.__class__)
        onehot[self.value] = 1
        return onehot
        
    def to_int(self) -> int:
        return self.value
        
    # Token type methods
    
    def is_standard(self) -> bool:
        ...
        
    def is_mask(self) -> bool:
        return self == self.Mask
        
    def is_gap(self) -> bool:
        return self == self.Gap
    
    def is_other(self) -> bool:
        return self == self.Other
        
    def is_special(self) -> bool:
        return self.is_mask() or self.is_gap() or self.is_other()
        
    def is_any(self) -> bool:
        ...
        
    def is_ambiguous(self) -> bool:
        ...


class Nucleotide(BioToken):
    """A nucleotide token.
    
    Attributes
    ----------
    name : str
        The name of the token. This is used to identify the token and may be also be its string representations.
    value : int
        The integer value of the token. This is used for indexing into onehot vectors.
    """
    A = 1
    C = 2
    G = 3
    T = 4
    # Non-standard symbols
    U = 5
    # Degenerate symbols
    R = 6  # A or G, purine
    Y = 7 # C or T, pyrimidine
    S = 8  # G or C, strong
    W = 9  # A or T, weak
    K = 10  # G or T, keto
    M = 11  # A or C, amino
    B = 12  # C or G or T, not A
    D = 13  # A or G or T, not C
    H = 14  # A or C or T, not G
    V = 15  # A or C or G, not T
    N = 16  # A or C or G or T, any
    # Special tokens
    Gap = 0  # Gap, -
    Other = 17  # Other unspecified but valid character, @
    Mask = 18  # Masked, #

    def __str__(self) -> str:
        return TO_NUCLEOTIDE_ONE_LETTER_TOKEN[self.name]

    # Factory methods
    
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
    def from_int(cls, i: int) -> 'Nucleotide':
        """Convert an integer to a Nucleotide.
        
        Parameters
        ----------
        i : int
            An integer representing a nucleotide.
            
        Returns
        -------
        Nucleotide
            A Nucleotide object.
            
        Raises
        ------
        ValueError
            If the integer does not represent a nucleotide.
            
        Examples
        --------
        >>> Nucleotide.from_int(1)
        A
        >>> Nucleotide.from_int(0)
        -
        >>> Nucleotide.from_int(18)
        #
        >>> Nucleotide.from_int(19)
        Traceback (most recent call last):
        ValueError: Invalid index value for Nucleotide: 19
        """
        return super().from_int(i)

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
        #
        >>> Nucleotide.from_onehot([0, 0, 0, 0, 0, 0, 0, 1])
        Traceback (most recent call last):
        ValueError: Invalid onehot vector for Nucleotide: [0, 0, 0, 0, 0, 0, 0, 1]
        >>> Nucleotide.from_onehot([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        Traceback (most recent call last):
        ValueError: Invalid onehot vector for Nucleotide: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        >>> Nucleotide.from_onehot([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        Traceback (most recent call last):
        ValueError: Invalid onehot vector for Nucleotide: [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        """
        return super().from_onehot(onehot)

    # Conversion methods
    
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
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
        """
        return super().to_onehot()
    
    def to_int(self) -> int:
        """Convert a Nucleotide to an integer.
        
        Returns
        -------
        int
            An integer representing the nucleotide.
            
        Examples
        --------
        >>> Nucleotide.A.to_int()
        1
        >>> Nucleotide.Gap.to_int()
        0
        >>> Nucleotide.N.to_int()
        16
        """
        return super().to_int()
    
    # Nucleotide-only conversions
    
    def to_complement(self) -> 'Nucleotide':
        """Convert a Nucleotide to its complement.
        
        Returns
        -------
        Nucleotide
            The complement of the nucleotide.
            
        Examples
        --------
        >>> Nucleotide.A.to_complement()
        T
        >>> Nucleotide.Gap.to_complement()
        -
        >>> Nucleotide.N.to_complement()
        N
        """
        complement_str = COMPLEMENTARY_NUCLEOTIDES[self.name]
        return Nucleotide[complement_str]
    
    # Token type methods

    def is_standard(self):
        """Return True if the nucleotide is a standard nucleotide.

        Returns
        -------
        bool
            True if the nucleotide is a standard nucleotide, False otherwise.
            
        See Also
        --------
        is_special : Checks if it is a special token
        
        Examples
        --------
        >>> Nucleotide.A.is_standard()
        True
        >>> Nucleotide.Gap.is_standard()
        False
        """
        return self.name in ['A', 'C', 'G', 'T']

    def is_mask(self):
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
        >>> Nucleotide.Mask.is_mask()
        True
        >>> Nucleotide.A.is_mask()
        False
        >>> Nucleotide.Gap.is_mask()
        False
        """
        return super().is_mask()
    
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
        >>> Nucleotide.Gap.is_gap()
        True
        >>> Nucleotide.A.is_gap()
        False
        >>> Nucleotide.Mask.is_gap()
        False
        """
        return super().is_gap()
    
    def is_other(self):
        """Return True if the nucleotide is a valid character but is not specified
        in the enumeration.
        
        Returns
        -------
        bool
            True if the nucleotide is an unspecified but valid character, False otherwise.
            
        See Also
        --------
        is_special : Checks if it is a special token
        
        Examples
        --------
        >>> Nucleotide.Other.is_other()
        True
        >>> Nucleotide.A.is_other()
        False
        >>> Nucleotide.N.is_other()
        False
        >>> Nucleotide.Gap.is_other()
        False
        """
        return super().is_other()

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
        >>> gap = Nucleotide.Gap
        >>> gap.is_special()
        True
        >>> mask = Nucleotide.Mask
        >>> mask.is_special()
        True
        >>> other = Nucleotide.Other
        >>> other.is_special()
        True
        >>> a = Nucleotide.A
        >>> a.is_special()
        False
        >>> n = Nucleotide.N
        >>> n.is_special()
        False
        """
        return super().is_special()
    
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


class AminoAcid(BioToken):
    """Amino acid tokens.
    
    Attributes
    ----------
    name : str
        The name of the amino acid in three-letter code.
        
        Other unspecified amino acids are named 'Xaa'.
        
        Special tokens representing a gap, mask, other special non-amino acid token, 
        or termination signal are named 'Gap', 'Mask', 'Other', and 'Stop', respectively.
    value : int
        The integer value of the amino acid.
        
        Standard amino acids are assigned values from 1 to 20.
        Sec and Pyl are assigned values 21 and 22, respectively.
        Ambiguous amino acids are assigned values 23 to 26.
        Special non-amino acid tokens are assigned values 0, and 27 to 29 inclusive.
    """
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
    Gap = 0  # Gap, -
    Stop = 27 # Terminator, *
    Other = 28  # Other unspecified but valid character, @
    Mask = 29  # Masked, ?
    
    def __str__(self) -> str:
        return NAME_TO_AMINO_ACID_ONE_LETTER_TOKEN[self.value]
    
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
            raise ValueError(f'Invalid amino acid string: {s}')
    
    @classmethod
    def from_int(cls, i: int) -> 'AminoAcid':
        """Convert an integer to an AminoAcid object.
        
        Parameters
        ----------
        i : int
            An integer representing an amino acid.
            
        Returns
        -------
        AminoAcid
            An AminoAcid object.
            
        Raises
        ------
        ValueError
            If the integer is not a valid amino acid.
            
        Examples
        --------
        >>> AminoAcid.from_int(1)
        A
        >>> AminoAcid.from_int(21)
        U
        >>> AminoAcid.from_int(0)
        -
        >>> AminoAcid.from_int(29)
        #
        >>> AminoAcid.from_int(30)
        Traceback (most recent call last):
        ValueError: Invalid index value for AminoAcid: 30
        """
        return super().from_int(i)
    
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
        #
        >>> AminoAcid.from_onehot([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        Traceback (most recent call last):
        ValueError: Invalid onehot vector for AminoAcid: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        """
        return super().from_onehot(onehot)
    
    # Factory methods unique to AminoAcid
    
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
            return cls[AMINO_ACID_ONE_LETTER_TOKEN_TO_NAME[one_letter_code.upper()]]
        except KeyError:
            raise ValueError(f'Invalid one-letter code: {one_letter_code}')
    
    @classmethod
    def from_three_letter(cls, three_letter_code: str) -> 'AminoAcid':
        f"""Convert a three-letter code to an AminoAcid object.
        
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
        ValueError: Invalid three-letter code for Amino Acid: XXX
        """
        try:
            return cls[AMINO_ACID_THREE_LETTER_TOKEN_TO_NAME[three_letter_code.capitalize()]]
        except KeyError:
            raise ValueError(f'Invalid three-letter code for {cls.__name__}: {three_letter_code}')
    
    # Conversion methods
    
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
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
        """
        return super().to_onehot()
    
    def to_int(self) -> int:
        """Convert an AminoAcid object to an integer representation.
        
        Returns
        -------
        int
            An integer representing the amino acid.
            
        Examples
        --------
        >>> aa = AminoAcid.Ala
        >>> aa.to_int()
        1
        >>> gap = AminoAcid.Gap
        >>> gap.to_int()
        0
        >>> unknown_aa = AminoAcid.Xaa
        >>> unknown_aa.to_int()
        26
        >>> masked_aa = AminoAcid.Mask
        >>> masked_aa.to_int()
        29
        """
        return super().to_int()
    
    # Conversion methods unique to AminoAcid
    
    def to_one_letter(self) -> str:
        """Convert an AminoAcid object to a one-letter code.
        
        Returns
        -------
        str
            A one-letter code representing the amino acid.
        
        Examples
        --------
        >>> aa = AminoAcid.from_str('A')
        >>> aa.to_one_letter()
        'A'
        >>> gap = AminoAcid.from_str('-')
        >>> gap.to_one_letter()
        '-'
        >>> unknown = AminoAcid.from_str('X')
        >>> unknown.to_one_letter()
        'X'
        """
        return str(self)
    
    def to_three_letter(self) -> str:
        """Convert an AminoAcid object to a three-letter code.
        
        Returns
        -------
        str
            A three-letter code representing the amino acid.
            
        Examples
        --------
        >>> aa = AminoAcid.from_str('A')
        >>> aa.to_three_letter()
        'Ala'
        >>> gap = AminoAcid.from_str('-')
        >>> gap.to_three_letter()
        'Gap'
        >>> unknown = AminoAcid.from_str('X')
        >>> unknown.to_three_letter()
        'Xaa'
        """
        return NAME_TO_AMINO_ACID_THREE_LETTER_TOKEN[self.name]
    
    # Token type methods
    
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
    
    def is_mask(self) -> bool:
        """Return True if masked.
        
        Returns
        -------
        bool
            True if the amino acid object is masked.
        
        Examples
        --------
        >>> AminoAcid['Mask'].is_mask()
        True
        >>> AminoAcid['Ala'].is_mask()
        False
        >>> AminoAcid['Gap'].is_mask()
        False
        """
        return super().is_mask()
    
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
        return super().is_gap()
    
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
        >>> AminoAcid['Stop'].is_special()
        True
        >>> AminoAcid['Other'].is_special()
        True
        >>> AminoAcid['Ala'].is_special()
        False
        >>> AminoAcid['Xaa'].is_special()
        False
        >>> AminoAcid['Sec'].is_special()
        False
        """
        return super().is_special() or self.is_stop()
    
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
        return self == self.Xaa
    
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
            self.is_any() or 
            self == self.Asx or 
            self == self.Glx or 
            self == self.Xle)
        
    # Token type methods unique to AminoAcid
    
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
        return self.name == 'Sec' or self.name == 'Pyl'
        
    def is_stop(self) -> bool:
        """Return True if it represents a terminator sequence.
        
        Returns
        -------
        bool
            True if the amino acid token is a terminator sequence.
            
        Examples
        --------
        >>> AminoAcid['Stop'].is_stop()
        True
        >>> AminoAcid['Ala'].is_stop()
        False
        >>> AminoAcid['Gap'].is_stop()
        False
        """
        return self == self.Stop

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


class Codon(BioToken):
    """A codon token.
    
    Attributes
    ----------
    name : str
        Three-letter sequence for the codon (e.g. 'AAA') or the designated name 
        of the special token.
    value : int
        The integer value of the codon.
        
        Standard codons are assigned values 1-64.
        Ambiguous codon NNN is assigned value 65.
        Special tokens 'Gap', 'Other', and 'Mask' are assigned values 0, 66, and 67, respectively.
    """
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
    # Ambiguous codons
    NNN = 65
    # Special tokens
    Gap = 0
    Other = 66
    Mask = 67
    
    def __str__(self) -> str:
        return NAME_TO_TRIPLET_TOKEN[self.name]
    
    # Factory methods
    
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
        ValueError: Invalid codon string: A-A
        """
        if len(s) != 3:
            raise ValueError(f'Invalid string length: {len(s)} ({s})')
        try:
            return cls[TRIPLET_TOKEN_TO_NAME[s.upper()]]
        except KeyError:
            try:
                return cls[TRIPLET_TOKEN_TO_NAME[s.capitalize()]]
            except KeyError:
                raise ValueError(f'Invalid codon string: {s}')

    @classmethod
    def from_int(cls, i: int) -> 'Codon':
        """Convert an integer to a Codon object.
        
        Parameters
        ----------
        i : int
            The integer to convert.
            
        Returns
        -------
        Codon
            The codon corresponding to the given integer.
            
        Raises
        ------
        ValueError
            If the integer is not a valid codon index.
        """
        return super().from_int(i)

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
        return super().from_onehot(onehot)

    # Factory methods unique to Codon
        
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
            raise ValueError(f'Invalid one-hot value for {cls.__name__}: {onehot}')

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

    # Conversion methods
    
    def to_onehot(self) -> Sequence[int]:
        """Convert the codon to a one-hot vector.
        
        Returns
        -------
        Sequence[int]
            The one-hot vector.
        """
        return super().to_onehot()
    
    def to_int(self) -> int:
        """Convert the codon to an integer representation.
        
        Returns
        -------
        int
            The integer.
        """
        return super().to_int()
    
    # Conversion methods unique to Codon
    
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
        return [Nucleotide.from_str(n) for n in str(self)]
    
    def to_anticodon(self) -> 'Codon':
        """Return the anticodon of the codon.
        
        Returns
        -------
        Codon
            Complimentary sequence of the codon.
            
        Examples
        --------
        >>> Codon.AAA.to_anticodon()
        TTT
        >>> Codon.TAG.to_anticodon()  # stop codon
        ATC
        >>> Codon.Gap.to_anticodon()
        ---
        >>> Codon.Mask.to_anticodon()
        ###
        """
        if self.is_special():
            return self.__class__(self.value)
        anticodn_str = ''.join([COMPLEMENTARY_NUCLEOTIDES[n] for n in str(self)])
        return self.__class__[anticodn_str]
        
    # Token type methods
    
    def is_standard(self) -> bool:
        """Return True if the codon is one of the 64 standard codon.
        
        Returns
        -------
        bool
            True if the codon is a standard codon.
            
        Examples
        --------
        >>> Codon.AAA.is_standard()
        True
        >>> Codon.ATG.is_standard()  # start codon
        True
        >>> Codon.TAG.is_standard()  # Stop codon
        True
        >>> Codon.Gap.is_standard()
        False
        >>> Codon.Mask.is_standard()
        False
        """
        return self.name in STANDARD_CODONS

    def is_mask(self) -> bool:
        """Return True if the codon is a mask codon.
        
        Returns
        -------
        bool
            True if the codon is a mask codon, False otherwise.
        
        Examples
        --------
        >>> Codon.Mask.is_mask()
        True
        >>> Codon.AAA.is_mask()
        False
        >>> Codon.ATG.is_mask()
        False
        >>> Codon.TAG.is_mask()
        False
        >>> Codon.Gap.is_mask()
        False
        """
        return super().is_mask()

    def is_gap(self) -> bool:
        """Return True if the codon is a gap codon.
        
        Returns
        -------
        bool
            True if the codon is a gap codon, False otherwise.
            
        Examples
        --------
        >>> Codon.Gap.is_gap()
        True
        >>> Codon.AAA.is_gap()
        False
        >>> Codon.ATG.is_gap()
        False
        >>> Codon.TAG.is_gap()
        False
        >>> Codon.Mask.is_gap()
        False
        """
        return super().is_gap()
    
    def is_other(self) -> bool:
        """Return True if the codon is a valid but unspecified codon.
        
        Returns
        -------
        bool
            True if the codon is a valid but unspecified codon, False otherwise.
            
        Examples
        --------
        >>> Codon.Other.is_other()
        True
        >>> Codon.AAA.is_other()
        False
        >>> Codon.ATG.is_other()
        False
        >>> Codon.TAG.is_other()
        False
        >>> Codon.Gap.is_other()
        False
        """
        return super().is_other()
    
    def is_special(self) -> bool:
        """Return True if the codon is a mask, gap, or other codon.
        
        Returns
        -------
        bool
            True if the codon is a mask, gap, or other codon, False otherwise.
            
        Examples
        --------
        >>> Codon.Mask.is_special()
        True
        >>> Codon.Gap.is_special()
        True
        >>> Codon.Other.is_special()
        True
        >>> Codon.AAA.is_special()
        False
        >>> Codon.ATG.is_special()
        False
        >>> Codon.TAG.is_special()
        False
        """
        return super().is_special()
    
    def is_any(self) -> bool:
        """Return True if the codon is NNN.
        
        Returns
        -------
        bool
            True if the codon is NNN, False otherwise.
            
        Examples
        --------
        >>> Codon.NNN.is_any()
        True
        >>> Codon.Gap.is_any()
        False
        >>> Codon.Mask.is_any()
        False
        >>> Codon.AAA.is_any()
        False
        """
        return self.name == 'NNN'
    
    def is_ambiguous(self) -> bool:
        """Return True if the codon is ambiguous.
        
        Returns
        -------
        bool
            True if the codon is ambiguous, False otherwise.
            
        Examples
        --------
        >>> Codon.NNN.is_ambiguous()
        True
        >>> Codon.Gap.is_ambiguous()
        False
        >>> Codon.Mask.is_ambiguous()
        False
        >>> Codon.AAA.is_ambiguous()
        False
        """
        return self.is_any()
    
    # Codon-specific token type methods
        
    def is_start_codon(self) -> bool:
        """Return True if the codon is the start codon.
        
        Returns
        -------
        bool
            True if the codon is the start codon, False otherwise.
            
        Examples
        --------
        >>> Codon.ATG.is_start_codon()
        True
        >>> Codon.AAA.is_start_codon()
        False
        >>> Codon.Mask.is_start_codon()
        False
        """
        return self.name == 'ATG'
    
    def is_stop_codon(self) -> bool:
        """Return True if the codon is a stop codon.

        Returns
        -------
        bool
            True if the codon is a stop codon, False otherwise.
        
        Examples
        --------
        >>> Codon.TAG.is_stop_codon()
        True
        >>> Codon.TAA.is_stop_codon()
        True
        >>> Codon.TGA.is_stop_codon()
        True
        >>> Codon.ATG.is_stop_codon()
        False
        >>> Codon.AAA.is_stop_codon()
        False
        >>> Codon.Mask.is_stop_codon()
        False
        """
        return self.name in STOP_CODONS

    # Properties
    
    def is_twofold_degenerate(self) -> bool:
        """Return True if the codon is twofold degenerate.
        
        Returns
        -------
        bool
            True if the codon is twofold degenerate, False otherwise.
        
        Examples
        --------
        >>> Codon.AAA.is_twofold_degenerate()  # other codon is AAG
        True
        >>> Codon.ATG.is_twofold_degenerate()  # Met is has no degenerate sites
        False
        >>> Codon.TAG.is_twofold_degenerate()
        False
        """
        return self.name in TWOFOLD_DEGENERATE_CODONS
    
    def is_threefold_degenerate(self) -> bool:
        """Return True if the codon is threefold degenerate.
        
        Returns
        -------
        bool 
            True if the codon is threefold degenerate, False otherwise.
            
        Examples
        --------
        >>> Codon.ATA.is_threefold_degenerate()  # Ile
        True
        >>> Codon.TAG.is_threefold_degenerate()  # Stop
        False
        >>> Codon.AAA.is_threefold_degenerate()
        False
        """
        return self.name in THREEFOLD_DEGENERATE_CODONS
    
    def is_fourfold_degenerate(self) -> bool:
        """Return True if the codon is fourfold degenerate.
        
        Returns
        -------
        bool
            True if the codon is fourfold degenerate, False otherwise.
            
        Examples
        --------
        >>> Codon.GTA.is_fourfold_degenerate()  # Val
        True
        >>> Codon.TTT.is_fourfold_degenerate()  # Phe
        False
        """
        return self.name in FOURFOLD_DEGENERATE_CODONS

    def is_sixfold_degenerate(self) -> bool:
        """Return True if the codon is sixfold degenerate.
        
        Returns
        -------
        bool
            True if the codon is sixfold degenerate, False otherwise.
            
        Examples
        --------
        >>> Codon.CGA.is_sixfold_degenerate()  # Arg
        True
        >>> Codon.TTT.is_sixfold_degenerate()  # Phe
        False
        """
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
            return AminoAcid.Stop
        try:
            return AminoAcid[TRANSLATION_TABLE[self.name]]
        except KeyError:
            raise ValueError(f'Cannot be translated: {self}')
