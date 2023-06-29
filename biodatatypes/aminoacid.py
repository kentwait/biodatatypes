from typing import Union, List
from collections.abc import Sequence
from enum import Enum

from biodatatypes.sequence import SeqMixin, MaskMixin, GapMixin


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
            return cls[FROM_ONE_LETTER_TOKEN[one_letter_code.upper()]]
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
            return cls[FROM_THREE_LETTER_TOKEN[three_letter_code.capitalize()]]
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
        if len(onehot) != len(cls):
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
        return TO_ONE_LETTER_TOKEN[self.name]
    
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


class AminoAcidSequence(GapMixin, MaskMixin, SeqMixin, Sequence):
    unit = AminoAcid
    
    def __init__(self, 
            sequence: Sequence[AminoAcid],
            is_standard: bool = None, 
            is_degenerate: bool = None, 
            is_gapped: bool = None, 
            is_masked: bool = None):
        self._sequence = list(sequence)
        self._is_standard = is_standard
        self._is_degenerate = is_degenerate
        self._is_gapped = is_gapped
        self._is_masked = is_masked
            
    @classmethod
    def from_str(cls, sequence: str) -> 'AminoAcidSequence':
        """Create an AminoAcidSequence from a string of amino acid tokens.
        
        Parameters
        ----------
        sequence : str
            A string of amino acid tokens in one-letter code.
            
        Returns
        -------
        AminoAcidSequence
            An AminoAcidSequence object.
            
        Examples
        --------
        >>> AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        ARNDCEQGHILKMFPSTWYV
        """
        return cls([AminoAcid.from_str(s) for s in sequence])
    
    @classmethod
    def from_onehot(cls, sequence: Sequence[Sequence[int]]) -> 'AminoAcidSequence':
        """Create an AminoAcidSequence from a one-hot encoded sequence.
        
        Parameters
        ----------
        sequence : Sequence[Sequence[int]]
            A sequence of one-hot encoded amino acid tokens.
        
        Returns
        -------
        AminoAcidSequence
            An AminoAcidSequence object.
        """        
        return cls([AminoAcid.from_onehot(s) for s in sequence])
    
    def to_str(self) -> str:
        """Return a string of amino acid tokens in one-letter code.
        
        Returns
        -------
        str
            A string of amino acid tokens in one-letter code.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.to_str()
        'ARNDCEQGHILKMFPSTWYV'
        >>> gapped_seq = AminoAcidSequence.from_str('ARNDC--EQGHILKMFPSTWYV')
        >>> gapped_seq.to_str()
        'ARNDC--EQGHILKMFPSTWYV'
        """
        return SeqMixin.to_str(self)
        
    def to_onehot(self) -> Sequence[Sequence[int]]:
        return SeqMixin.to_onehot(self)

    def startswith(self, seq: Union[str, Sequence[AminoAcid], 'AminoAcidSequence']) -> bool:
        """Return True if the sequence starts with the given sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[AminoAcid], AminoAcidSequence]
            The sequence to check.
            
        Returns
        -------
        bool
            True if the sequence starts with the given sequence.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.startswith('ARN')
        True
        >>> seq.startswith('A')
        True
        >>> seq.startswith('ARNDCEQGHILKMFPSTWYVX')
        False
        >>> seq.startswith('-')
        False
        """
        if isinstance(seq, str):
            seq = AminoAcidSequence.from_str(seq)
        return SeqMixin.startswith(self, seq)

    def endswith(self, seq: Union[str, Sequence[AminoAcid], 'AminoAcidSequence']) -> bool:
        """Return True if the sequence ends with the given sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[AminoAcid], AminoAcidSequence]
            The sequence to check.
            
        Returns
        -------
        bool
            True if the sequence ends with the given sequence.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.endswith('V')
        True
        >>> seq.endswith('WYV')
        True
        >>> seq.endswith('ARNDC--EQGHILKMFPSTWYV')
        False
        >>> seq.endswith('-')
        False
        """
        if isinstance(seq, str):
            seq = AminoAcidSequence.from_str(seq)
        return SeqMixin.endswith(self, seq)

    def find(self, seq: Union[str, Sequence[AminoAcid], 'AminoAcidSequence']) -> int:
        """Return the index of the first occurrence of the given sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[AminoAcid], AminoAcidSequence]
            The sequence to find.
            
        Returns
        -------
        int
            The index of the first occurrence of the given sequence.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.find('ARN')
        0
        >>> seq.find('RND')
        1
        >>> seq.find('V')
        19
        >>> seq.find('')
        0
        >>> seq.find('X')
        -1
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, AminoAcidSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], AminoAcid):
            seq = ''.join([str(s) for s in seq])
        return SeqMixin.find(self, seq)

    def rfind(self, seq: Union[str, Sequence[AminoAcid], 'AminoAcidSequence']) -> int:
        """Return the index of the last occurrence of the given sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[AminoAcid], AminoAcidSequence]
            The sequence to find.
            
        Returns
        -------
        int
            The index of the last occurrence of the given sequence.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.rfind('ARN')
        0
        >>> seq.rfind('RND')
        1
        >>> seq.rfind('V')
        19
        >>> seq.rfind('')
        20
        >>> seq.rfind('-')
        -1
        >>> seq.rfind('X')
        -1
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, AminoAcidSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], AminoAcid):
            seq = ''.join([str(s) for s in seq])
        return SeqMixin.rfind(self, seq)
    
    def count(self, seq: Union[str, Sequence[AminoAcid], 'AminoAcidSequence']) -> int:
        """Return the number of occurrences of the given sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[AminoAcid], AminoAcidSequence]
            The sequence to count.
            
        Returns
        -------
        int
            The number of occurrences of the given sequence.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.count('ARN')
        1
        >>> seq.count('R')
        1
        >>> seq.count('V')
        1
        >>> seq.count('')
        21
        >>> seq.count('X')
        0
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, AminoAcidSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], AminoAcid):
            seq = ''.join([str(s) for s in seq])
        return SeqMixin.count(self, seq)
    
    def mask(self, positions: Union[int, Sequence[int]]) -> 'AminoAcidSequence':
        """Return a new sequence with the given positions masked.
        
        Parameters
        ----------
        positions : Union[int, Sequence[int]]
            The positions to mask.
            
        Returns
        -------
        AminoAcidSequence
            A new sequence with the given positions masked.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.mask(0)
        #RNDCEQGHILKMFPSTWYV
        >>> seq.mask([0, 1, -1])
        ##NDCEQGHILKMFPSTWY#
        """
        return MaskMixin.mask(self, positions)
            
    def masked_positions(self) -> List[int]:
        """Return the positions that are masked.
        
        Returns
        -------
        List[int]
            The positions that are masked.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.masked_positions()
        []
        >>> seq.mask(0).masked_positions()
        [0]
        >>> seq.mask([0, 1, -1]).masked_positions()
        [0, 1, 19]
        """
        return MaskMixin.masked_positions(self)
    
    def count_masked(self) -> int:
        """Return the number of masked positions.
        
        Returns
        -------
        int
            The number of masked positions.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.count_masked()
        0
        >>> seq.mask(0).count_masked()
        1
        >>> seq.mask([0, 1, -1]).count_masked()
        3
        """
        return MaskMixin.count_masked(self)
    
    def gapped_positions(self) -> List[int]:
        """Return the positions that are gapped.
        
        Returns
        -------
        List[int]
            The positions that are gapped.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.gapped_positions()
        []
        >>> gapped_seq = AminoAcidSequence.from_str('AR-NDCEQGHILKMFPST-W-YV')
        >>> gapped_seq.gapped_positions()
        [2, 18, 20]
        """
        return GapMixin.gapped_positions(self)
    
    def count_gaps(self) -> int:
        """Return the number of gapped positions.
        
        Returns
        -------
        int
            The number of gapped positions.
            
        Examples
        --------
        >>> seq = AminoAcidSequence.from_str('ARNDCEQGHILKMFPSTWYV')
        >>> seq.count_gaps()
        0
        >>> gapped_seq = AminoAcidSequence.from_str('AR-NDCEQGHILKMFPST-W-YV')
        >>> gapped_seq.count_gaps()
        3
        """
        return GapMixin.count_gaps(self)

    
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
