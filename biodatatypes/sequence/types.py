from typing import Union, Optional, List, Any, Iterator, Iterable, cast, Type
from collections.abc import Sequence

from biodatatypes.unit.types import Nucleotide, AminoAcid, Codon
from biodatatypes.unit.base import BioToken, NucleotideEnum, AminoAcidEnum, CodonEnum
from biodatatypes.constants.nucleotide import COMPLEMENTARY_NUCLEOTIDES


class BioSequence(Sequence):
    """
    A BioSequence is a generic sequence composed objects inheriting from the BioToken class
    such as Nucleotide, AminoAcid, or Codon. It is a wrapper around a list of BioTokens
    that provides some useful methods for working with sequences.
    
    This class is meant to be inherited from and extended to provide more specific
    functionality for different types of sequences.
    
    Attributes
    ----------
    sequence : Sequence
        A list of BioTokens.
    is_standard : bool
        Whether or not the sequence is composed of standard tokens.
    is_degenerate : bool
        Whether or not the sequence contains any ambiguous tokens.
    is_gapped : bool
        Whether or not the sequence contains any gap tokens.
    is_masked : bool
        Whether or not the sequence contains any mask tokens.
    
    """    
    def __init__(self, 
            sequence: Sequence[BioToken], 
            is_standard: Optional[bool] = None,
            is_ambiguous: Optional[bool] = None, 
            is_gapped: Optional[bool] = None, 
            is_masked: Optional[bool] = None):
        self.unit = type(sequence[0])
        self._sequence = list(sequence)
        self._is_standard = is_standard
        self._is_degenerate = is_ambiguous
        self._is_gapped = is_gapped
        self._is_masked = is_masked
    
    def __getitem__(self, index: int) -> Any:
        return self.sequence[index]
    
    def __len__(self) -> int:
        return len(self.sequence)
    
    def __str__(self) -> str:
        return ''.join(list(map(str, self.sequence)))
    
    def __repr__(self) -> str:
        return self.__str__()
    
    # def __eq__(self, other: Union['BioSequence', str]) -> bool:
    #     if isinstance(other, str):
    #         return str(self) == other
    #     elif isinstance(other, BioSequence):
    #         return str(self) == str(other)
    #     else:
    #         return False

    @property
    def is_standard(self) -> bool:
        if self._is_standard is None:
            self._is_standard = all([s.is_standard() for s in self.sequence])
        return self._is_standard
    
    @property
    def is_ambiguous(self) -> bool:
        if self._is_degenerate is None:
            self._is_degenerate = any([s.is_ambiguous() for s in self.sequence])
        return self._is_degenerate
    
    @property
    def is_gapped(self) -> bool:
        if self._is_gapped is None:
            self._is_gapped = any([s.is_gap() for s in self.sequence])  # type: ignore
        return self._is_gapped
    
    @property
    def is_masked(self) -> bool:
        if self._is_masked is None:
            self._is_masked = any([s.is_mask() for s in self.sequence])  # type: ignore
        return self._is_masked
    
    @property
    def sequence(self) -> Sequence[BioToken]:
        return self._sequence
    
    @classmethod
    def from_str(cls, 
            sequence: Iterable[str],
            unit: Type[BioToken] = BioToken) -> 'BioSequence':
        """Create an BioSequence from a string of tokens.
        
        Parameters
        ----------
        sequence : Iterable[str]
            A string or iterable of tokens.
            
        Returns
        -------
        NucleotideSequence
            An BioSequence object.
        """
        return cls([unit.from_str(s) for s in sequence])
    
    @classmethod
    def from_onehot(cls, 
            sequence: Sequence[Sequence[int]], 
            unit: Type[BioToken] = BioToken) -> 'BioSequence':
        """Create an NucleotideSequence from a one-hot encoded sequence.
        
        Parameters
        ----------
        sequence : Sequence[Sequence[int]]
            A sequence of one-hot encoded tokens.
        
        Returns
        -------
        BioSequence
            A BioSequence object.
        """        
        return cls([unit.from_onehot(s) for s in sequence])
    
    def to_str(self) -> str:
        """Return a string of amino acid tokens in one-letter code.
        
        Returns
        -------
        str
            A string of amino acid tokens in one-letter code.
        """
        return str(self)
        
    def to_onehot(self) -> Sequence[Sequence[int]]:
        return [s.to_onehot() for s in self.sequence]

    def startswith(self, seq: Sequence) -> bool:
        """Return True if the sequence starts with the given sequence.
        
        Parameters
        ----------
        substr : Sequence
            The sequence to check.
            
        Returns
        -------
        bool
            True if the sequence starts with the given sequence.
        """
        return self.sequence[:len(seq)] == seq

    def endswith(self, seq: Sequence) -> bool:
        """Return True if the sequence ends with the given sequence.
        
        Parameters
        ----------
        other : Sequence
            The sequence to check.
            
        Returns
        -------
        bool
            True if the sequence ends with the given sequence.
        """
        return self.sequence[-len(seq):] == seq

    def find(self, substr: str) -> int:
        """Return the index of the first occurrence of the given sequence.
        
        Parameters
        ----------
        substr : str
            The sequence to find.
            
        Returns
        -------
        int
            The index of the first occurrence of the given sequence.
        """
        if len(substr) == 0:
            return 0
        elif len(substr) > len(self):
            return -1
        return str(self).find(substr)

    def rfind(self, substr: str) -> int:
        """Return the index of the last occurrence of the given sequence.
        
        Parameters
        ----------
        substr : str
            The sequence to find.
            
        Returns
        -------
        int
            The index of the last occurrence of the given sequence.
        """
        if len(substr) == 0:
            return len(self)
        elif len(substr) > len(self):
            return -1
        return str(self).rfind(substr)
    
    def count(self, substr: str) -> int:
        """Return the number of occurrences of the given sequence.
        
        Parameters
        ----------
        substr : str
            The sequence to count.
            
        Returns
        -------
        int
            The number of occurrences of the given sequence.
        """
        if len(substr) == 0:
            return len(self) + 1
        elif len(substr) > len(self):
            return 0
        return str(self).count(substr)

    def mask(self, positions: Union[int, Sequence[int]]) -> Sequence:
        """Return a new sequence with the given positions masked.
        
        Parameters
        ----------
        positions : Union[int, Sequence[int]]
            The positions to mask.
            
        Returns
        -------
        Sequence
            A new sequence with the given positions masked.
        """
        if isinstance(positions, int):
            positions = [positions]
        sequence = list(self._sequence)
        for i in positions:
            sequence[i] = self.unit['Mask']
        return self.__class__(sequence, is_masked=(len(positions) > 0))
    
    def masked_positions(self) -> List[int]:
        """Return the positions that are masked.
        
        Returns
        -------
        List[int]
            The positions that are masked.
        """
        masked_pos = [i for i, n in enumerate(self._sequence) 
                      if n.name == 'Mask']
        if self._is_masked:
            self._is_masked = len(masked_pos) > 0
        return masked_pos

    def count_masked(self) -> int:
        """Return the number of masked positions.
        
        Returns
        -------
        int
            The number of masked positions.
        """
        return len(self.masked_positions())

    def gapped_positions(self) -> List[int]:
        """Return the positions that are gapped.
        
        Returns
        -------
        List[int]
            The positions that are gapped.
        """
        gapped_pos = [i for i, n in enumerate(self._sequence) 
                      if n.name == 'Gap']
        if self._is_gapped is None:
            self._is_gapped = len(gapped_pos) > 0
        return gapped_pos
    
    def count_gaps(self) -> int:
        """Return the number of gapped positions.
        
        Returns
        -------
        int
            The number of gapped positions.
        """
        return len(self.gapped_positions())


class NucleotideSequence(BioSequence):
    """
    A sequence of nucleotides.
    
    Parameters
    ----------
    sequence : Sequence[NucleotideEnum]
        A sequence of nucleotides as Nucleotide objects.
    is_standard : bool, optional
        Whether the sequence is standard (i.e. only contains standard nucleotides).
    is_ambiguous : bool, optional
        True if the sequence contains ambiguous nucleotides
    is_gapped : bool, optional
        True if the sequence is gapped (i.e. contains at least one gap).
    is_masked : bool, optional
        True if the sequence is masked (i.e. contains at least one the masked position).
    """    
    def __init__(self, 
            sequence: Sequence[NucleotideEnum], 
            is_standard: Optional[bool] = None,
            is_ambiguous: Optional[bool] = None, 
            is_gapped: Optional[bool] = None, 
            is_masked: Optional[bool] = None):
        super().__init__(sequence, is_standard, is_ambiguous, is_gapped, is_masked)
        
    @classmethod
    def from_str(cls, sequence: str, unit: Type = Nucleotide) -> 'NucleotideSequence':
        """Create an NucleotideSequence from a string of nucleotide tokens.
        
        Parameters
        ----------
        sequence : str
            A string of nucleotide tokens in one-letter code.
            
        Returns
        -------
        NucleotideSequence
            An NucleotideSequence object.
            
        Examples
        --------
        >>> NucleotideSequence.from_str('ATGCCGTATGAATGA')
        ATGCCGTATGAATGA
        >>> NucleotideSequence.from_str('ATG-A-CCGTATGAA---TGA')
        ATG-A-CCGTATGAA---TGA
        """
        return cast(NucleotideSequence, super().from_str(sequence, unit))
    
    @classmethod
    def from_onehot(cls, sequence: Sequence[Sequence[int]], unit: Type = Nucleotide) -> 'NucleotideSequence':
        """Create an NucleotideSequence from a one-hot encoded sequence.
        
        Parameters
        ----------
        sequence : Sequence[Sequence[int]]
            A sequence of one-hot encoded nucleotide tokens.
        
        Returns
        -------
        NucleotideSequence
            An NucleotideSequence object.
        """        
        return cast(NucleotideSequence, super().from_onehot(sequence, unit))
    
    def to_str(self) -> str:
        """Return a string of nucleotide tokens in one-letter code.
        
        Returns
        -------
        str
            A string of nucleotide tokens in one-letter code.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.to_str()
        'ATGCCGTATGAATGA'
        >>> gapped_seq = NucleotideSequence.from_str('ATG-A-CCGTATGAA---TGA')
        >>> gapped_seq.to_str()
        'ATG-A-CCGTATGAA---TGA'
        """
        return super().to_str()
        
    def to_onehot(self) -> Sequence[Sequence[int]]:
        return super().to_onehot()
    
    def startswith(self, seq: Union[str, Sequence[NucleotideEnum], 'NucleotideSequence']) -> bool:
        """Return True if the sequence starts with the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[NucleotideEnum]]
            A nucleotide or sequence of nucleotides.
        
        Returns
        -------
        bool
            True if the sequence starts with the given nucleotide or sequence, False otherwise.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.startswith('A')
        True
        >>> seq.startswith([Nucleotide.A])
        True
        >>> other_seq = NucleotideSequence.from_str('A')
        >>> seq.startswith(other_seq)
        True
        >>> seq.startswith('ATG')
        True
        >>> seq.startswith([Nucleotide.A, Nucleotide.T, Nucleotide.G])
        True
        >>> seq.startswith('-')
        False
        >>> seq.startswith('A--')
        False
        """
        if isinstance(seq, str):
            seq = cast(Sequence[NucleotideEnum], self.from_str(seq, self.unit).sequence)
        elif isinstance(seq, self.__class__):
            seq = cast(Sequence[NucleotideEnum], seq.sequence)
        return super().startswith(seq)

    def endswith(self, seq: Union[str, Sequence[NucleotideEnum], 'NucleotideSequence']) -> bool:
        """Return True if the sequence ends with the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[NucleotideEnum]]
            A nucleotide or sequence of nucleotides.
        
        Returns
        -------
        bool
            True if the sequence ends with the given nucleotide or sequence, False otherwise.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.endswith('A')
        True
        >>> seq.endswith([Nucleotide.A])
        True
        >>> other_seq = NucleotideSequence.from_str('A')
        >>> seq.endswith(other_seq)
        True
        >>> seq.endswith('TGA')
        True
        >>> seq.endswith('-')
        False
        >>> seq.endswith('A--')
        False
        """
        if isinstance(seq, str):
            seq = cast(Sequence[NucleotideEnum], self.from_str(seq, self.unit).sequence)
        elif isinstance(seq, self.__class__):
            seq = cast(Sequence[NucleotideEnum], seq.sequence)
        return super().endswith(seq)

    def find(self, seq: Union[str, Sequence[NucleotideEnum], 'NucleotideSequence']) -> int:
        """Return the index of the first occurrence of the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[NucleotideEnum]]
            A nucleotide or sequence of nucleotides.
        
        Returns
        -------
        int
            The index of the first occurrence of the given nucleotide or sequence.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.find('A')
        0
        >>> seq.find('TGA')
        8
        >>> seq.find('')
        0
        >>> seq.find('-')
        -1
        >>> seq.find('A--')
        -1
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, NucleotideSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], NucleotideEnum):
            seq = ''.join([str(n) for n in seq])
        return super().find(cast(str, seq))

    def rfind(self, seq: Union[str, Sequence[NucleotideEnum], 'NucleotideSequence']) -> int:
        """Return the index of the last occurrence of the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[NucleotideEnum]]
            A nucleotide or sequence of nucleotides.
        
        Returns
        -------
        int
            The index of the last occurrence of the given nucleotide or sequence.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.rfind('A')
        14
        >>> seq.rfind('TGA')
        12
        >>> seq.rfind('')
        15
        >>> seq.rfind('-')
        -1
        >>> seq.rfind('A--')
        -1
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, NucleotideSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], NucleotideEnum):
            seq = ''.join([str(n) for n in seq])
        return super().rfind(cast(str, seq))

    def count(self, seq: Union[str, Sequence[NucleotideEnum], 'NucleotideSequence']) -> int:
        """Return the number of occurrences of the given nucleotide or sequence.
        
        Parameters
        ----------
        nucl : Union[str, Sequence[NucleotideEnum]]
            A nucleotide or sequence of nucleotides.
        
        Returns
        -------
        int
            The number of occurrences of the given nucleotide or sequence.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.count('A')
        5
        >>> seq.count('TGA')
        2
        >>> seq.count('')
        16
        >>> seq.count('-')
        0
        >>> seq.count('A--')
        0
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, NucleotideSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], NucleotideEnum):
            seq = ''.join([str(n) for n in seq])
        return super().count(cast(str, seq))

    def mask(self, positions: Union[int, Sequence[int]]) -> 'NucleotideSequence':
        """Return a copy of the sequence with the given positions masked.
        
        Parameters
        ----------
        positions : Union[int, Sequence[int]]
            The positions to mask.
            
        Returns
        -------
        NucleotideSequence
            A copy of the sequence with the given positions masked.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.mask(3)
        ATG#CGTATGAATGA
        """
        return cast(NucleotideSequence, super().mask(positions))
    
    def masked_positions(self) -> List[int]:
        """Return the positions of the masked nucleotides in the sequence.
        
        Returns
        -------
        List[int]
            The positions of the masked nucleotides in the sequence.
            
        Examples
        --------
        >>> masked_seq = NucleotideSequence.from_str('ATG#CGTATGAATGA')
        >>> masked_seq.masked_positions()
        [3]
        """
        return super().masked_positions()
        
    def count_masked(self) -> int:
        """Return the number of masked nucleotides in the sequence.
        
        Returns
        -------
        int
            The number of masked nucleotides in the sequence.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.count_masked()
        0
        >>> masked_seq = NucleotideSequence.from_str('ATG#CGTATGAATGA')
        >>> masked_seq.count_masked()
        1
        """
        return super().count_masked()

    def gapped_positions(self) -> List[int]:
        """Return the positions of the gaps in the sequence.
        
        Returns
        -------
        List[int]
            The positions of the gaps in the sequence.
            
        Examples
        --------
        >>> gapped_seq = NucleotideSequence.from_str('ATG-A-CCGTATGAA---TGA')
        >>> gapped_seq.gapped_positions()
        [3, 5, 15, 16, 17]
        """
        return super().gapped_positions()
    
    def count_gaps(self) -> int:
        """Return the number of gaps in the sequence.
        
        Returns
        -------
        int
            The number of gaps in the sequence.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.count_gaps()
        0
        >>> gapped_seq = NucleotideSequence.from_str('ATG-A-CCGTATGAA---TGA')
        >>> gapped_seq.count_gaps()
        5
        """
        return super().count_gaps()

    # NucleotideSequence-specific methods
    
    def to_reverse_complement(self, unit: Type = Nucleotide) -> 'NucleotideSequence':
        """Return the reverse complement of the sequence.
        
        Returns
        -------
        NucleotideSequence
            The reverse complement of the sequence.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.to_reverse_complement()
        TCATTCATACGGCAT
        """
        rc_list = [cast(unit, s).to_complement() for s in self._sequence[::-1]]
        return self.__class__(
            rc_list, 
            self._is_standard, 
            self._is_degenerate, 
            self._is_gapped, 
            self._is_masked)

    def to_codon_sequence(self, unit: Type = Codon) -> 'CodonSequence':
        """Return the sequence as a CodonSequence.
        
        Returns
        -------
        CodonSequence
            The sequence as a CodonSequence.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq.to_codon_sequence()
        ATG CCG TAT GAA TGA
        >>> invalid_seq = NucleotideSequence.from_str('ATGCCGTATGAATG')
        >>> invalid_seq.to_codon_sequence()
        Traceback (most recent call last):
        ValueError: Cannot convert to CodonSequence. Sequence length must be a multiple of 3: 14
        >>> gapped_seq = NucleotideSequence.from_str('ATG---TATGAATGA')
        >>> gapped_seq.to_codon_sequence()
        ATG --- TAT GAA TGA
        >>> invalid_gapped_seq = NucleotideSequence.from_str('AT---CTATGAATGA')
        >>> invalid_gapped_seq.to_codon_sequence()
        Traceback (most recent call last):
        ValueError: Invalid nucleotide sequence triplet: [A, T, -]
        """
        if len(self) % 3 != 0:
            raise ValueError(f'Cannot convert to CodonSequence. Sequence length must be a multiple of 3: {len(self)}')
        return CodonSequence.from_nucleotides(self, unit)


class AminoAcidSequence(BioSequence):
    """
    A sequence of amino acids.
    
    Parameters
    ----------
    sequence : Sequence[AminoAcid]
        The sequence of amino acids as AminoAcid objects.
    is_standard : bool, optional
        Whether the sequence is standard (i.e. only contains standard amino acids).
    is_ambiguous : bool, optional
        True if the sequence contains ambiguous nucleotides
    is_gapped : bool, optional
        True if the sequence is gapped (i.e. contains at least one gap).
    is_masked : bool, optional
        True if the sequence is masked (i.e. contains at least one the masked position).
    """    
    def __init__(self, 
            sequence: Sequence[AminoAcidEnum],
            is_standard: Optional[bool] = None,
            is_ambiguous: Optional[bool] = None, 
            is_gapped: Optional[bool] = None, 
            is_masked: Optional[bool] = None):
        super().__init__(sequence, is_standard, is_ambiguous, is_gapped, is_masked)
            
    @classmethod
    def from_str(cls, sequence: str, unit: Type = AminoAcid) -> 'AminoAcidSequence':
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
        return cast(AminoAcidSequence, super().from_str(sequence, unit))
    
    @classmethod
    def from_onehot(cls, sequence: Sequence[Sequence[int]], unit: Type = AminoAcid) -> 'AminoAcidSequence':
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
        return cast(AminoAcidSequence, super().from_onehot(sequence, unit))
    
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
        return super().to_str()
        
    def to_onehot(self) -> Sequence[Sequence[int]]:
        return super().to_onehot()

    def startswith(self, seq: Union[str, Sequence[AminoAcidEnum], 'AminoAcidSequence']) -> bool:
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
        >>> seq.startswith('A')
        True
        >>> seq.startswith([AminoAcid.Ala])
        True
        >>> seq.startswith('ARN')
        True
        >>> other_seq = AminoAcidSequence.from_str('ARN')
        >>> seq.startswith(other_seq)
        True
        >>> seq.startswith('ARNDCEQGHILKMFPSTWYVX')
        False
        >>> seq.startswith('-')
        False
        """
        if isinstance(seq, str):
            seq = cast(Sequence[AminoAcidEnum], self.from_str(seq, self.unit).sequence)
        elif isinstance(seq, self.__class__):
            seq = cast(Sequence[AminoAcidEnum], seq.sequence)
        return super().startswith(seq)

    def endswith(self, seq: Union[str, Sequence[AminoAcidEnum], 'AminoAcidSequence']) -> bool:
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
        >>> seq.endswith([AminoAcid.Val])
        True
        >>> seq.endswith('WYV')
        True
        >>> other_seq = AminoAcidSequence.from_str('WYV')
        >>> seq.endswith(other_seq)
        True
        >>> seq.endswith('ARNDC--EQGHILKMFPSTWYV')
        False
        >>> seq.endswith('-')
        False
        """
        if isinstance(seq, str):
            seq = cast(Sequence[AminoAcidEnum], self.from_str(seq, self.unit).sequence)
        elif isinstance(seq, self.__class__):
            seq = cast(Sequence[AminoAcidEnum], seq.sequence)
        return super().endswith(seq)

    def find(self, seq: Union[str, Sequence[AminoAcidEnum], 'AminoAcidSequence']) -> int:
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], AminoAcidEnum):
            seq = ''.join([str(s) for s in seq])
        return super().find(cast(str, seq))

    def rfind(self, seq: Union[str, Sequence[AminoAcidEnum], 'AminoAcidSequence']) -> int:
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], AminoAcidEnum):
            seq = ''.join([str(s) for s in seq])
        return super().rfind(cast(str, seq))
    
    def count(self, seq: Union[str, Sequence[AminoAcidEnum], 'AminoAcidSequence']) -> int:
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], AminoAcidEnum):
            seq = ''.join([str(s) for s in seq])
        return super().count(cast(str, seq))
    
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
        return cast(AminoAcidSequence, super().mask(positions))
            
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
        return super().masked_positions()
    
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
        return super().count_masked()
    
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
        return super().gapped_positions()
    
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
        return super().count_gaps()


class CodonSequence(BioSequence):
    """A sequence of codons.
    
    Parameters
    ----------
    sequence : Sequence[Codon]
        The sequence of codons as Codon objects.
    is_standard : bool, optional
        Whether the sequence is standard (i.e. only contains standard codons).
    is_ambiguous : bool, optional
        True if the sequence contains ambiguous codons.
    is_gapped : bool, optional
        True if the sequence is gapped (i.e. contains at least one gap).
    is_masked : bool, optional
        True if the sequence is masked (i.e. contains at least one the masked position).
    """
    def __init__(self, 
            sequence: Sequence[CodonEnum],
            is_standard: Optional[bool] = None,
            is_ambiguous: Optional[bool] = None, 
            is_gapped: Optional[bool] = None, 
            is_masked: Optional[bool] = None):
        super().__init__(sequence, is_standard, is_ambiguous, is_gapped, is_masked)

    def __repr__(self) -> str:
        return ' '.join([str(c) for c in self._sequence])
    
    @staticmethod
    def unit_iterator(sequence) -> Iterator:
        return (sequence[i:i+3] for i in range(0, len(sequence), 3))
    
    @classmethod
    def from_str(cls, sequence: str, unit: Type = Codon) -> 'CodonSequence':
        """Create a CodonSequence from a string of triplet nucleotides
        
        Parameters
        ----------
        sequence : str
            The string of triplet nucleotides.
            
        Returns
        -------
        CodonSequence
            The CodonSequence.
        
        Examples
        --------
        >>> CodonSequence.from_str('ATGAAATAG')
        ATG AAA TAG
        """
        if len(sequence) % 3 != 0:
            raise ValueError('Sequence length must be a multiple of 3')
        return cast(CodonSequence, 
                    super().from_str(cls.unit_iterator(sequence), unit))

    @classmethod
    def from_onehot(cls, sequence: Sequence[Sequence[int]], unit: Type = Codon) -> 'CodonSequence':
        """Create a CodonSequence from a one-hot encoded sequence.
        
        Parameters
        ----------
        sequence : Sequence[Sequence[int]]
            The one-hot encoded sequence.
            
        Returns
        -------
        CodonSequence
            The CodonSequence.
        """
        return cast(CodonSequence, 
                    super().from_onehot(sequence, unit))

    def startswith(self, seq: Union[str, Sequence[CodonEnum], 'CodonSequence']) -> bool:
        """Return True if the sequence starts with the given sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Codon], CodonSequence]
            The sequence to check.
            
        Returns
        -------
        bool
            True if the sequence starts with the given sequence.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAATAG')
        >>> seq.startswith('ATG')
        True
        >>> seq.startswith([Codon.ATG])
        True
        >>> seq.startswith(CodonSequence([Codon.ATG]))
        True
        >>> seq.startswith('ATGAAA')
        True
        >>> seq.startswith('ATGAAATAGTTT')
        False
        >>> seq.startswith('TGA')
        False
        """
        if isinstance(seq, str):
            seq = cast(Sequence[CodonEnum], self.from_str(seq, self.unit).sequence)
        elif isinstance(seq, self.__class__):
            seq = cast(Sequence[CodonEnum], seq.sequence)
        return super().startswith(seq)

    def endswith(self, seq: Union[str, Sequence[CodonEnum], 'CodonSequence']) -> bool:
        """Return True if the sequence ends with the given sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Codon], CodonSequence]
            The sequence to check.
            
        Returns
        -------
        bool
            True if the sequence ends with the given sequence.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAATAG')
        >>> seq.endswith('TAG')
        True
        >>> seq.endswith([Codon.TAG])
        True
        >>> seq.endswith(CodonSequence([Codon.TAG]))
        True
        >>> seq.endswith('AAATAG')
        True
        >>> seq.endswith('TTTATGAAATAG')
        False
        >>> seq.endswith('TGA')
        False
        """
        if isinstance(seq, str):
            seq = cast(Sequence[CodonEnum], self.from_str(seq, self.unit).sequence)
        elif isinstance(seq, self.__class__):
            seq = cast(Sequence[CodonEnum], seq.sequence)
        return super().endswith(seq)
    
    def find(self, seq: Union[str, Sequence[CodonEnum], 'CodonSequence']) -> int:
        """Return the lowest index in the sequence where the given subsequence is found.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Codon], CodonSequence]
            The sequence to find.
            
        Returns
        -------
        int
            The lowest index in the sequence where the given subsequence is found.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.find('ATG')
        0
        >>> seq.find(CodonSequence([Codon.ATG]))
        0
        >>> seq.find('AAA')
        1
        >>> seq.find('TAG')
        3
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, CodonSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], CodonEnum):
            seq = ''.join([str(s) for s in seq])
        return super().find(cast(str, seq)) // 3
        
    def rfind(self, seq: Union[str, Sequence[CodonEnum], 'CodonSequence']) -> int:
        """Return the highest index in the sequence where the given subsequence is found.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Codon], CodonSequence]
            The sequence to find.
            
        Returns
        -------
        int
            The highest index in the sequence where the given subsequence is found.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.rfind('ATG')
        0
        >>> seq.rfind(CodonSequence([Codon.ATG]))
        0
        >>> seq.rfind('AAA')
        2
        >>> seq.find('TAG')
        3
        """
        if isinstance(seq, str):
            pass
        elif isinstance(seq, CodonSequence):
            seq = str(seq)
        elif isinstance(seq, Sequence) and isinstance(seq[0], CodonEnum):
            seq = ''.join([str(s) for s in seq])
        return super().rfind(cast(str, seq)) // 3

    def count(self, seq: Union[str, Sequence[CodonEnum], 'CodonSequence']) -> int:
        """Return the number of non-overlapping occurrences of the given subsequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Codon], CodonSequence]
            The sequence to count.
            
        Returns
        -------
        int
            The number of non-overlapping occurrences of the given subsequence.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.count(CodonSequence([Codon.ATG]))
        1
        >>> seq.count('AAA')
        2
        >>> seq.count([Codon.TAG])
        1
        """
        if isinstance(seq, str):
            if len(seq) == 0:
                return len(self) + 1
            seq = CodonSequence.from_str(seq)
        elif isinstance(seq, CodonSequence):
            pass
        elif isinstance(seq, Sequence) and isinstance(seq[0], CodonEnum):
            seq = CodonSequence(seq)
        seq = cast(CodonSequence, seq)
        return sum(self.sequence[i:i+len(seq)] == seq.sequence
                   for i in range(0, len(self)-len(seq)+1))

    def mask(self, positions: Union[int, Sequence[int]]) -> 'CodonSequence':
        """Return a new sequence with the given positions masked.
        
        Parameters
        ----------
        positions : Union[int, Sequence[int]]
            The positions to mask.
            
        Returns
        -------
        CodonSequence
            A new sequence with the given positions masked.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.mask(0)
        ### AAA AAA TAG
        >>> seq.mask([0, 1, 2])
        ### ### ### TAG
        """
        return cast(CodonSequence, super().mask(positions))
    
    def masked_positions(self) -> List[int]:
        """Return a list of masked positions.
        
        Returns
        -------
        List[int]
            List of positions that are masked.
            
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.masked_positions()
        []
        >>> seq.mask(0).masked_positions()
        [0]
        >>> seq.mask([0, 1, 2]).masked_positions()
        [0, 1, 2]
        """
        return super().masked_positions()

    def count_masked(self) -> int:
        """Return the number of masked positions.
        
        Returns
        -------
        int
            The number of masked positions.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.count_masked()
        0
        >>> seq.mask(0).count_masked()
        1
        >>> seq.mask([0, 1, 2]).count_masked()
        3
        """
        return super().count_masked()

    def gapped_positions(self) -> List[int]:
        """Return a list of gapped positions.
        
        Returns
        -------
        List[int]
            List of positions that are gapped.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.gapped_positions()
        []
        >>> gapped_seq = CodonSequence.from_str('ATG---AAA---TAG')
        >>> gapped_seq.gapped_positions()
        [1, 3]
        """
        return super().gapped_positions()

    def count_gaps(self) -> int:
        """Return the number of gapped positions.
        
        Returns
        -------
        int
            The number of gapped positions.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.count_gaps()
        0
        >>> gapped_seq = CodonSequence.from_str('ATG---AAA---TAG')
        >>> gapped_seq.count_gaps()
        2
        """
        return super().count_gaps()

    # CodonSequence-specific methods

    @classmethod
    def from_nucleotides(cls, 
            sequence: Sequence[NucleotideEnum], 
            unit: Type = Codon) -> 'CodonSequence':
        """Create a CodonSequence from a NucleotideSequence.
        
        Parameters
        ----------
        sequence : Sequence[NucleotideEnum]
            Sequence of Nucleotide objects.
            
        Returns
        -------
        CodonSequence
            The CodonSequence.
            
        Examples
        --------
        >>> seq = [Nucleotide.A, Nucleotide.T, Nucleotide.G, Nucleotide.A, Nucleotide.A, Nucleotide.A, Nucleotide.T, Nucleotide.A, Nucleotide.G]
        >>> CodonSequence.from_nucleotides(seq)
        ATG AAA TAG
        """
        if len(sequence) % 3 != 0:
            raise ValueError('Sequence length must be a multiple of 3')
        return cls([unit.from_nucleotides(s) for s in cls.unit_iterator(sequence)])

    @classmethod
    def from_nucleotide_onehot(cls, 
            sequence: Sequence[Sequence[int]], 
            unit: Type = Codon) -> 'CodonSequence':
        """Create a CodonSequence from seqeunce of Nucleotide one-hot vectors.
        
        Parameters
        ----------
        sequence : Sequence[Sequence[int]]
            The one-hot encoded sequence.
            
        Returns
        -------
        CodonSequence
            The CodonSequence.
        """
        if len(sequence) % 3 != 0:
            raise ValueError('Sequence length must be a multiple of 3')
        nucl_sequence = [NucleotideEnum.from_onehot(s) for s in sequence]
        return cls([unit.from_nucleotides(s) for s in cls.unit_iterator(nucl_sequence)])

    def translate(self) -> AminoAcidSequence:
        """Return the amino acid sequence translated from the codon sequence.
        
        Returns
        -------
        AminoAcidSequence
            The amino acid sequence translated from the codon sequence.
            
        Examples
        --------
        >>> seq = CodonSequence.from_str('ATGAAAAAATAG')
        >>> seq.translate()
        MKK*
        >>> gapped_seq = CodonSequence.from_str('ATGAAA---AAATAG')
        >>> gapped_seq.translate()
        MK-K*
        >>> masked_seq = CodonSequence.from_str('ATGAAA###AAATAG')
        >>> masked_seq.translate()
        MK#K*
        """
        return AminoAcidSequence([c.translate() for c in self._sequence])  # type: ignore

    def to_reverse_complement(self) -> 'CodonSequence':
        """Return the reverse complement of the sequence.
        
        Returns
        -------
        NucleotideSequence
            The reverse complement of the sequence.
            
        Examples
        --------
        >>> str_seq = 'ATGCCGTATGAATGA'
        >>> seq = CodonSequence.from_str(str_seq)
        >>> rc = seq.to_reverse_complement()
        >>> rc
        TCA TTC ATA CGG CAT
        >>> rc_rc = rc.to_reverse_complement()
        >>> rc_rc
        ATG CCG TAT GAA TGA
        >>> str_seq == str(rc_rc)
        True
        """
        rc_nucl_str_list = [COMPLEMENTARY_NUCLEOTIDES[s] for s in str(self)[::-1]]
        rc_list = [
            self.unit.from_str(''.join(rc_nucl_str_list[i:i+3])) 
            for i in range(0, len(rc_nucl_str_list), 3)]
        return CodonSequence(
            rc_list,  # type: ignore
            self._is_standard, 
            self._is_degenerate, 
            self._is_gapped, 
            self._is_masked)
    
    def to_amino_acid_sequence(self) -> 'AminoAcidSequence':
        """Translates the codon sequence to an amino acid sequence.
        
        Returns
        -------
        AminoAcidSequence
            The translated amino acid sequence.
        
        See also
        --------
        translate
        """
        return self.translate()
        
    def to_nucleotide_sequence(self) -> 'NucleotideSequence':
        """Return the sequence as a NuclotideSequence.
        
        Returns
        -------
        NuclotideSequence
            The sequence as a NuclotideSequence.
            
        Examples
        --------
        >>> str_seq = 'ATGCCG---TATGAATGA'
        >>> seq = CodonSequence.from_str(str_seq)
        >>> nucl_seq = seq.to_nucleotide_sequence()
        >>> nucl_seq
        ATGCCG---TATGAATGA
        >>> str_seq == str(nucl_seq)
        True
        """
        return NucleotideSequence.from_str(str(self))
