from typing import Union, List
from collections.abc import Sequence
from enum import Enum

from biodatatypes.sequence import SeqMixin, MaskMixin, GapMixin
from biodatatypes.constants.aminoacid import *
from biodatatype.units import AminoAcid


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
