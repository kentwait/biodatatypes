from typing import Union, List
from collections.abc import Sequence
from enum import Enum

from biodatatypes.codon import CodonSequence
from biodatatypes.sequence import SeqMixin, MaskMixin, GapMixin
from biodatatypes.constants.nucleotide import *
from biodatatype.units import Nucleotide


class NucleotideSequence(GapMixin, MaskMixin, SeqMixin, Sequence):
    unit = Nucleotide
    
    def __init__(self, 
            sequence: Sequence[Nucleotide], 
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
    def from_str(cls, sequence: str) -> 'NucleotideSequence':
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
        return cls([Nucleotide.from_str(s) for s in sequence])
    
    @classmethod
    def from_onehot(cls, sequence: Sequence[Sequence[int]]) -> 'NucleotideSequence':
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
        return cls([Nucleotide.from_onehot(s) for s in sequence])
    
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
        return SeqMixin.to_str(self)
        
    def to_onehot(self) -> Sequence[Sequence[int]]:
        return SeqMixin.to_onehot(self)
    
    def to_codon_sequence(self) -> 'CodonSequence':
        """Converts the nucleotide sequence to a codon sequence.
        
        Returns
        -------
        CodonSequence
            A CodonSequence object.
            
        Examples
        --------
        >>> seq = NucleotideSequence.from_str('ATGCCGTATGAATGA')
        >>> seq
        ATGCCGTATGAATGA
        >>> seq.to_codon_sequence()
        ATG CCG TAT GAA TGA
        """
        if len(self) % 3 != 0:
            raise ValueError('Sequence length must be a multiple of 3.')
        return CodonSequence(self)
    
    def startswith(self, seq: Union[str, Sequence[Nucleotide], 'NucleotideSequence']) -> bool:
        """Return True if the sequence starts with the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Nucleotide]]
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
        >>> seq.startswith('ATG')
        True
        >>> seq.startswith('-')
        False
        >>> seq.startswith('A--')
        False
        """
        if isinstance(seq, str):
            seq = NucleotideSequence.from_str(seq)
        return SeqMixin.startswith(self, seq)

    def endswith(self, seq: Union[str, Sequence[Nucleotide], 'NucleotideSequence']) -> bool:
        """Return True if the sequence ends with the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Nucleotide]]
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
        >>> seq.endswith('TGA')
        True
        >>> seq.endswith('-')
        False
        >>> seq.endswith('A--')
        False
        """
        if isinstance(seq, str):
            seq = NucleotideSequence.from_str(seq)
        return SeqMixin.endswith(self, seq)

    def find(self, seq: Union[str, Sequence[Nucleotide], 'NucleotideSequence']) -> int:
        """Return the index of the first occurrence of the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Nucleotide]]
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], Nucleotide):
            seq = ''.join([str(n) for n in seq])
        return SeqMixin.find(self, seq)

    def rfind(self, seq: Union[str, Sequence[Nucleotide], 'NucleotideSequence']) -> int:
        """Return the index of the last occurrence of the given nucleotide or sequence.
        
        Parameters
        ----------
        seq : Union[str, Sequence[Nucleotide]]
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], Nucleotide):
            seq = ''.join([str(n) for n in seq])
        return SeqMixin.rfind(self, seq)

    def count(self, seq: Union[str, Sequence[Nucleotide], 'NucleotideSequence']) -> int:
        """Return the number of occurrences of the given nucleotide or sequence.
        
        Parameters
        ----------
        nucl : Union[str, Sequence[Nucleotide]]
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], Nucleotide):
            seq = ''.join([str(n) for n in seq])
        return SeqMixin.count(self, seq)

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
        return MaskMixin.mask(self, positions)
    
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
        return MaskMixin.masked_positions(self)
        
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
        return MaskMixin.count_masked(self)

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
        return GapMixin.gapped_positions(self)
    
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
        return GapMixin.count_gaps(self)
