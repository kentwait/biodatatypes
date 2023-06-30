from typing import List, Union, Any
from enum import Enum
from collections.abc import Sequence

from biodatatypes import Nucleotide, AminoAcid, AminoAcidSequence

from biodatatypes.sequence import SeqMixin, MaskMixin, GapMixin
from biodatatypes.constants.codon import *
from biodatatype.units import Codon


class CodonSequence(GapMixin, MaskMixin, SeqMixin, Sequence):
    unit = Codon

    def __init__(self, 
            sequence: Sequence[Codon],
            is_standard: bool = None, 
            is_degenerate: bool = None, 
            is_gapped: bool = None, 
            is_masked: bool = None):
        self._sequence = list(sequence)
        self._is_standard = is_standard
        self._is_degenerate = is_degenerate
        self._is_gapped = is_gapped
        self._is_masked = is_masked

    def __repr__(self) -> str:
        return ' '.join([str(c) for c in self._sequence])
    
    @classmethod
    def from_str(cls, sequence: str) -> 'CodonSequence':
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
        return cls([Codon.from_str(sequence[i:i+3]) for i in range(0, len(sequence), 3)])
    
    @classmethod
    def from_nucleotides(cls, sequence: Sequence[Nucleotide]) -> 'CodonSequence':
        """Create a CodonSequence from a NucleotideSequence.
        
        Parameters
        ----------
        sequence : Sequence[Nucleotide]
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
        return cls([Codon.from_nucleotides(sequence[i:i+3]) for i in range(0, len(sequence), 3)])    

    @classmethod
    def from_onehot(cls, sequence: Sequence[Sequence[int]]) -> 'CodonSequence':
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
        return cls([Codon.from_onehot(s) for s in sequence])

    @classmethod
    def from_nucleotide_onehot(cls, sequence: Sequence[Sequence[int]]) -> 'CodonSequence':
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
        sequence = [Nucleotide.from_onehot(s) for s in sequence]
        return cls([Codon.from_nucleotides(sequence[i:i+3]) for i in range(0, len(sequence), 3)])

    def startswith(self, seq: Union[str, Sequence[Codon], 'CodonSequence']) -> bool:
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
        >>> seq.startswith(CodonSequence([Codon.ATG]))
        True
        >>> seq.startswith('ATGAAA')
        True
        """
        if isinstance(seq, str):
            seq = CodonSequence.from_str(seq)
        return SeqMixin.startswith(self, seq)

    def endswith(self, seq: Union[str, Sequence[Codon], 'CodonSequence']) -> bool:
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
        >>> seq.endswith(CodonSequence([Codon.TAG]))
        True
        >>> seq.endswith('AAATAG')
        True
        """
        if isinstance(seq, str):
            seq = CodonSequence.from_str(seq)
        return SeqMixin.endswith(self, seq)
    
    def find(self, seq: Union[str, Sequence[Codon], 'CodonSequence']) -> int:
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], Codon):
            seq = ''.join([str(s) for s in seq])
        return SeqMixin.find(self, seq) // 3
        
    def rfind(self, seq: Union[str, Sequence[Codon], 'CodonSequence']) -> int:
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], Codon):
            seq = ''.join([str(s) for s in seq])
        return SeqMixin.rfind(self, seq) // 3

    def count(self, seq: Union[str, Sequence[Codon], 'CodonSequence']) -> int:
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
        elif isinstance(seq, Sequence) and isinstance(seq[0], Codon):
            seq = CodonSequence(seq)
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
        return MaskMixin.mask(self, positions)
    
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
        return MaskMixin.masked_positions(self)

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
        return MaskMixin.count_masked(self)

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
        return GapMixin.gapped_positions(self)

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
        return GapMixin.count_gaps(self)

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
        return AminoAcidSequence([c.translate() for c in self._sequence])
