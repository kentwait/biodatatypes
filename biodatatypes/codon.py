from typing import List, Union, Any
from enum import Enum
from collections.abc import Sequence

from biodatatypes import Nucleotide, AminoAcid, AminoAcidSequence

from biodatatypes.sequence import SeqMixin, MaskMixin, GapMixin
from biodatatypes.constants.codon import *


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
