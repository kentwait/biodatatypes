from typing import Union, List, Any, Iterator, Iterable
from collections.abc import Sequence

from biodatatypes.units import Nucleotide, AminoAcid, Codon


class BioSequence(Sequence):
    def __init__(self, 
            sequence: Sequence, 
            is_standard: bool = None,
            is_degenerate: bool = None, 
            is_gapped: bool = None, 
            is_masked: bool = None):
        self._sequence = list(sequence)
        self._is_standard = is_standard
        self._is_degenerate = is_degenerate
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

    @property
    def is_standard(self) -> bool:
        if self._is_standard is None:
            self._is_standard = all([s.is_standard() for s in self.sequence])
        return self._is_standard
    
    @property
    def is_degenerate(self) -> bool:
        if self._is_degenerate is None:
            self._is_degenerate = any([s.is_degenerate() for s in self.sequence])
        return self._is_degenerate
    
    @property
    def is_gapped(self) -> bool:
        if self._is_gapped is None:
            self._is_gapped = any([s.is_gap() for s in self.sequence])
        return self._is_gapped
    
    @property
    def is_masked(self) -> bool:
        if self._is_masked is None:
            self._is_masked = any([s.is_masked() for s in self.sequence])
        return self._is_masked
    
    @property
    def sequence(self) -> Sequence:
        return self._sequence
    
    @classmethod
    def from_str(cls, sequence: Iterable[str]) -> 'BioSequence':
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
        return cls([cls.unit.from_str(s) for s in sequence])
    
    @classmethod
    def from_onehot(cls, sequence: Sequence[Sequence[int]]) -> 'BioSequence':
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
        return cls([cls.unit.from_onehot(s) for s in sequence])
    
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
        return self.sequence[:len(seq)] == seq.sequence

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
        return self.sequence[-len(seq):] == seq.sequence

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
    unit = Nucleotide
    
    def __init__(self, 
            sequence: Sequence[Nucleotide], 
            is_standard: bool = None, 
            is_degenerate: bool = None, 
            is_gapped: bool = None, 
            is_masked: bool = None):
        super().__init__(sequence, is_standard, is_degenerate, is_gapped, is_masked)
        
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
        return super().from_str(sequence)
    
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
        return super().from_onehot(sequence)
    
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
        return CodonSequence.from_nucleotides(self)
    
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
        return super().startswith(seq)

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
        return super().endswith(seq)

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
        return super().find(seq)

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
        return super().rfind(seq)

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
        return super().count(seq)

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
        return super().mask(positions)
    
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


class AminoAcidSequence(BioSequence):
    unit = AminoAcid
    
    def __init__(self, 
            sequence: Sequence[AminoAcid],
            is_standard: bool = None, 
            is_degenerate: bool = None, 
            is_gapped: bool = None, 
            is_masked: bool = None):
        super().__init__(sequence, is_standard, is_degenerate, is_gapped, is_masked)
            
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
        return super().from_str(sequence)
    
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
        return super().from_onehot(sequence)
    
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
        return super().startswith(seq)

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
        return super().endswith(seq)

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
        return super().find(seq)

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
        return super().rfind(seq)
    
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
        return super().count(seq)
    
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
        return super().mask(positions)
            
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
    unit = Codon

    def __init__(self, 
            sequence: Sequence[Codon],
            is_standard: bool = None, 
            is_degenerate: bool = None, 
            is_gapped: bool = None, 
            is_masked: bool = None):
        super().__init__(sequence, is_standard, is_degenerate, is_gapped, is_masked)

    def __repr__(self) -> str:
        return ' '.join([str(c) for c in self._sequence])
    
    @staticmethod
    def unit_iterator(sequence: Sequence[Codon]) -> Iterator[Codon]:
        return (sequence[i:i+3] for i in range(0, len(sequence), 3))
    
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
        return super().from_str(cls.unit_iterator(sequence))
    
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
        return cls([Codon.from_nucleotides(s) for s in cls.unit_iterator(sequence)])

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
        return super().from_onehot(cls.unit_iterator(sequence))

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
        return cls([Codon.from_nucleotides(s) for s in cls.unit_iterator(sequence)])

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
        return super().startswith(seq)

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
        return super().endswith(seq)
    
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
        return super().find(seq) // 3
        
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
        return super().rfind(seq) // 3

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
        return super().mask(positions)
    
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

