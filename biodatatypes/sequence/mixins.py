from typing import List, Sequence, Union, Any


class SeqMixin:
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
    

class MaskMixin:
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


class GapMixin:
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
