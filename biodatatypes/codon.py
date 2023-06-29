from enum import Enum
from collections.abc import Sequence

from biodatatypes import Nucleotide, AminoAcid


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
    
    def __eq__(self, other):
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


# Constants

TRANSLATION_TABLE = {
    'AAA': 'Lys',
    'AAC': 'Asn',
    'AAG': 'Lys',
    'AAT': 'Asn',
    'ACA': 'Thr',
    'ACC': 'Thr',
    'ACG': 'Thr',
    'ACT': 'Thr',
    'AGA': 'Arg',
    'AGC': 'Ser',
    'AGG': 'Arg',
    'AGT': 'Ser',
    'ATA': 'Ile',
    'ATC': 'Ile',
    'ATG': 'Met',
    'ATT': 'Ile',
    'CAA': 'Gln',
    'CAC': 'His',
    'CAG': 'Gln',
    'CAT': 'His',
    'CCA': 'Pro',
    'CCC': 'Pro',
    'CCG': 'Pro',
    'CCT': 'Pro',
    'CGA': 'Arg',
    'CGC': 'Arg',
    'CGG': 'Arg',
    'CGT': 'Arg',
    'CTA': 'Leu',
    'CTC': 'Leu',
    'CTG': 'Leu',
    'CTT': 'Leu',
    'GAA': 'Glu',
    'GAC': 'Asp',
    'GAG': 'Glu',
    'GAT': 'Asp',
    'GCA': 'Ala',
    'GCC': 'Ala',
    'GCG': 'Ala',
    'GCT': 'Ala',
    'GGA': 'Gly',
    'GGC': 'Gly',
    'GGG': 'Gly',
    'GGT': 'Gly',
    'GTA': 'Val',
    'GTC': 'Val',
    'GTG': 'Val',
    'GTT': 'Val',
    'TAA': 'Ter',
    'TAC': 'Tyr',
    'TAG': 'Ter',
    'TAT': 'Tyr',
    'TCA': 'Ser',
    'TCC': 'Ser',
    'TCG': 'Ser',
    'TCT': 'Ser',
    'TGA': 'Ter',
    'TGC': 'Cys',
    'TGG': 'Trp',
    'TGT': 'Cys',
    'TTA': 'Leu',
    'TTC': 'Phe',
    'TTG': 'Leu',
    'TTT': 'Phe',
    # Special triplet tokens
    'Mask': 'Mask',
    'Gap': 'Gap',
    'Special': 'Special',
    # Ambiguous triplet tokens
    'NNN': 'Xaa'   
}
FROM_TRIPLET_TOKEN = {
    'AAA': 'AAA',
    'AAC': 'AAC',
    'AAG': 'AAG',
    'AAT': 'AAT',
    'ACA': 'ACA',
    'ACC': 'ACC',
    'ACG': 'ACG',
    'ACT': 'ACT',
    'AGA': 'AGA',
    'AGC': 'AGC',
    'AGG': 'AGG',
    'AGT': 'AGT',
    'ATA': 'ATA',
    'ATC': 'ATC',
    'ATG': 'ATG',
    'ATT': 'ATT',
    'CAA': 'CAA',
    'CAC': 'CAC',
    'CAG': 'CAG',
    'CAT': 'CAT',
    'CCA': 'CCA',
    'CCC': 'CCC',
    'CCG': 'CCG',
    'CCT': 'CCT',
    'CGA': 'CGA',
    'CGC': 'CGC',
    'CGG': 'CGG',
    'CGT': 'CGT',
    'CTA': 'CTA',
    'CTC': 'CTC',
    'CTG': 'CTG',
    'CTT': 'CTT',
    'GAA': 'GAA',
    'GAC': 'GAC',
    'GAG': 'GAG',
    'GAT': 'GAT',
    'GCA': 'GCA',
    'GCC': 'GCC',
    'GCG': 'GCG',
    'GCT': 'GCT',
    'GGA': 'GGA',
    'GGC': 'GGC',
    'GGG': 'GGG',
    'GGT': 'GGT',
    'GTA': 'GTA',
    'GTC': 'GTC',
    'GTG': 'GTG',
    'GTT': 'GTT',
    'TAA': 'TAA',
    'TAC': 'TAC',
    'TAG': 'TAG',
    'TAT': 'TAT',
    'TCA': 'TCA',
    'TCC': 'TCC',
    'TCG': 'TCG',
    'TCT': 'TCT',
    'TGA': 'TGA',
    'TGC': 'TGC',
    'TGG': 'TGG',
    'TGT': 'TGT',
    'TTA': 'TTA',
    'TTC': 'TTC',
    'TTG': 'TTG',
    'TTT': 'TTT',
    '---': 'Gap',
    '###': 'Mask',
    '@@@': 'Special',
    'NNN': 'NNN',
}
STOP_CODONS = ['TAA', 'TAG', 'TGA']

# Degenerate codons
# Unique codons
ONEFOLD_DEGENERATE_CODONS = [
    'ATG',  # Met
    'TGG',  # Trp
]
# 2 codons code for the same amino acid
TWOFOLD_DEGENERATE_CODONS = [
    'AAT', 'AAC',  # Asn
    'GAT', 'GAC',  # Asp
    'TGT', 'TGC',  # Cys
    'CAA', 'CAG',  # Gln
    'GAA', 'GAG',  # Glu
    'CAT', 'CAC',  # His
    'AAA', 'AAG',  # Lys
    'TTT', 'TTC',  # Phe
    'TAT', 'TAC',  # Tyr
]
# 3 codons code for the same amino acid
THREEFOLD_DEGENERATE_CODONS = [
    'ATT', 'ATC', 'ATA',  # Ile
]
# 4 codons code for the same amino acid
FOURFOLD_DEGENERATE_CODONS = [
    'GCT', 'GCC', 'GCA', 'GCG',  # Ala
    'GGT', 'GGC', 'GGA', 'GGG',  # Gly
    'CCT', 'CCC', 'CCA', 'CCG',  # Pro
    'ACT', 'ACC', 'ACA', 'ACG',  # Thr
    'GTT', 'GTC', 'GTA', 'GTG',  # Val
]
# 6 codons code for the same amino acid
SIXFOLD_DEGENERATE_CODONS = [
    'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG',  # Leu
    'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',  # Arg
    'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',  # Ser
]