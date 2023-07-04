from typing import cast, Type
from collections.abc import Sequence
from enum import Enum

from biodatatypes.constants.nucleotide import *
from biodatatypes.constants.aminoacid import *
from biodatatypes.constants.codon import *


class BioToken(Enum):
    """Base class for all biological sequence unit tokens as Enum members.

    This subclasses the Enum datatype but does not define any members.
    This class is intended to be subclassed for each type of unit token: 
    Nucleotide, AminoAcidEnum, Codon. 
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

    #region Methods that must be implemented by subclasses
    #region Magic methods
    def __str__(self) -> str:
        raise NotImplementedError('BioToken must be subclassed to implement __str__()')
    
    #endregion
    
    #region Factory methods
    @classmethod
    def from_str(cls, s: str) -> 'BioToken':
        raise NotImplementedError('BioToken must be subclassed to implement from_str()')
    
    #endregion
     
    #region Token type methods
    def is_standard(self) -> bool:
        raise NotImplementedError('BioToken must be subclassed to implement is_standard()')
    
    def is_any(self) -> bool:
        raise NotImplementedError('BioToken must be subclassed to implement is_any()')
        
    def is_ambiguous(self) -> bool:
        raise NotImplementedError('BioToken must be subclassed to implement is_ambiguous()')

    #endregion
    #endregion
    
    #region Methods with default implementations
    #region Magic methods
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
    
    #endregion
    
    #region Factory methods
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
    
    #endregion
    
    #region Conversion methods
    def to_onehot(self) -> Sequence[int]:
        onehot = [0] * len(self.__class__)
        onehot[self.value] = 1
        return onehot
        
    def to_int(self) -> int:
        return self.value
    
    #endregion
    #endregion


class NucleotideEnum(BioToken):
    #region Methods that require implementation
    def __str__(self) -> str:
        return TO_NUCLEOTIDE_ONE_LETTER_TOKEN[self.name]

    @classmethod
    def from_str(cls, s: str) -> 'NucleotideEnum':
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
        """
        try:
            return cls[FROM_NUCLEOTIDE_ONE_LETTER_TOKEN[s.upper()]]
        except KeyError:
            raise ValueError(f'Invalid nucleotide: {s}')

    def is_standard(self):
        """Return True if the nucleotide is a standard nucleotide.

        Returns
        -------
        bool
            True if the nucleotide is a standard nucleotide, False otherwise.
            
        See Also
        --------
        is_special : Checks if it is a special token
        """
        return self.name in STANDARD_NUCLEOTIDES

    def is_any(self):
        """Return True if it represents any nucleotide.
        
        Expects 'N" to be the token name representing any nucleotide.
        
        Returns
        -------
        bool
            True if the nucleotide token represents any nucleotide.
        
        See Also
        --------
        is_ambiguous : Checks if it is an ambiguous token
        """
        return self.name == 'N'
   
    def is_ambiguous(self):
        """Return True if it represents an ambiguous nucleotide.
        
        Subclasses should override this method to check if the token represents
        an ambiguous nucleotide.
        
        Returns
        -------
        bool
            True if the nucleotide token represents an ambiguous nucleotide.
        
        See Also
        --------
        is_any : Checks if the token represents any nucleotide
        """
        return self.name in AMBIGUOUS_NUCLEOTIDES
 
    #endregion

    #region Override methods with default implementations
    
    # Factory methods
    @classmethod
    def from_int(cls, i: int) -> 'NucleotideEnum':
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
        """
        return cast('NucleotideEnum', super().from_int(i))

    @classmethod
    def from_onehot(cls, onehot:Sequence[int]) -> 'NucleotideEnum':
        """Convert a onehot encoding to a Nucleotide.
        
        Parameters
        ----------
        onehot : Sequence[int]
            A onehot encoding of a nucleotide.
            
        Returns
        -------
        Nucleotide
            A Nucleotide object.
        """
        return cast('NucleotideEnum', super().from_onehot(onehot))
    
    # Conversion methods
    def to_onehot(self) -> Sequence[int]:
        """Convert a Nucleotide to a onehot encoding.
        
        Returns
        -------
        Sequence[int]
            A onehot encoding of the nucleotide.
        """
        return super().to_onehot()
    
    def to_int(self) -> int:
        """Convert a Nucleotide to an integer.
        
        Returns
        -------
        int
            An integer representing the nucleotide.
        """
        return super().to_int()
    
    #endregion
    
    #region NucleotideEnum-only methods

    # Properties
    
    def is_purine(self) -> bool:
        """Return True if the nucleotide is a purine.
        
        Returns
        -------
        bool
            True if the nucleotide is a purine, False otherwise.
        """
        return self.name == 'A' or self.name == 'G'

    def is_pyrimidine(self) -> bool:
        """Return True if the nucleotide is a pyrimidine.
        
        Returns
        -------
        bool
            True if the nucleotide is a pyrimidine, False otherwise.
        """
        return self.name == 'C' or self.name == 'T'

    def is_strong(self) -> bool:
        """Return True if the nucleotide is a nucleotide involved in strong triple H-bond interaction.
        
        Returns
        -------
        bool
            True if the nucleotide is C or G, False otherwise.
        """
        return self.name == 'C' or self.name == 'G'

    def is_weak(self) -> bool:
        """Return True if the nucleotide is a nucleotide involved in weak double H-bond interaction.
        
        Returns
        -------
        bool
            True if the nucleotide is A or T, False otherwise.
        """
        return self.name == 'A' or self.name == 'T'

    def is_amino(self) -> bool:
        return self.name == 'A' or self.name == 'C'

    def is_keto(self) -> bool:
        return self.name == 'G' or self.name == 'T'
    
    # Conversion
    
    def to_complement(self) -> 'NucleotideEnum':
        """Convert a Nucleotide to its complement.
        
        Returns
        -------
        Nucleotide
            The complement of the nucleotide.
        """
        complement_str = COMPLEMENTARY_NUCLEOTIDES[self.name]
        return self.from_str(complement_str)
    
    #endregion


class AminoAcidEnum(BioToken):
    #region Methods that must be implemented by subclasses
    def __str__(self) -> str:
        return NAME_TO_AMINO_ACID_ONE_LETTER_TOKEN[self.value]
    
    @classmethod
    def from_str(cls, s: str) -> 'AminoAcidEnum':
        """Convert a string to an AminoAcidEnum object.
        
        Parameters
        ----------
        s : str
            A string representing an amino acid. It can be a one-letter code or
            a three-letter code.

        Returns
        -------
        AminoAcidEnum
            An AminoAcidEnum object.
            
        Raises
        ------
        ValueError
            If the string is not a valid amino acid.
        """
        if len(s) == 1:
            return cls.from_one_letter(s)
        elif len(s) == 3:
            return cls.from_three_letter(s)
        else:
            raise ValueError(f'Invalid amino acid string: {s}')
    
    def is_standard(self) -> bool:
        """Return True if it is a standard amino acid.
            
        Returns
        -------
        bool
            True if the amino acid token represents a standard amino acid.
        """
        return self.name in STANDARD_AA
    
    def is_any(self):
        """Return True if it represents any amino acid.
        
        Expects 'Xaa' to be the enum name representing any amino acid.
        
        Returns
        -------
        bool
            True if the amino acid token represents any amino acid.
        """
        return self.name == 'Xaa'
    
    def is_ambiguous(self) -> bool:
        """Return True if it represents an ambiguous amino acid token.
        
        Subclasses should override this method if they have ambiguous amino acid tokens.
        
        Returns
        -------
        bool
            True if the amino acid token represents an ambiguous amino acid token.
        """
        return self.is_any()

    #endregion
    
    #region Override methods with default implementations 
    
    # Factory methods
    @classmethod
    def from_int(cls, i: int) -> 'AminoAcidEnum':
        """Convert an integer to an AminoAcidEnum object.
        
        Parameters
        ----------
        i : int
            An integer representing an amino acid.
            
        Returns
        -------
        AminoAcidEnum
            An AminoAcidEnum object.
            
        Raises
        ------
        ValueError
            If the integer is not a valid amino acid.
        """
        return cast('AminoAcidEnum', super().from_int(i))
    
    @classmethod
    def from_onehot(cls, onehot: Sequence[int]) -> 'AminoAcidEnum':
        """Convert a onehot vector to an AminoAcidEnum object.
        
        Parameters
        ----------
        onehot : Sequence[int]
            A onehot vector representing an amino acid.
            
        Returns
        -------
        AminoAcidEnum
            An AminoAcidEnum object.
            
        Raises
        ------
        ValueError
            If the onehot vector is not a valid amino acid.
        """
        return cast('AminoAcidEnum', super().from_onehot(onehot))

    # Conversion methods
    def to_onehot(self) -> Sequence[int]:
        """Convert an AminoAcidEnum object to a onehot vector.
        
        Returns
        -------
        Sequence[int]
            A onehot vector representing the amino acid.
        """
        return super().to_onehot()
    
    def to_int(self) -> int:
        """Convert an AminoAcidEnum object to an integer representation.
        
        Returns
        -------
        int
            An integer representing the amino acid.
        """
        return super().to_int()
    
    #endregion
    
    #region Methods unique to AminoAcidEnum
    
    # Amino acid properties
    def is_stop(self) -> bool:
        """Return True if it represents a terminator sequence.
        
        Returns
        -------
        bool
            True if the amino acid token is a terminator sequence.
        """
        return self.name == 'Stop'
    
    def is_polar(self) -> bool:
        """Return True if the amino acid identity is known and is polar.
        
        Returns
        -------
        bool
            True if the amino acid token represents a polar amino acid.
        """
        return self.name in POLAR_AA

    def is_nonpolar(self) -> bool:
        """Return True if the amino acid identity is known and is nonpolar.
        
        Returns
        -------
        bool
            True if the amino acid token represents a nonpolar amino acid.
        """
        return self.name in NONPOLAR_AA
    
    def is_acidic(self) -> bool:
        """Return True if the amino acid identity is known and is acidic.
        
        Returns
        -------
        bool
            True if the amino acid token represents an acidic amino acid.
        """
        return self.name in ACIDIC_AA
    
    def is_basic(self) -> bool:
        """Return True if the amino acid identity is known and is basic.
        
        Returns
        -------
        bool
            True if the amino acid token represents a basic amino acid.
        """
        return self.name in BASIC_AA
    
    def is_aromatic(self) -> bool:
        """Return True if the amino acid identity is known and is aromatic.
        
        Returns
        -------
        bool
            True if the amino acid token represents an aromatic amino acid.
        """
        return self.name in AROMATIC_AA
    
    def is_hydrophobic(self) -> bool:
        """Return True if the amino acid identity is known and is hydrophobic.
        
        Returns
        -------
        bool
            True if the amino acid token represents a hydrophobic amino acid.
        """
        return self.name in HYDROPHOBIC_AA
    
    def has_sulfur(self) -> bool:
        """Return True if the amino acid identity is known and contains sulfur.
        
        Returns
        -------
        bool
            True if the amino acid token represents an amino acid with sulfur.
        """
        return self.name in SULFUR_AA
    
    def has_amide(self) -> bool:
        """Return True if the amino acid identity is known and contains an amide.
        
        Returns
        -------
        bool
            True if the amino acid token represents an amino acid with an amide.
        """
        return self.name in AMIDE_AA
    
    def has_hydroxyl(self) -> bool:
        """Return True if the amino acid identity is known and contains a hydroxyl.
        
        Returns
        -------
        bool
            True if the amino acid token represents an amino acid with a hydroxyl.
        """
        return self.name in HYDROXYL_AA
        
    # Factory methods unique to AminoAcidEnum
    @classmethod
    def from_one_letter(cls, one_letter_code: str) -> 'AminoAcidEnum':
        """Convert a one-letter code to an AminoAcidEnum object.
        
        Parameters
        ----------
        one_letter_code : str
            A one-letter code representing an amino acid.
        
        Returns
        -------
        AminoAcidEnum
            An AminoAcidEnum object.
            
        Raises
        ------
        ValueError
            If the one-letter code is not a valid amino acid.
        """
        try:
            return cls[AMINO_ACID_ONE_LETTER_TOKEN_TO_NAME[one_letter_code.upper()]]
        except KeyError:
            raise ValueError(f'Invalid one-letter code: {one_letter_code}')
    
    @classmethod
    def from_three_letter(cls, three_letter_code: str) -> 'AminoAcidEnum':
        """Convert a three-letter code to an AminoAcidEnum object.
        
        Parameters
        ----------
        three_letter_code : str
            A three-letter code representing an amino acid.
            
        Returns
        -------
        AminoAcidEnum
            An AminoAcidEnum object.
            
        Raises
        ------
        ValueError
            If the three-letter code is not a valid amino acid.
        """
        try:
            return cls[AMINO_ACID_THREE_LETTER_TOKEN_TO_NAME[three_letter_code.capitalize()]]
        except KeyError:
            raise ValueError(f'Invalid three-letter code for {cls.__name__}: {three_letter_code}')

    # Conversion methods unique to AminoAcidEnum
    def to_one_letter(self) -> str:
        """Convert an AminoAcidEnum object to a one-letter code.
        
        Returns
        -------
        str
            A one-letter code representing the amino acid.
        """
        return str(self)
    
    def to_three_letter(self) -> str:
        """Convert an AminoAcidEnum object to a three-letter code.
        
        Returns
        -------
        str
            A three-letter code representing the amino acid.
        """
        return NAME_TO_AMINO_ACID_THREE_LETTER_TOKEN[self.name]
    
    #endregion


class CodonEnum(BioToken):
    #region Methods that require implementation
    def __str__(self) -> str:
        return NAME_TO_TRIPLET_TOKEN[self.name]
    
    @classmethod
    def from_str(cls, s: str) -> 'CodonEnum':
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

    def is_standard(self) -> bool:
        """Return True if the codon is one of the 64 standard codon.
        
        Returns
        -------
        bool
            True if the codon is a standard codon.
        """
        return self.name in STANDARD_CODONS

    def is_any(self) -> bool:
        """Return True if the codon is NNN.
        
        Returns
        -------
        bool
            True if the codon is NNN, False otherwise.
        """
        return self.name == 'NNN'
    
    def is_ambiguous(self) -> bool:
        """Return True if the codon is ambiguous.
        
        Returns
        -------
        bool
            True if the codon is ambiguous, False otherwise.
        """
        return self.is_any()
    
    #endregion
    
    #region Override methods with default implementations
    
    # Factory methods
    @classmethod
    def from_int(cls, i: int) -> 'CodonEnum':
        """Convert an integer to a Codon object.
        
        Parameters
        ----------
        i : int
            The integer to convert.
            
        Returns
        -------
        Codon
            The codon corresponding to the given integer.
        """
        return cast('CodonEnum', super().from_int(i))

    @classmethod
    def from_onehot(cls, onehot: Sequence[int]) -> 'CodonEnum':
        """Convert a onehot vector to a Codon object.
        
        Parameters
        ----------
        onehot : Sequence[int]
            The one-hot vector.
            
        Returns
        -------
        Codon
            The codon corresponding to the given one-hot vector.
        """
        return cast('CodonEnum', super().from_onehot(onehot))

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
    
    #endregion
    
    #region Methods unique to CodonEnum
    # Factory methods
    @classmethod
    def from_nucleotides(cls, seq: Sequence[NucleotideEnum]) -> 'CodonEnum':
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
        """
        try:
            str_triplet = ''.join([str(n) for n in seq])
            return cls.from_str(str_triplet)
        except ValueError:
            raise ValueError(f'Invalid nucleotide sequence triplet: {seq}')

    @classmethod
    def from_nucleotide_onehot(cls, onehot: Sequence[Sequence[int]]) -> 'CodonEnum':
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
            return cls.from_nucleotides([NucleotideEnum.from_onehot(onehot[i]) for i in range(3)])
        except ValueError:
            raise ValueError(f'Invalid one-hot value for {cls.__name__}: {onehot}')

    @classmethod
    def start_codon(cls) -> 'CodonEnum':
        """Return the start codon.
        
        Returns
        -------
        Codon
            The start codon ATG.
        """
        return cls['ATG']

    # Conversion methods
    
    @property
    def nucleotide_class(self) -> Type[NucleotideEnum]:
        raise NotImplementedError("CodonEnum must be subclassed to implement this property")
    
    @property
    def aminoacid_class(self) -> Type[AminoAcidEnum]:
        raise NotImplementedError("CodonEnum must be subclassed to implement this property")
    
    
    def to_nucleotides(self) -> Sequence[NucleotideEnum]:
        """Convert the codon to a sequence of nucleotides.
        
        Returns
        -------
        Sequence[Nucleotide]
            The sequence of nucleotides.
        """
        return [self.nucleotide_class.from_str(n) for n in str(self)]
    
    def to_anticodon(self) -> 'CodonEnum':
        """Return the anticodon of the codon.
        
        Returns
        -------
        Codon
            Complimentary sequence of the codon.
        """
        try:
            anticodn_str = ''.join([COMPLEMENTARY_NUCLEOTIDES[n] for n in str(self)])
            return self.__class__[anticodn_str]
        except KeyError:
            return self.__class__(self.value)
    
    def translate(self) -> AminoAcidEnum:
        """Translate the codon to an amino acid.
        
        Returns
        -------
        AminoAcid
            The amino acid corresponding to the codon.
            
        Raises
        ------
        ValueError
            If the codon is not a valid codon.
        """
        if self.is_stop_codon():
            return self.aminoacid_class['Stop']
        try:
            return self.aminoacid_class[TRANSLATION_TABLE[self.name]]
        except KeyError:
            raise ValueError(f'Cannot be translated: {self}')

    # Properties
    def is_start_codon(self) -> bool:
        """Return True if the codon is the start codon.
        
        Returns
        -------
        bool
            True if the codon is the start codon, False otherwise.
        """
        return self.name == 'ATG'
    
    def is_stop_codon(self) -> bool:
        """Return True if the codon is a stop codon.

        Returns
        -------
        bool
            True if the codon is a stop codon, False otherwise.
        """
        return self.name in STOP_CODONS
    
    def is_twofold_degenerate(self) -> bool:
        """Return True if the codon is twofold degenerate.
        
        Returns
        -------
        bool
            True if the codon is twofold degenerate, False otherwise.
        """
        return self.name in TWOFOLD_DEGENERATE_CODONS
    
    def is_threefold_degenerate(self) -> bool:
        """Return True if the codon is threefold degenerate.
        
        Returns
        -------
        bool 
            True if the codon is threefold degenerate, False otherwise.
        """
        return self.name in THREEFOLD_DEGENERATE_CODONS
    
    def is_fourfold_degenerate(self) -> bool:
        """Return True if the codon is fourfold degenerate.
        
        Returns
        -------
        bool
            True if the codon is fourfold degenerate, False otherwise.
        """
        return self.name in FOURFOLD_DEGENERATE_CODONS

    def is_sixfold_degenerate(self) -> bool:
        """Return True if the codon is sixfold degenerate.
        
        Returns
        -------
        bool
            True if the codon is sixfold degenerate, False otherwise.
        """
        return self.name in SIXFOLD_DEGENERATE_CODONS

    #endregion
    