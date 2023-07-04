from enum import Enum


class MaskTokenMixin(Enum):
    def is_mask(self) -> bool:
        """Returns True if the token is a mask token.
        Expects the Mask enum name to be 'Mask'.
        
        Returns
        -------
        bool
            True if the token is a mask token, False otherwise.
        """
        return self.name == 'Mask'


class GapTokenMixin(Enum):
    def is_gap(self) -> bool:
        """Returns True if the token is a gap token.
        Expects the Gap enum name to be 'Gap'.
        
        Returns
        -------
        bool
            True if the token is a gap token, False otherwise.
        """
        return self.name == 'Gap'


class OtherTokenMixin(Enum):
    def is_other(self) -> bool:
        """Returns True if the token is an other token.
        Expects the Other enum name to be 'Other'.
        
        Returns
        -------
        bool
            True if the token is an other token, False otherwise.
        """
        return self.name == 'Other'
    
    
class SpecialTokenMixin(MaskTokenMixin, GapTokenMixin, OtherTokenMixin):
    def is_special(self) -> bool:
        """Returns True if the token is a special token.
        
        Returns
        -------
        bool
            True if the token is a special token, False otherwise.
        """
        return self.is_mask() or self.is_gap() or self.is_other()
