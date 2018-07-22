from abc import abstractproperty
from collections import OrderedDict
from functools import lru_cache; cache = lru_cache()
from typing import List

class AbstractPlateLayout:
    """
    Represents the grid layout of a well plate for biological research.
    """

    @abstractproperty
    def size(self):
        """Size of plate"""

    @abstractproperty
    def end_char(self):
        """Last character in sequence of Y-coordinates"""

    @abstractproperty
    def end_int(self):
        """Last integer in sequence of X-coordinates"""

    def __init__(self, contents: dict) -> None:
        self.chars = [chr(i) for i in range(ord('A'), ord(self.end_char) + 1)]
        self.ints = list(range(1, self.end_int + 1))
        assert len(self.chars) * len(self.ints) == self.size

        if not isinstance(contents, OrderedDict):
            contents = OrderedDict(contents.items())
        assert(len(contents) <= self.size)
        self.contents = contents

    @property
    @cache
    def well_positions(self) -> List[str]:
        """
        >>> pp = Plate96Layout({}).well_positions
        >>> pp[:3] == ['A1', 'A2', 'A3']
        True
        >>> pp[-3:] == ['H10', 'H11', 'H12']
        True
        >>> len(pp) == Plate96Layout.size
        True
        """
        return [c + str(i) for c in self.chars for i in self.ints]

    @property
    @cache
    def well_names(self) -> OrderedDict:
        """
        >>> pc = Plate96Layout({'xc': 'TCTACCTCTGTTGCACAGGC TGG'}).well_names
        >>> pc['A1'] == 'xc'
        True
        >>> pc['H12'] == None
        True
        >>> len(pc) == Plate96Layout.size
        True
        """
        names = list(self.contents.keys())
        return OrderedDict((pos, names[i] if i < len(names) else None)
                           for i, pos in enumerate(self.well_positions))

    @property
    @cache
    def well_seqs(self):
        """
        >>> pc = Plate96Layout({'xc': 'TCTACCTCTGTTGCACAGGC TGG'}).well_seqs
        >>> pc['A1'] == 'TCTACCTCTGTTGCACAGGC TGG'
        True
        >>> pc['H12'] == None
        True
        >>> len(pc) == Plate96Layout.size
        True
        """
        values = list(self.contents.values())
        return OrderedDict((pos, values[i] if i < len(values) else None)
                           for i, pos in enumerate(self.well_positions))


class Plate96Layout(AbstractPlateLayout):
    """
    96-well plate.
    """

    size = 96
    end_char = 'H'
    end_int = 12


class Plate384Layout(AbstractPlateLayout):
    """
    384-well plate.

    >>> pp = Plate384Layout({})
    >>> pp.size == len(pp.well_positions)
    True
    """

    size = 384
    end_char = 'P'
    end_int = 24


if __name__ == '__main__':
    import doctest
    doctest.testmod()
