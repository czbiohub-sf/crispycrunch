import functools
import re


@functools.total_ordering
class ChrLoc:
    """
    Represents a valid chromosome location for crispycrunch.

    `start` and `end` are in the same coordinate system as the input:
    one-based, inclusive.

    >>> chr_loc = ChrLoc('chr7:5569177-5569415')
    >>> str(chr_loc)
    'chr7:5569177-5569415'
    >>> repr(chr_loc)
    "ChrLoc('chr7:5569177-5569415')"

    >>> chr_loc = ChrLoc('asdf')
    Traceback (most recent call last):
    ...
    ValueError: Cannot parse chromosome location from "asdf"

    >>> ChrLoc('chr1:11,130,540-11,230,751')
    Traceback (most recent call last):
    ...
    ValueError: 100212 is longer than the max length of 2000 for chr1:11,130,540-11,230,751

    >>> ChrLoc('chr7:5569177-5569178')
    Traceback (most recent call last):
    ...
    ValueError: 2 is shorter than the min length of 20 for chr7:5569177-5569178.

    >>> ChrLoc('chr7:0-1000')
    Traceback (most recent call last):
    ...
    ValueError: Position must start at one not zero: "chr7:0-1000"

    Strands.

    >>> ChrLoc('chr7:5569177-5569415:+')
    ChrLoc('chr7:5569177-5569415:+')
    >>> ChrLoc('chr7:5569177-5569415:-')
    ChrLoc('chr7:5569177-5569415:-')

    >>> ChrLoc('chr7:5569177-5569415:+').as_strand_direction
    'chr7:5569177-5569415'
    >>> ChrLoc('chr7:5569177-5569415:-').as_strand_direction
    'chr7:5569415-5569177'
    """
    # See also CHR_REGEX in conversions.py
    CHR_REGEX = r'^chr([0-9XY]+):([0-9,]+)-([0-9,]+[0-9])(:[+\-1])?$'

    max_length = 2000
    min_length = 20

    def __init__(self, value: str) -> None:
        if not isinstance(value, str):
            raise ValueError('Input value must be a string')

        matches = re.match(self.CHR_REGEX, value.strip())
        if not matches:
            raise ValueError('Cannot parse chromosome location from "{}"'.format(value))

        self.chr = matches[1]  # remember chrX
        assert self.chr in ('X', 'Y') or int(self.chr) in range(1, 23)

        self.start = int(matches[2].replace(',', ''))
        self.end = int(matches[3].replace(',', ''))
        if self.start == 0 or self.end == 0:
            raise ValueError('Position must start at one not zero: "{}"'.format(value))

        assert self.start <= self.end, (self.start, self.end)

        if len(self) > self.max_length:
            raise ValueError('{} is longer than the max length of {} for {}'.format(
                len(self), self.max_length, value))

        if len(self) < self.min_length:
            raise ValueError('{} is shorter than the min length of {} for {}.'.format(
                len(self), self.min_length, value))

        self.strand = matches[4][1] if matches.group(4) else ''
        assert self.strand in ('+', '-', '')

    def __len__(self):
        return self.end - self.start + 1  # inclusive range

    def __str__(self):
        return 'chr{}:{}-{}{}'.format(
            self.chr, self.start, self.end,
            ':' + self.strand if self.strand else '')

    @property
    def opposite_strand(self):
        return {'+': '-', '-': '+'}[self.strand]

    @property
    def as_strand_direction(self):
        return 'chr{}:{}-{}'.format(
            self.chr,
            self.start if self.strand == '+' else self.end,
            self.end if self.strand == '+' else self.start,
        )

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return f"ChrLoc('{self}')"

    def __eq__(self, other):
        """
        >>> ChrLoc('chr5:1-20') == ChrLoc('chr5:1-20')
        True
        >>> ChrLoc('chr5:1-40') == ChrLoc('chr5:1-20')
        False
        >>> ChrLoc('chr5:1-40') == 'chr5:1-40'
        True
        >>> ChrLoc('chr5:1-20:-') == ChrLoc('chr5:1-20:-')
        True
        >>> ChrLoc('chr5:1-20') == ChrLoc('chr5:1-20:-')
        True
        >>> ChrLoc('chr5:1-20:+') == ChrLoc('chr5:1-20:-')
        False
        """
        if isinstance(other, str):
            # TODO (gdingle): wise idea to promote on comparison?
            other = ChrLoc(other)
        return (self.chr == other.chr and
                self.start == other.start and
                self.end == other.end
                # Only compare if both have a strand specified
                and (self.strand == other.strand  # noqa
                     if self.strand and other.strand
                     else True))

    def __lt__(self, other):
        """
        # TODO (gdingle): does this make sense?

        >>> ChrLoc('chr5:1-20') < ChrLoc('chr5:1-20')
        False
        >>> ChrLoc('chr5:1-20') < ChrLoc('chr5:2-21')
        True
        >>> ChrLoc('chr4:1-20') < ChrLoc('chr5:1-20')
        True
        """
        if isinstance(other, str):
            # TODO (gdingle): wise idea to promote on comparison?
            other = ChrLoc(other)
        return (self.chr < other.chr or
                self.start < other.start)

    def __contains__(self, other):
        """
        >>> ChrLoc('chr5:1-20') in ChrLoc('chr5:1-40')
        True
        >>> ChrLoc('chr5:1-40') in ChrLoc('chr5:1-20')
        False
        """
        if isinstance(other, str):
            # TODO (gdingle): wise idea to promote on comparison?
            other = ChrLoc(other)
        return (self.chr == other.chr and
                self.start <= other.start and
                self.end >= other.end)

    def copy(self, **kwargs) -> 'ChrLoc':
        """
        >>> chr_loc = ChrLoc('chr7:5569177-5569415')
        >>> chr_loc.copy(chr=14)
        ChrLoc('chr14:5569177-5569415')
        """
        chr_loc = ChrLoc(str(self))
        for k, v in kwargs.items():
            setattr(chr_loc, k, v)
        return chr_loc


class GuideChrLoc(ChrLoc):
    """
    Represents structure of guide accounting for strand.
    Based on https://czi.quip.com/YbAhAbOV4aXi/.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        assert len(self) == 20, (args, len(self))
        assert self.strand

    @property
    def pam(self) -> ChrLoc:
        """
        PAM sequence is adjacent to the 20bp protospacer guide.
        """
        if self.strand == '-':
            return self.copy(
                start=self.start - 3,
                end=self.start - 1
            )
        else:
            return self.copy(
                start=self.end + 1,
                end=self.end + 3
            )

    @property
    def cut(self):
        """
        Cut location is assumed to be always in between the 3rd and 4th
        nucleotide in the guide. ChrLoc always starts at 1.
        """
        if self.strand == '-':
            return self.copy(
                start=self.start + 2,
                end=self.start + 3
            )
        else:
            return self.copy(
                start=self.end - 3,
                end=self.end - 2
            )


def get_guide_loc(
        target_loc: ChrLoc,
        guide_offset: int,
        guide_len: int = 20,
        guide_direction: str = '+') -> GuideChrLoc:
    """
    Example forward:
    AAGATAGGTGATGAAGGAGGGTCCCCAGG
          GGTGATGAAGGAGGGTCCCC

    >>> get_guide_loc(ChrLoc('chr7:1-30'), 26)
    ChrLoc('chr7:7-26:+')

    Example reverse:
    CCATGGCTGAGCTGGATCCGTTCGGC
       TGGCTGAGCTGGATCCG

    >>> get_guide_loc(ChrLoc('chr7:1-26'), 0, 20, '-')
    ChrLoc('chr7:4-23:-')
    """
    pam = target_loc.start + guide_offset
    return GuideChrLoc(
        'chr{}:{}-{}:{}'.format(
            target_loc.chr,
            pam + 3 if guide_direction == '-' else pam - guide_len,
            pam + guide_len + 2 if guide_direction == '-' else pam - 1,
            guide_direction)
    )


def get_insert(
        target_loc: ChrLoc,
        hdr_tag: str = 'start_codon') -> int:
    """
    Get the desired insert point assuming the start or stop codon
    is in the expected place. See protospacex._get_start_end.

    # TODO (gdingle): we shouldn't need to recompute this... refactor with protospacex

    >>> get_insert(ChrLoc('chr5:1-30'))
    4
    >>> get_insert(ChrLoc('chr5:1-60'))
    4
    >>> get_insert(ChrLoc('chr5:1-30'), 'stop_codon')
    13
    >>> get_insert(ChrLoc('chr5:1-60'), 'stop_codon')
    28
    """
    # Insert happens to the left of the integer position.
    if hdr_tag == 'start_codon':
        # Insert position is assumed to be always one codon in.
        return target_loc.start + 3
    elif hdr_tag == 'stop_codon':
        assert len(target_loc) % 2 == 0, len(target_loc)
        mid = len(target_loc) // 2
        return target_loc.start + mid - 3
    else:
        assert False


def get_guide_cut_to_insert(
        target_loc: ChrLoc,
        guide_loc: GuideChrLoc,
        hdr_tag: str = 'start_codon') -> int:
    """
    Based on example from https://czi.quip.com/YbAhAbOV4aXi/

    >>> get_guide_cut_to_insert(ChrLoc('chr5:1-40'),
    ... GuideChrLoc('chr5:1-20:+'))
    14

    >>> get_guide_cut_to_insert(ChrLoc('chr5:1-40'),
    ... GuideChrLoc('chr5:2-21:-'))
    1

    # 30bp target, stop codon ends at 16, insert at 13
    # 20bp guide, end of cut is at 18
    >>> get_guide_cut_to_insert(ChrLoc('chr5:1-30'),
    ... GuideChrLoc('chr5:1-20:+'), 'stop_codon')
    5
    """
    assert guide_loc in target_loc
    cut = guide_loc.cut.end

    insert = get_insert(target_loc, hdr_tag)
    return cut - insert


def get_primer_loc(
        primer_product: str,
        guide_seq: str,
        guide_loc: GuideChrLoc) -> ChrLoc:
    """
    Returns the chr loc of a primer product seq based on the known position of
    the guide within it.

    >>> get_primer_loc('NNNAATACAAGACTGTACACTGTNNN',
    ... 'AATACAAGACTGTACACTGT', GuideChrLoc('chr2:11-30:+'))
    ChrLoc('chr2:8-33')
    """
    assert guide_seq in primer_product
    assert len(guide_loc) == len(guide_seq)

    primer_start = guide_loc.start - primer_product.index(guide_seq)
    primer_end = primer_start + len(primer_product) - 1
    assert primer_end - primer_start == len(primer_product) - 1

    return ChrLoc(
        'chr{}:{}-{}'.format(guide_loc.chr, primer_start, primer_end)
    )


if __name__ == '__main__':
    import doctest
    doctest.testmod()
