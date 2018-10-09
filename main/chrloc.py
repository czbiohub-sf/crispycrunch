import re


class ChrLoc:
    """
    Represents a valid chromosome location.

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
    """
    # See also CHR_REGEX in conversions.py
    CHR_REGEX = r'^chr([0-9XY]+):([0-9,]+)-([0-9,]+[0-9])$'

    max_length = 2000
    # TODO (gdingle): set to min of Crispor?
    min_length = 20

    def __init__(self, value: str) -> None:
        if not isinstance(value, str):
            raise ValueError('Input value must be a string')

        matches = re.match(self.CHR_REGEX, value.strip())
        if not matches:
            raise ValueError('Cannot parse chromosome location from "{}"'.format(value))

        self.chr = matches[1]  # remember chrX
        assert int(self.chr) in range(1, 24) or self.chr in ('X', 'Y')

        self.start = int(matches[2].replace(',', ''))
        self.end = int(matches[3].replace(',', ''))
        assert self.start <= self.end, (self.start, self.end)

        if len(self) > self.max_length:
            raise ValueError('{} is longer than the max length of {} for {}'.format(
                len(self), self.max_length, value))

        if len(self) < self.min_length:
            raise ValueError('{} is shorter than the min length of {} for {}.'.format(
                len(self), self.min_length, value))

    def __len__(self):
        return self.end - self.start + 1  # inclusive range

    def __str__(self):
        return 'chr{}:{}-{}'.format(
            self.chr, self.start, self.end)

    def __repr__(self):
        return f"ChrLoc('{self}')"

    def __eq__(self, other):
        return self.chr == other.chr and self.start == other.start and self.end == other.end


def get_guide_loc(target_loc: ChrLoc, guide_offset: int, guide_len=20
                  ) -> ChrLoc:
    """
    # TODO (gdingle): is this actually correct def for guide_loc?
    >>> get_guide_loc(ChrLoc('chr7:5569177-5569415'), 191)
    ChrLoc('chr7:5569348-5569367')
    """
    start = target_loc.start + guide_offset
    # TODO (gdingle): is this correct for reverse guides?
    # Guide goes backwards from pam, right to left
    # Minus one, for length inclusive
    return ChrLoc(
        'chr{}:{}-{}'.format(target_loc.chr, start - guide_len, start - 1)
    )


def get_primer_loc(
        primer_product: str,
        guide_seq: str,
        guide_loc: ChrLoc) -> ChrLoc:
    """
    Returns the chr loc of a primer product seq based on the known position of
    the guide within it.

    >>> get_primer_loc('''TGCTGGCTGGCCATTTCTAAACTTCCATTTGAATTTAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ... NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ... NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ... NNNNNNNNNNNNNNNNNNNNNNNTATTACTTTTGTCTTCTACTAGCCAAAAGAATGTCAACAGAAATCAGAACATAACAC
    ... TAAGTAAGTTTAACATGTACTTTTATTAACAACTTAATACAAGACTGTACACTGTAGGTGCTGAAATCAACCCACTCCT''',
    ... 'AATACAAGACTGTACACTGTAGG', ChrLoc('chr2:136114380-136114402'))
    ChrLoc('chr2:136114025-136114423')
    """
    primer_product = primer_product.replace('\n', '')
    # TODO (gdingle): time to start using seq objects? Biopython?
    # validate_seq(primer_product)
    # validate_seq(guide_seq)
    assert guide_seq in primer_product
    assert len(guide_loc) == len(guide_seq)

    primer_start = guide_loc.start - primer_product.index(guide_seq)
    primer_end = primer_start + len(primer_product) - 1
    assert primer_end - primer_start <= 500, 'Primers should always be less than 500 bp for paired reads'
    assert primer_end - primer_start == len(primer_product) - 1

    return ChrLoc(
        'chr{}:{}-{}'.format(guide_loc.chr, primer_start, primer_end)
    )


def get_guide_cut_to_insert(target_loc: ChrLoc, guide_loc: ChrLoc) -> int:
    """
    Based on https://czi.quip.com/YbAhAbOV4aXi/

    >>> get_guide_cut_to_insert(ChrLoc('chr5:134649077-134649174'), ChrLoc('chr5:134649061-134649080'))
    -2
    >>> get_guide_cut_to_insert(ChrLoc('chr5:134649077-134649174'), ChrLoc('chr5:134649077-134649096'))
    14
    """
    # insert location is assumed to be always one codon past start codon
    # TODO (gdingle): fix me for stop_codon
    insert_loc = target_loc.start + 3

    # cut location is assumed to be always in between the 3rd and 4th nucleotide
    # away from the PAM site
    cut_loc = guide_loc.end - 3
    # TODO (gdingle): is this correct for reverse guides?
    # Plus one for inclusive range
    return cut_loc - insert_loc + 1


if __name__ == '__main__':
    import doctest
    doctest.testmod()
