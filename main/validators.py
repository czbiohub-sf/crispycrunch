"""
Custom validators of crispycrunch specific types,
plus some util functions that operate on those types.

# TODO (gdingle): move util functions, or refactor all into objects
"""

import doctest
import re

from django.core.exceptions import ValidationError

# See also CHR_REGEX in conversions.py
CHR_REGEX = r'^chr([0-9XY]+):([0-9,]+)-([0-9,]+[0-9])$'
# See https://www.genenames.org/about/guidelines
# And see https://www.biostars.org/p/60118/ .
GENE_REGEX = r'^[A-Z0-9-]+$|^C[0-9XY]+orf[0-9]+$'


def validate_fastq(filename: str) -> None:
    """
    >>> validate_fastq('crispresso/fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq')
    >>> validate_fastq('crispresso/fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz')
    >>> validate_fastq('crispresso/fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fa')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['"crispresso/fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fa" is not a valid fastq file']
    """
    if (not filename.endswith('.fastq') and not filename.endswith('.fastq.gz')):
        raise ValidationError('"{}" is not a valid fastq file'.format(filename))


def validate_seq(value: str) -> None:
    """
    See https://en.wikipedia.org/wiki/Nucleic_acid_sequence.

    Also allows a trailing space separated PAM.

    >>> validate_seq('gtca')
    >>> validate_seq('CACTGCAACCTTGGCCTCCC GGG')
    >>> validate_seq('asdf')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['"asdf" is not a nucleic acid sequence']
    """
    if re.match(r'^[ACGTRYKMSWBDHVN]+( [ACGTRYKMSWBDHVN]{3})?$', value.upper()) is None:
        raise ValidationError('"{}" is not a nucleic acid sequence'.format(value))


def is_seq(value: str) -> bool:
    """
    >>> is_seq('gtca')
    True
    >>> is_seq('asdf')
    False
    """
    try:
        validate_seq(value)
    except ValidationError:
        return False
    else:
        return True


def validate_chr(value: str) -> None:
    """
    >>> validate_chr('chr1:11,130,540-11,130,751')

    >>> validate_chr('chrX:153701031-153701090')

    >>> validate_chr('chr1:11,130,540-11,130,751,')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['"chr1:11,130,540-11,130,751," is not a chromosome location']

    >>> validate_chr('chrA:11,130,540-11,130,751')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['"chrA:11,130,540-11,130,751" is not a chromosome location']
    """
    if re.match(CHR_REGEX, value) is None:
        raise ValidationError('"{}" is not a chromosome location'.format(value))


def validate_chr_length(value: str, max_length: int = 2000, min_length=23) -> None:
    """
    >>> validate_chr_length('chr1:11,130,540-11,130,751')
    >>> validate_chr_length('chr1:11,230,540-11,130,751')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['99789 is longer than the max length of 2000 for chr1:11,230,540-11,130,751']
    """
    matches = [m.replace(',', '') for m in re.match(CHR_REGEX, value).groups()]  # type: ignore
    length = abs(int(matches[1]) - sum(int(m) for m in matches[2:]))
    if length > max_length:
        raise ValidationError('{} is longer than the max length of {} for {}'.format(
            length, max_length, value))
    if length < min_length:
        raise ValidationError('{} is shorter than the min length of {} for {}.'.format(
            length, min_length, value))
    return None


def is_chr(value: str) -> bool:
    """
    >>> is_chr('chr1:11,130,540-11,130,751')
    True
    >>> is_chr('gtca')
    False
    """
    try:
        validate_chr(value)
    except ValidationError:
        return False
    return True


def validate_ensemble_transcript(value: str) -> None:
    """
    >>> validate_ensemble_transcript('EENST00000330949')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['"EENST00000330949" is not a Ensembl transcript ID']
    """
    if re.match(r'^ENST[0-9]+$', value) is None:
        raise ValidationError('"{}" is not a Ensembl transcript ID'.format(value))


def is_ensemble_transcript(value: str) -> bool:
    """
    >>> is_ensemble_transcript('ENST00000330949')
    True
    >>> is_ensemble_transcript('ENST00000398844')
    True
    >>> is_ensemble_transcript('EENST00000330949')
    False
    """
    try:
        validate_ensemble_transcript(value)
    except ValidationError:
        return False
    else:
        return True


def validate_chr_or_seq_or_enst_or_gene(value: str) -> None:
    if not any((
            is_chr(value),
            is_seq(value),
            is_ensemble_transcript(value),
            is_gene(value))):
        raise ValidationError(
            '"{}" is not a chromosome location or nucleic acid sequence or a Ensembl transcript ID or a HGNC gene name'.format(value))
    if is_chr(value):
        validate_chr_length(value)
        # TODO (gdingle): length of other types


def validate_gene(value: str) -> None:
    """
    >>> validate_gene('ATL2') is None
    True
    >>> validate_gene('atl2')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['"atl2" is not a valid HGNC gene name']
    """
    if re.match(GENE_REGEX, value) is None:
        raise ValidationError('"{}" is not a valid HGNC gene name'.format(value))


def is_gene(value: str) -> bool:
    """
    >>> is_gene('ATL3')
    True
    >>> is_gene('ATL_3')
    False
    """
    try:
        validate_gene(value)
    except ValidationError:
        return False
    else:
        return True


def validate_num_wells(value: dict, max: int = 96) -> None:
    total = sum(len(seqs) for seqs in value.values())
    if total > max:
        raise ValidationError(
            '{} items do not fit in a 96-well plate'.format(total))


def get_guide_loc(target_loc: str, guide_offset: int, guide_len=20) -> str:
    """
    # TODO (gdingle): is this actually correct def for guide_loc?
    >>> get_guide_loc('chr7:5569177-5569415', 191)
    'chr7:5569348-5569367'
    """
    validate_chr(target_loc)
    matches = re.match(CHR_REGEX, target_loc).groups()  # type: ignore
    start = int(matches[1]) + guide_offset
    # TODO (gdingle): is this correct for reverse guides?
    # Guide goes backwards from pam, right to left
    # Minus one, for length inclusive
    return 'chr{}:{}-{}'.format(matches[0], start - guide_len, start - 1)


def get_guide_cut_to_insert(target_loc: str, guide_loc: str) -> int:
    """
    Based on https://czi.quip.com/YbAhAbOV4aXi/

    >>> get_guide_cut_to_insert('chr5:134649077-134649174', 'chr5:134649061-134649080')
    -2
    >>> get_guide_cut_to_insert('chr5:134649077-134649174', 'chr5:134649077-134649096')
    14
    """
    validate_chr(target_loc)
    validate_chr(guide_loc)
    target_matches = re.match(CHR_REGEX, target_loc).groups()  # type: ignore
    # insert location is assumed to be always one codon past start codon
    insert_loc = int(target_matches[1]) + 3

    guide_matches = re.match(CHR_REGEX, guide_loc).groups()  # type: ignore
    # cut location is assumed to be always in between the 3rd and 4th nucleotide
    # away from the PAM site
    cut_loc = int(guide_matches[2]) - 3
    # TODO (gdingle): is this correct for reverse guides?
    # Plus one for inclusive range
    return cut_loc - insert_loc + 1


def get_hdr_template(target_seq: str, hdr_seq: str, hdr_tag: str = 'start_codon') -> str:
    """
    Inserts HDR sequencee in correct position in target sequence.
    Based on https://czi.quip.com/YbAhAbOV4aXi/.

    >>> get_hdr_template('ATGTCCCAGCCGGGAAT', 'NNN')
    'ATGNNNTCCCAGCCGGGAAT'
    """
    validate_seq(target_seq)
    validate_seq(hdr_seq)
    # TODO (gdingle): stop_codon
    assert hdr_tag == 'start_codon', 'stop_codon not implemented'
    first_codon = target_seq[0:3]
    assert first_codon == 'ATG'
    return first_codon + hdr_seq + target_seq[3:]


def get_hdr_primer(primer_product: str, hdr_seq: str, hdr_tag: str = 'start_codon') -> str:
    """
    Locates target codon in primer product then inserts HDR sequence.
    >>> get_hdr_primer('ATGTCCCAGCCGGGAAT', 'NNN')
    'ATGNNNTCCCAGCCGGGAAT'
    """
    validate_seq(primer_product)
    validate_seq(hdr_seq)
    # TODO (gdingle): stop_codon
    assert hdr_tag == 'start_codon', 'stop_codon not implemented'
    codon_index = primer_product.find('ATG')
    if codon_index == -1:
        # TODO (gdingle): good return value?
        return 'start_codon not found'
    assert primer_product[codon_index:codon_index + 3] == 'ATG'

    return primer_product[:codon_index] + \
        get_hdr_template(primer_product[codon_index:], hdr_seq)


def get_primer_loc(primer_product: str, guide_seq: str, guide_loc: str) -> str:
    """
    Returns the chr loc of a primer product seq based on the known position of
    the guide within it.

    >>> get_primer_loc('''TGCTGGCTGGCCATTTCTAAACTTCCATTTGAATTTAANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ... NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ... NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
    ... NNNNNNNNNNNNNNNNNNNNNNNTATTACTTTTGTCTTCTACTAGCCAAAAGAATGTCAACAGAAATCAGAACATAACAC
    ... TAAGTAAGTTTAACATGTACTTTTATTAACAACTTAATACAAGACTGTACACTGTAGGTGCTGAAATCAACCCACTCCT''',
    ... 'AATACAAGACTGTACACTGTAGG', 'chr2:136114380-136114402')
    'chr2:136114025-136114423'
    """
    primer_product = primer_product.replace('\n', '')
    # TODO (gdingle): time to start using seq objects? Biopython?
    validate_seq(primer_product)
    validate_seq(guide_seq)
    validate_chr(guide_loc)
    assert guide_seq in primer_product
    chr_num, start, end = re.match(CHR_REGEX, guide_loc).groups()  # type: ignore
    start, end = int(start), int(end)
    assert end - start == len(guide_seq) - 1  # inclusive range

    primer_start = start - primer_product.index(guide_seq)
    primer_end = primer_start + len(primer_product) - 1
    assert primer_end - primer_start <= 500, 'Primers should always be less than 500 bp for paired reads'
    assert primer_end - primer_start == len(primer_product) - 1

    return 'chr{}:{}-{}'.format(chr_num, primer_start, primer_end)


def mutate_guide_seq(guide_seq: str) -> str:
    """
    Silently mutates input sequence by substituing a different
    codon that encodes the same amino acid wherever possible.

    The input is assumed to start with a codon of 3bp.

    Based on https://czi.quip.com/YbAhAbOV4aXi/.

    Data from http://biopython.org/DIST/docs/api/Bio.SeqUtils.CodonUsage-pysrc.html

    >>> mutate_guide_seq('TGTTGCGATGAC')
    TGCTGTGACGAT
    """
    synonymous = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN': ['AAC', 'AAT'],
        'PRO': ['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'STOP': ['TAG', 'TGA', 'TAA'],
        'THR': ['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE': ['TTT', 'TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY': ['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE': ['ATC', 'ATA', 'ATT'],
        'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
        'HIS': ['CAT', 'CAC'],
        'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP': ['TGG'],
        'VAL': ['GTA', 'GTC', 'GTG', 'GTT'],
        'GLU': ['GAG', 'GAA'],
        'TYR': ['TAT', 'TAC'],
    }
    synonymous_index = dict(
        (codon, aa)
        for aa, codons in synonymous.items()
        for codon in codons
    )
    validate_seq(guide_seq)
    new_guide = ''
    for i in range(0, len(guide_seq), 3):
        codon = guide_seq[i:i + 3].upper()
        syns = set(synonymous[synonymous_index[codon]])
        if len(syns) > 1:
            syns.remove(codon)
        # TODO (gdingle): better to choose random syn?
        new_guide += syns.pop().lower()
    assert len(new_guide) == len(guide_seq)
    assert new_guide != guide_seq
    return new_guide


if __name__ == '__main__':
    doctest.testmod()
