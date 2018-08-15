import doctest
import re

from django.core.exceptions import ValidationError

# See also CHR_REGEX in conversions.py
CHR_REGEX = r'^chr([0-9XY]+):([0-9,]+)-([0-9,]+[0-9])$'
# See https://www.genenames.org/about/guidelines
# And see https://www.biostars.org/p/60118/ .
GENE_REGEX = r'^[A-Z0-9-]+$|^C[0-9XY]+orf[0-9]+$'


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
    if re.match(r'^[ACGTRYKMSWBDHV]+( [ACGTRYKMSWBDHV]{3})?$', value.upper()) is None:
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
            '{}" is not a chromosome location or nucleic acid sequence or a Ensembl transcript ID or a HGNC gene name'.format(value))


def get_guide_loc(target_loc: str, guide_offset: int, guide_len=20) -> str:
    """
    # TODO (gdingle): is this actually correct def for guide_loc?
    >>> get_guide_loc('chr7:5569177-5569415', 191)
    'chr7:5569368-5569388'
    """
    validate_chr(target_loc)
    matches = re.match(CHR_REGEX, target_loc)
    start = int(matches[2]) + guide_offset
    return 'chr{}:{}-{}'.format(matches[1], start, start + guide_len)


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


if __name__ == '__main__':
    doctest.testmod()
