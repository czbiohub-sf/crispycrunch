import doctest
import re

from django.core.exceptions import ValidationError

CHR_REGEX = r'^chr([0-9]+):([0-9,]+)-([0-9,]+[0-9])$'


def validate_seq(value: str) -> None:
    """
    See https://en.wikipedia.org/wiki/Nucleic_acid_sequence.

    >>> validate_seq('gtca')
    >>> validate_seq('asdf')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['asdf is not a nucleic acid sequence']
    """
    if re.match(r'^[ACGTRYKMSWBDHV]+$', value.upper()) is None:
        raise ValidationError('{} is not a nucleic acid sequence'.format(value))


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

    >>> validate_chr('chr1:11,130,540-11,130,751,')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['chr1:11,130,540-11,130,751, is not a chromosome location']

    >>> validate_chr('chrA:11,130,540-11,130,751')
    Traceback (most recent call last):
    ...
    django.core.exceptions.ValidationError: ['chrA:11,130,540-11,130,751 is not a chromosome location']
    """
    if re.match(CHR_REGEX, value) is None:
        raise ValidationError('{} is not a chromosome location'.format(value))


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
    django.core.exceptions.ValidationError: ['EENST00000330949 is not a Ensembl transcript ID']
    """
    if re.match(r'^ENST[0-9]+$', value) is None:
        raise ValidationError('{} is not a Ensembl transcript ID'.format(value))


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


def validate_chr_or_seq_or_enst(value: str) -> None:
    if not any((is_chr(value), is_seq(value), is_ensemble_transcript(value))):
        raise ValidationError(
            '{} is not a chromosome location or nucleic acid sequence or a Ensembl transcript ID'.format(value))


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


if __name__ == '__main__':
    doctest.testmod()
