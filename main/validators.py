import re

from django.core.exceptions import ValidationError


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
    if re.match(r'^chr[0-9]+:[0-9,-]*[0-9]$', value) is None:
        raise ValidationError('{} is not a chromosome location'.format(value))


def validate_chr_or_seq(value: str) -> None:
    try:
        validate_seq(value)
    except ValidationError:
        try:
            validate_chr(value)
        except ValidationError:
            raise ValidationError('{} is not a chromosome location or nucleic acid sequence'.format(value))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
