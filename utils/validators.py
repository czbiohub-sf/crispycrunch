"""
Custom validators of specialized types.
"""

import doctest
import re

from django.core.exceptions import ValidationError

try:
    from utils.chrloc import ChrLoc
except ImportError:
    from chrloc import ChrLoc  # type: ignore


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
    try:
        ChrLoc(value)
    except Exception:
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
    # skip optional "N" or "C" trailing
    value = value.split(',')[0].strip()
    if not any((
            is_chr(value),
            is_seq(value),
            is_ensemble_transcript(value),
            is_gene(value))):
        raise ValidationError(
            '"{}" is not a chromosome location or nucleic acid sequence or a Ensembl transcript ID or a HGNC gene name'.format(value))


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


def validate_unique_set(value: list) -> None:
    items = list(value)  # make copy
    for v in set(items):
        items.remove(v)
    if items:
        raise ValidationError('Duplicate targets: {}'.format(items))


if __name__ == '__main__':
    doctest.testmod()
