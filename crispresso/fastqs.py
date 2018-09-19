import doctest
import gzip

from pathlib import Path
from typing import Tuple

"""
Matches fastq files to designed guides and primers so we can avoid relying
on brittle file naming conventions or mutable sample sheets.

The matching works on the assumption that primer sequences appear always at the beginning
of sequence lines in a fastq, guides are somewhere following, and there will always be a "high"
number of such matches in a "matching" file.

An added benefit is validating fastqs before full alignment by Crispresso.
"""


def in_fastq(fastq: str, primer_seq: str, guide_seq: str) -> bool:
    """
    Determines whether a fastq contains a primer sequence and guide sequence in the expected locations.

    >>> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >>> primer_seq = 'CGAGGAGATACAGGCGGAG'
    >>> guide_seq = 'AATCGGTACAAGATGGCGGA'

    >>> in_fastq(r1, primer_seq, guide_seq)
    True
    >>> in_fastq(r1, guide_seq, primer_seq)
    False

    >>> r1gz = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz'
    >>> in_fastq(r1gz, primer_seq, guide_seq)
    True
    """
    file = gzip.open(fastq, 'rt') if fastq.endswith('.gz') else open(fastq)
    with file:
        seq_lines = [line for line in file.readlines()
                     if line[0] in 'ACGT']
    # TODO (gdingle): how long should primer_seq be?
    primer_matches = [line for line in seq_lines
                      if line.startswith(primer_seq)]
    # TODO (gdingle): should we always expect a full guide? what about break location?
    guide_matches = [line for line in primer_matches
                     if guide_seq in line[len(primer_seq):]]
    return (
        # TODO (gdingle): are these thresholds good?
        len(primer_matches) > len(seq_lines) * 0.5 and
        len(guide_matches) > len(primer_matches) * 0.9
    )


def matches_fastq_pair(
        fastq_r1: str,
        fastq_r2: str,
        primer_seq_fwd: str,
        primer_seq_rev: str,
        guide_seq: str) -> bool:
    """
    Determines whether a pair of fastq files, r1 and r2, contain the given primers and guide.

    >>> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >>> r2 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq'
    >>> primer_seq_fwd = 'CGAGGAGATACAGGCGGAG'
    >>> primer_seq_rev = 'GTGGACGAGACGTGGTTAA'
    >>> guide_seq = 'AATCGGTACAAGATGGCGGA'

    >>> matches_fastq_pair(r1, r2, primer_seq_fwd, primer_seq_rev, guide_seq)
    True
    >>> matches_fastq_pair(r2, r1, primer_seq_fwd, primer_seq_rev, guide_seq)
    False
    """
    in_r1 = in_fastq(fastq_r1, primer_seq_fwd, guide_seq)
    in_r2 = in_fastq(fastq_r2, primer_seq_rev, reverse_complement(guide_seq))
    return in_r1 and in_r2


def find_matching_pair(
        fastq_dir: str,
        primer_seq_fwd: str,
        primer_seq_rev: str,
        guide_seq: str,
        file_suffix='fastq') -> Tuple[str, ...]:
    """
    Find matching pair of fastq files in a dir based on primers and guide.

    >>> primer_seq_fwd = 'CGAGGAGATACAGGCGGAG'
    >>> primer_seq_rev = 'GTGGACGAGACGTGGTTAA'
    >>> guide_seq = 'AATCGGTACAAGATGGCGGA'
    >>> find_matching_pair('fastqs', primer_seq_fwd, primer_seq_rev, guide_seq)
    ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq')

    >>> find_matching_pair('fastqs', primer_seq_fwd, primer_seq_rev, guide_seq, 'fastq.gz')
    ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq.gz')
    """
    fastq_r1s = Path(fastq_dir).glob('*_R1_*.' + file_suffix)
    fastq_r2s = Path(fastq_dir).glob('*_R2_*.' + file_suffix)
    matches = [
        (str(r1), str(r2)) for r1, r2 in zip(fastq_r1s, fastq_r2s)
        if matches_fastq_pair(str(r1), str(r2), primer_seq_fwd, primer_seq_rev, guide_seq)]
    assert len(matches) <= 1, 'More than one match'
    if matches:
        return matches[0]
    else:
        return tuple()


def reverse_complement(seq: str) -> str:
    """
    >>> seq_in = 'AATCGGTACAAGATGGCGGA'
    >>> seq_out = 'TCCGCCATCTTGTACCGATT'
    >>> reverse_complement(seq_in) == seq_out
    True
    >>> reverse_complement(seq_out) == seq_in
    True
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[base] for base in reversed(seq))


if __name__ == '__main__':
    doctest.testmod()
