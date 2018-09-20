import doctest
import gzip

from functools import lru_cache
from pathlib import Path
from typing import Iterable, List, Tuple

"""
Matches fastq files to designed guides and primers so we can avoid relying on
brittle file naming conventions or mutable sample sheets.

The matching works on the assumption that primer sequences appear always at the
beginning of sequence lines in a fastq, guides are somewhere following, and
there will always be a "high" number of such matches in a "matching" file.

An added benefit is validating fastqs before full alignment by Crispresso.
"""


def in_fastq(fastq: str, primer_seq: str, guide_seq: str
             ) -> Tuple[int, int, int]:
    """
    Counts lines of a fastq that contain a primer sequence and guide sequence
    in the expected locations.

    >>> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >>> primer_seq = 'CGAGGAGATACAGGCGGAG'
    >>> guide_seq = 'AATCGGTACAAGATGGCGGA'

    >>> in_fastq(r1, primer_seq, guide_seq)
    (18897, 11490, 11019)
    >>> in_fastq(r1, guide_seq, primer_seq)
    (18897, 0, 0)

    >>> r1gz = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz'
    >>> in_fastq(r1gz, primer_seq, guide_seq)
    (18897, 11490, 11019)
    """
    seq_lines = _get_seq_lines(fastq)
    primer_matches = [line for line in seq_lines
                      if line.startswith(primer_seq)]
    # Only the first 16 chars have matches because of cut site
    # TODO (gdingle): verify CRISPR logic
    guide_matches = [line for line in primer_matches
                     if guide_seq[0:16] in line]

    return (len(seq_lines), len(primer_matches), len(guide_matches))


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

    return (
        # TODO (gdingle): are these thresholds good?
        in_r1[1] + in_r1[1] > (in_r1[0] + in_r2[0]) * 0.4 and
        # TODO (gdingle): why does guide_seq appear so infrequently sometimes?
        in_r1[2] + in_r2[2] > (in_r1[1] + in_r1[1]) * 0.0001
    )


def find_matching_pair_from_dir(
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
    >>> find_matching_pair_from_dir('fastqs', primer_seq_fwd, primer_seq_rev, guide_seq)
    ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq')

    >>> find_matching_pair_from_dir('fastqs', primer_seq_fwd, primer_seq_rev, guide_seq, 'fastq.gz')
    ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq.gz')
    """
    fastq_r1s = Path(fastq_dir).glob('*_R1_*.' + file_suffix)
    fastq_r2s = Path(fastq_dir).glob('*_R2_*.' + file_suffix)
    return find_matching_pair(fastq_r1s, fastq_r2s, primer_seq_fwd, primer_seq_rev, guide_seq)


def find_matching_pair(
        fastq_r1s: Iterable,
        fastq_r2s: Iterable,
        primer_seq_fwd: str,
        primer_seq_rev: str,
        guide_seq: str) -> Tuple[str, ...]:
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


@lru_cache(maxsize=None)
def _get_seq_lines(fastq: str) -> List[str]:
    file = gzip.open(fastq, 'rt') if fastq.endswith('.gz') else open(fastq)
    with file:
        seq_lines = [line for line in file if line[0] in 'ACGT']
    return seq_lines


if __name__ == '__main__':
    # print(reverse_complement('TTGCATAGGAAGTTCCCAAAGTACCAGTTTGCCACGGCATCAACTGCCCAGAAGGGAAGCGTGATGACAAAGAGGAGGTCGGCCACTGACAGGTGCAGCCTGTACTTGTCCGTCATGCTTCTCAGTTTCTTCTGGTAACCCATGACCAGGATGACCAATCCATTGCCCACAATGCCAGTTAAGAAGATGATGGAGTAGATGGTGGG'))
    doctest.testmod()
    # import cProfile
    # from pstats import SortKey
    # cProfile.run('doctest.testmod()', sort=SortKey.TIME)
