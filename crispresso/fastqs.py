import doctest
import gzip
import logging

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache, partial
from itertools import islice
from pathlib import Path
from typing import Iterable, List, Mapping, Sequence, Set, Tuple

logger = logging.getLogger(__name__)

"""
Matches fastq files to designed guides and primers so we can avoid relying on
brittle file naming conventions or mutable sample sheets.

The matching works on the assumption that primer sequences appear always at the
beginning of sequence lines in a fastq, guides are somewhere following, and
there will always be a "high" number of such matches in a "matching" file.

An added benefit is validating fastqs before full alignment by Crispresso.

NOTE: vectorzied string matching with pandas is actually slower here.
See https://stackoverflow.com/questions/49112552/.

# TODO (gdingle): consider also
https://bergvca.github.io/2017/10/14/super-fast-string-matching.html
"""


def in_fastq(fastq: str, primer_seq: str) -> Tuple[int, int]:
    """
    Counts lines of a fastq that contain a primer sequence and guide sequence
    in the expected locations.

    >>> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >>> primer_seq = 'CGAGGAGATACAGGCGGAG'

    >>> in_fastq(r1, primer_seq)
    (12022, 11790)
    >>> in_fastq(r1, reverse_complement(primer_seq))
    (12022, 0)

    >>> r1gz = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz'
    >>> in_fastq(r1gz, primer_seq)
    (12022, 11790)
    """
    seq_lines = _get_seq_lines(fastq)
    primer_matches = [line for line in seq_lines
                      # 60 length seen to miss only 1 in 10,000
                      if primer_seq in line[:60]]
    return (len(seq_lines), len(primer_matches))


def matches_fastq_pair(
        primer_seq_fwd: str,
        primer_seq_rev: str,
        fastq_r1: str,
        fastq_r2: str) -> bool:
    """
    Determines whether a pair of fastq files, r1 and r2, contain the given primers and guide.

    >>> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >>> r2 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq'
    >>> primer_seq_fwd = 'CGAGGAGATACAGGCGGAG'
    >>> primer_seq_rev = 'GTGGACGAGACGTGGTTAA'

    >>> matches_fastq_pair(primer_seq_fwd, primer_seq_rev, r1, r2)
    True
    """
    fastq_r1 = str(fastq_r1)
    fastq_r2 = str(fastq_r2)
    assert fastq_r1.replace('_R1_', '') == fastq_r2.replace('_R2_', ''), \
        'FastQ filenames should match: {} {}'.format(fastq_r1, fastq_r2)

    in_r1 = in_fastq(fastq_r1, primer_seq_fwd)
    in_r2 = in_fastq(fastq_r2, primer_seq_rev)

    return (
        # The lowest seen so far has been 29% ... for a single correct file
        in_r1[1] + in_r1[1] > (in_r1[0] + in_r2[0]) * 0.25
    )


def find_matching_pair_from_dir(
        fastq_dir: str,
        primer_seq_fwd: str,
        primer_seq_rev: str,
        file_suffix='fastq') -> Tuple[str, ...]:
    """
    Find matching pair of fastq files in a dir based on primers and guide.

    >>> primer_seq_fwd = 'CGAGGAGATACAGGCGGAG'
    >>> primer_seq_rev = 'GTGGACGAGACGTGGTTAA'
    >>> find_matching_pair_from_dir('fastqs', primer_seq_fwd, primer_seq_rev)
    ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq')

    >>> find_matching_pair_from_dir('fastqs', primer_seq_fwd, primer_seq_rev, 'fastq.gz')
    ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq.gz')
    """
    fastq_r1s = list(Path(fastq_dir).glob('*_R1_*.' + file_suffix))
    fastq_r2s = list(Path(fastq_dir).glob('*_R2_*.' + file_suffix))
    return find_matching_pair(fastq_r1s, fastq_r2s, primer_seq_fwd, primer_seq_rev)


def find_matching_pairs(
        fastqs: Iterable,
        records: Iterable[Mapping[str, str]],
        parallelize: bool = False,
        demultiplex: bool = False,
) -> Sequence[Tuple[str, str]]:
    """
    >>> fastqs = ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq')
    >>> records = [{
    ... 'primer_seq_fwd': 'CGAGGAGATACAGGCGGAG',
    ... 'primer_seq_rev': 'GTGGACGAGACGTGGTTAA'}]
    >>> find_matching_pairs(fastqs, records) == [fastqs]
    True
    """
    seen: Set[str] = set()
    pairs: List[Tuple[str, str]] = []
    match_keys: Set[tuple] = set()

    if parallelize:
        pool = ProcessPoolExecutor()
    else:
        pool = None  # type: ignore

    if demultiplex:
        fastqs = _demultiplex(fastqs, records)

    for row in records:
        pair = find_matching_pair(
            [f for f in fastqs if '_R1_' in f and f not in seen],
            [f for f in fastqs if '_R2_' in f and f not in seen],
            row['primer_seq_fwd'].strip().upper(),
            row['primer_seq_rev'].strip().upper(),
            pool)
        if pair:
            pairs.append(pair)
            seen.add(pair[0])
            seen.add(pair[1])
        match_key = (row['primer_seq_fwd'], row['primer_seq_rev'])
        if match_key in match_keys:
            logger.warning('Duplicate detected: {}. Results may be unexpected.'.format(
                match_key))
        match_keys.add(match_key)

    if parallelize:
        pool.shutdown()

    return pairs


def _demultiplex(fastqs: Iterable,
                 records: Iterable[Mapping[str, str]],) -> Iterable[str]:
    """
    Split fastqs files into new files by prefix or suffix.

    >>> fastqs = ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq')
    >>> records = [{
    ... 'primer_seq_fwd': 'CGAGGAGATACAGGCGGAG',
    ... 'primer_seq_rev': 'GTGGACGAGACGTGGTTAA'}]
    >>> _demultiplex(fastqs, records)
    ['A1-ATL2-N-sorted-180212_S1_L001_R1_001-demultiplexed.fastq', 'A1-ATL2-N-sorted-180212_S1_L001_R2_001-demultiplexed.fastq']
    """
    new_fastqs = defaultdict(list)  # type: ignore
    for fastq in fastqs:
        path = Path(fastq)
        new_path = path.stem + '-demultiplexed' + path.suffix

        for read in _get_reads(fastq):
            # TODO (gdingle): do some work here
            new_fastqs[new_path].append(read)

    return list(new_fastqs.keys())


def find_matching_pair(
        fastq_r1s: Iterable,
        fastq_r2s: Iterable,
        primer_seq_fwd: str,
        primer_seq_rev: str,
        pool: ProcessPoolExecutor = None) -> Tuple[str, str]:

    if pool:
        bools = pool.map(
            partial(matches_fastq_pair, primer_seq_fwd, primer_seq_rev),
            fastq_r1s,
            fastq_r2s)
        matches = [
            (str(r1), str(r2)) for r1, r2, is_match in zip(fastq_r1s, fastq_r2s, bools)
            if is_match]
    else:
        matches = [
            (str(r1), str(r2)) for r1, r2 in zip(fastq_r1s, fastq_r2s)
            if matches_fastq_pair(primer_seq_fwd, primer_seq_rev, r1, r2)]

    if matches:
        if len(matches) > 1:
            logger.warning('More than one match: {}'.format(matches))
            # Return the first match on the assumption that inputs rows and files are
            # ordered similarly.
            # TODO (gdingle): deal with multi matches better
        return matches[0]
    else:
        raise ValueError(
            'Cannot find matching pair in {} candidate FastQ files for primer_seq_fwd {}'.format(
                len(list(fastq_r1s)), primer_seq_fwd))


# TODO (gdingle): move to conversions or somewhere?
def reverse_complement(seq: str) -> str:
    """
    >>> seq_in = 'AATCGGTACAAGATGGCGGA'
    >>> seq_out = 'TCCGCCATCTTGTACCGATT'
    >>> reverse_complement(seq_in) == seq_out
    True
    >>> reverse_complement(seq_out) == seq_in
    True
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return ''.join(complement[base] for base in reversed(seq))


@lru_cache(maxsize=1024)
def _get_seq_lines(fastq: str) -> List[str]:
    file = gzip.open(fastq, 'rt') if fastq.endswith('.gz') else open(fastq)
    first_line = next(file)
    assert first_line.startswith('@'), 'Expecting fastq format, not: ' + first_line
    with file:
        # Every fourth line
        seq_lines = [line for i, line in enumerate(file) if i % 4 == 0]
    return seq_lines


@lru_cache(maxsize=1024)
def _get_reads(fastq: str) -> Iterable[tuple]:
    file = gzip.open(fastq, 'rt') if fastq.endswith('.gz') else open(fastq)
    with file:
        while True:
            next_read = tuple(l.strip() for l in islice(file, 4))
            if not len(next_read) == 4:
                # TODO (gdingle): how does fastq file end?
                break
            # See https://en.wikipedia.org/wiki/FASTQ_format
            assert next_read[0].startswith('@'), next_read[0]
            assert next_read[1].startswith(tuple('AGCT')), next_read[1]
            assert next_read[2].startswith('+'), next_read[2]
            yield next_read


if __name__ == '__main__':
    doctest.testmod(optionflags=doctest.FAIL_FAST)
    # print(reverse_complement('CGGGCAGCGGGTCCATCGCG'))

    # import timeit

    # fastqs = (
    #     'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq',
    #     'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq',
    #     'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R1_001.fastq',
    #     'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R2_001.fastq',
    #     # 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq',
    #     # 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq',
    #     # 'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R1_001.fastq',
    #     # 'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R2_001.fastq',

    #     # 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq',
    #     # 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq',
    #     # 'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R1_001.fastq',
    #     # 'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R2_001.fastq',
    #     # 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq',
    #     # 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq',
    #     # 'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R1_001.fastq',
    #     # 'fastqs/C12-CLTA-N-sorted-180212_S36_L001_R2_001.fastq',

    # )
    # records = [{
    #     'primer_seq_fwd': 'CGAGGAGATACAGGCGGAG',
    #     'primer_seq_rev': 'GTGGACGAGACGTGGTTAA'}]
    # find_matching_pairs(fastqs, records) == [fastqs]

    # out = timeit.timeit(
    #     """find_matching_pairs(fastqs, records)""",
    #     number=100,
    #     globals=globals()
    # )
    # print(out, 's')

    # out = timeit.timeit(
    #     f"""find_matching_pairs(fastqs, records, True)""",
    #     number=100,
    #     globals=globals()
    # )
    # print(out, 's')
