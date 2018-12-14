import doctest
import gzip
import logging
import random
import shutil

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

# 60 length seen to miss only 1 in 10,000
PRIMER_IN_READ_LIMIT = 60

# For consistency of sampling
random.seed('CrispyCrunch')


def in_fastq(fastq: str, primer_seq: str) -> Tuple[int, int]:
    """
    Counts lines of a fastq that contain a primer sequence and guide sequence
    in the expected locations.

    >>> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >>> primer_seq = 'CGAGGAGATACAGGCGGAG'

    >>> in_fastq(r1, primer_seq)
    (1232, 1213)
    >>> in_fastq(r1, reverse_complement(primer_seq))
    (1232, 0)

    >>> r1gz = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz'
    >>> in_fastq(r1gz, primer_seq)
    (1199, 1178)
    """
    seq_lines = _get_random_seq_lines(fastq)
    primer_matches = [line for line in seq_lines
                      if primer_seq in line[:PRIMER_IN_READ_LIMIT]]
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

    # logger.warning((primer_seq_fwd, fastq_r1.split('/')[-1], in_r1, fastq_r2.split('/')
    #                 [-1], in_r2, in_r1[1] + in_r1[1] > (in_r1[0] + in_r2[0]) * 0.25))

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
    ... 'guide_loc': 'chr7:4-23:-',
    ... 'primer_seq_fwd': 'CGAGGAGATACAGGCGGAG',
    ... 'primer_seq_rev': 'GTGGACGAGACGTGGTTAA'}]
    >>> find_matching_pairs(fastqs, records) == [fastqs]
    True

    >>> find_matching_pairs(fastqs, records, demultiplex=True)
    [('fastqs/demultiplexed/chr7_4-23_-_R1_.fastq', 'fastqs/demultiplexed/chr7_4-23_-_R2_.fastq')]

    >>> fastqs = ['input/CrispyCrunch/mNGplate3_unsorted_A10_TAF1B-C_S10_R1_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A10_TAF1B-C_S10_R2_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A11_TAF1C-C_S11_R1_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A11_TAF1C-C_S11_R2_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A12_TAF1D-N_S12_R1_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A12_TAF1D-N_S12_R2_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A1_POLR1A-C_S1_R1_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A1_POLR1A-C_S1_R2_001.fastq.gz']
    >>> records = [{
    ... 'guide_loc': 'chr7:4-23:-',
    ... 'primer_seq_fwd': 'TGTACTGTCACTTGGA',
    ... 'primer_seq_rev': 'CTCAACACCCTGACAC'}]
    >>> find_matching_pairs(fastqs, records)
    [('input/CrispyCrunch/mNGplate3_unsorted_A1_POLR1A-C_S1_R1_001.fastq.gz', 'input/CrispyCrunch/mNGplate3_unsorted_A1_POLR1A-C_S1_R2_001.fastq.gz')]

    >>> records.append(records[0])
    >>> find_matching_pairs(fastqs, records)
    Traceback (most recent call last):
    ...
    ValueError: Primers should be unique for matching to fastqs
    """
    if len({(row['primer_seq_fwd'], row['primer_seq_rev']) for row in records}) != len(list(records)):
        raise ValueError('Primers should be unique for matching to fastqs')

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
        primer_seq_fwd = row['primer_seq_fwd'].strip().upper()
        primer_seq_rev = row['primer_seq_rev'].strip().upper()
        match_key = (primer_seq_fwd, primer_seq_rev)
        logger.info('Finding matching pair of files for {} ...'.format(match_key))
        pair = find_matching_pair(
            [f for f in fastqs if '_R1_' in f and f not in seen],
            [f for f in fastqs if '_R2_' in f and f not in seen],
            primer_seq_fwd,
            primer_seq_rev,
            pool)
        if pair:
            logger.info('Found matching pair of files {} for {}.'.format(pair, match_key))
            pairs.append(pair)
            seen.add(pair[0])
            seen.add(pair[1])
        if match_key in match_keys:
            logger.warning('Duplicate detected: {}. Results may be unexpected.'.format(
                match_key))
        match_keys.add(match_key)

    if parallelize:
        pool.shutdown()

    return pairs


def _demultiplex(fastqs: Iterable,
                 records: Iterable[Mapping[str, str]],
                 parallelize: bool = False) -> Iterable[str]:
    """
    Split fastqs files into new files by prefix or suffix.

    >>> fastqs = ('fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq', 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq')
    >>> records = [{
    ... 'guide_loc': 'chr7:4-23:-',
    ... 'primer_seq_fwd': 'CGAGGAGATACAGGCGGAG',
    ... 'primer_seq_rev': 'GTGGACGAGACGTGGTTAA'}]

    >>> _demultiplex(fastqs, records)
    ['fastqs/demultiplexed/chr7_4-23_-_R1_.fastq', 'fastqs/demultiplexed/chr7_4-23_-_R2_.fastq']

    >>> _demultiplex(fastqs, records, True)
    ['fastqs/demultiplexed/chr7_4-23_-_R1_.fastq', 'fastqs/demultiplexed/chr7_4-23_-_R2_.fastq']

    >>> records.append(records[0])
    >>> _demultiplex(fastqs, records)
    Traceback (most recent call last):
    ...
    ValueError: primers must be unique for demultiplexing
    """
    if not (len({row['primer_seq_fwd'] for row in records}) == len(list(records))
            and len({row['primer_seq_rev'] for row in records}) == len(list(records))):
        raise ValueError('primers must be unique for demultiplexing')

    # See _get_demux_path
    demux_dir = (Path(list(fastqs)[0]).parent / 'demultiplexed')
    shutil.rmtree(demux_dir, ignore_errors=True)
    demux_dir.mkdir()

    if parallelize:
        pool = ProcessPoolExecutor()
    else:
        pool = None  # type: ignore

    new_paths = set()  # type: ignore
    if parallelize:
        sets = pool.map(_demux_fastq, fastqs, (records for f in fastqs))
        for s in sets:
            new_paths = new_paths.union(s)
    else:
        for fastq in fastqs:
            new_paths = new_paths.union(_demux_fastq(fastq, records))

    if parallelize:
        pool.shutdown(wait=True)

    return sorted(new_paths)


def _demux_fastq(fastq: str, records) -> set:
    logging.info('Demultiplexing fastq {}...'.format(fastq))
    discarded = 0
    total = 0
    new_paths = set()  # type: ignore
    for read in _get_reads(fastq):
        line = read[1]
        new_path = _get_demux_path(line, records, Path(fastq), '.fastq')
        if new_path:
            with open(new_path, 'at') as file:
                file.write('\n'.join(read) + '\n')
            new_paths.add(new_path)
        else:
            discarded += 1
        total += 1
    logger.info('{} reads out of {} in {} were discarded because they could not be matched by primer'.format(
        discarded, total, Path(fastq).name))
    return new_paths


def _get_demux_path(
    line: str,
    records: Iterable[Mapping[str, str]],
    old_path: Path,
    suffix: str = '.fastq.gz',
) -> str:
    matches = 0
    read_file_marker = ''
    new_path = ''
    for row in records:
        # TODO (gdingle): should PRIMER_IN_READ_LIMIT be stricter for demuxing?
        if row['primer_seq_fwd'] in line[:PRIMER_IN_READ_LIMIT]:
            read_file_marker = '_R1_'
            matches += 1

        if row['primer_seq_rev'] in line[:PRIMER_IN_READ_LIMIT]:
            read_file_marker = '_R2_'
            matches += 1

        if matches:
            guide_loc = row['guide_loc'].replace(':', '_')
            new_path = '{}/demultiplexed/{}{}{}'.format(
                old_path.parent, guide_loc, read_file_marker, suffix)

    if matches > 1:
        # TODO (gdingle): how to rank more than one match?
        # TODO (gdingle): what about bias towards last match?
        logger.warning('More than one match for read: ' + line)

    return new_path


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
            'Cannot find match for primers {} in {} candidate FastQ file pairs'.format(
                (primer_seq_fwd, primer_seq_rev), len(list(fastq_r1s)), ))


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


@lru_cache(maxsize=96 * 2)
def _get_random_seq_lines(fastq: str, random_fraction: float = 0.1) -> List[str]:
    file = gzip.open(fastq, 'rt') if fastq.endswith('.gz') else open(fastq)
    first_line = next(file)
    assert first_line.startswith('@'), 'Expecting fastq format, not: ' + first_line
    with file:
        # Every fourth line
        seq_lines = [line for i, line in enumerate(file)
                     if i % 4 == 0 and random.random() < random_fraction]
    return seq_lines


def _get_reads(fastq: str) -> Iterable[tuple]:
    file = gzip.open(fastq, 'rt') if fastq.endswith('.gz') else open(fastq)
    with file:
        while True:
            next_read = tuple(l.strip() for l in islice(file, 4))
            if not len(next_read) == 4:
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
