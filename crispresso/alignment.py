import doctest

"""
Matches fastq files to designed guides and primers so we can avoid relying
on brittle file naming conventions or mutable sample sheets.

An added benefit is validating fastqs before full alignment by Crispresso.
"""


def in_fastq(fastq: str, primer_seq: str, guide_seq: str) -> bool:
    """
    Determines whether a fastq contains a primer sequence and guide sequence in the expected locations.

    # TODO (gdingle): should we always expect a full guide? what about break location?

    >> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >> primer_seq = 'CGAGGAGATACAGGCGGAGGGCGAGGAGATACAGGCGGAGGGCGAGGAGATACAGGCGGAGAGCG'
    >> guide_seq = 'AATCGGTACAAGATGGCGGA'
    >> in_fastq(r1, primer_seq, guide_seq)
    True
    >> in_fastq(r1, guide_seq, primer_seq)
    False
    """
    # # TODO (gdingle): full seq of guide? or subseq?
    pass


def matches_fastq_pair(
        fastq_r1: str,
        fastq_r2: str,
        primer_seq_fwd: str,
        primer_seq_rev: str,
        guide_seq: str) -> bool:
    """
    Determines whether a pair of fastq files, r1 and r2, contain the given primers and guide.

    >> r1 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq'
    >> r2 = 'fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq'
    >> primer_seq_fwd = 'CGAGGAGATACAGGCGGAGGGCGAGGAGATACAGGCGGAGGGCGAGGAGATACAGGCGGAGAGCG'
    # TODO (gdingle): verify these primers
    >> primer_seq_rev = 'GTGGACGAGACGTGGTTAACCGCGGCGCTTGGGTCGCTGGTCCGTCGCCGGCGCCACAGCCCCTG'
    >> guide_seq = 'AATCGGTACAAGATGGCGGA'
    >> in_fastq(r1, primer_seq, guide_seq)
    """
    # TODO (gdingle): reverse complement of guide?
    in_r1 = in_fastq(fastq_r1, primer_seq_fwd, guide_seq)
    in_r2 = in_fastq(fastq_r2, primer_seq_rev, reverse_complement(guide_seq))
    return in_r1 and in_r2


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
