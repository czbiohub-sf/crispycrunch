import doctest

import pysam

# TODO (gdingle): download hg38.fa and mm10.fa from public repo in docker build


def get_reference_amplicon(chr_loc, genome='hg38'):
    """
    >>> get_reference_amplicon('chr7:5569176-5569415')
    'ctgtctcaaaaaaaaaaaaaaaaaaaaaaGGATAAGAAGTGGAAAATACCAGTGCATCCTGGCTAGAACTGTCGTAAGGGCTTCTGTAACTGTGTAGAGGTGACAGGACTTCAGGTATGTATGATGAAGGCTCCAGCTGCATATCCCTGTGACTTCCTTAGCAGTGTCTTTCCCATGGGTCAGCTTTGTCTGTATTCATAGGTGGCCACATGCCTGCTTGAGGTCCTGGTGACTCCAGGT'
    """
    fasta = pysam.FastaFile(genome + '.fa')
    seq = fasta.fetch(region=chr_loc)
    return seq


if __name__ == '__main__':
    doctest.testmod()
