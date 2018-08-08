import doctest
import os

import pysam

# TODO (gdingle): download hg38.fa and mm10.fa from public repo in docker build
# TODO (gdingle): move this module to crispycrunch?


def get_reference_amplicon(chr_loc, genome='hg38'):
    """
    >>> get_reference_amplicon('chr7:5569176-5569415')
    'ctgtctcaaaaaaaaaaaaaaaaaaaaaaGGATAAGAAGTGGAAAATACCAGTGCATCCTGGCTAGAACTGTCGTAAGGGCTTCTGTAACTGTGTAGAGGTGACAGGACTTCAGGTATGTATGATGAAGGCTCCAGCTGCATATCCCTGTGACTTCCTTAGCAGTGTCTTTCCCATGGGTCAGCTTTGTCTGTATTCATAGGTGGCCACATGCCTGCTTGAGGTCCTGGTGACTCCAGGT'
    """
    path = os.path.join(os.path.dirname(__file__), genome + '.fa')
    fasta = pysam.FastaFile(path)
    seq = fasta.fetch(region=chr_loc)
    return seq


if __name__ == '__main__':
    doctest.testmod()
