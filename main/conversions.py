"""
Convert between DNA sequence formats using using togows.org and crispycrunch
defaults.
"""
import doctest
import re

import requests

# See also CHR_REGEX in validators.py
CHR_REGEX = r'chr([0-9XY]+):([0-9,]+)-([0-9,]+[0-9])'


def chr_loc_to_seq(chr_loc: str, genome: str = 'hg38') -> str:
    """
    >>> chr_loc_to_seq('chr1:12,345-12,500', 'hg38')
    'TCAGACCAGCCGGCTGGAGGGAGGGGCTCAGCAGGTCTGGCTTTGGCCCTGGGAGAGCAGGTGGAAGATCAGGCAGGCCATCGCTGCCACAGAACCCAGTGGATTGGCCTAGGTGGGATCTCTGAGCTCAACAAGCCCTCTCTGGGTGGTAGGTGC'

    >>> chr_loc_to_seq('chr7:5569176-5569415', 'hg38')
    'CTGTCTCAAAAAAAAAAAAAAAAAAAAAAGGATAAGAAGTGGAAAATACCAGTGCATCCTGGCTAGAACTGTCGTAAGGGCTTCTGTAACTGTGTAGAGGTGACAGGACTTCAGGTATGTATGATGAAGGCTCCAGCTGCATATCCCTGTGACTTCCTTAGCAGTGTCTTTCCCATGGGTCAGCTTTGTCTGTATTCATAGGTGGCCACATGCCTGCTTGAGGTCCTGGTGACTCCAGGT'

    >>> chr_loc_to_seq('chr2:38377154-38377424', 'hg38')
    'GTGGACGAGACGTGGTTAACCGCGGCGCTTGGGTCGCTGGTCCGTCGCCGGCGCCACAGCCCCTGGTGCGGTTGCTGCCCTCGCGCTGCCTCGTCCCCCTCCGCCATCTTGTACCGATTTAAAATTAACTCCCCACTGACGTCAAACGCCGCCGCCGCTTTATCAACCTCCCGGGAGCCGGCGGGGTGGCCGGCGGGTCCTAGCGCCGCTCTCCGCCTGTATCTCCTCGCCCTCCGCCTGTATCTCCTCGCCCTCCGCCTGTATCTCCTCG'
    """
    url = 'http://togows.org/api/ucsc/{}/{}.fasta'.format(
        genome, chr_loc)
    response = requests.get(url)
    response.raise_for_status()
    return _reformat_fasta(response.text)


def _reformat_fasta(fasta: str) -> str:
    r"""
    >>> fasta = '''>hg38:chr1:12,345-12,500\n
    ... TCAGACCAGCCGGCTGGAGGGAGGGGCTCAGCAGGTCTGGCTTTGGCCCTGGGAGAGCAG\n
    ... GTGGAAGATCAGGCAGGCCATCGCTGCCACAGAACCCAGTGGATTGGCCTAGGTGGGATC\n
    ... TCTGAGCTCAACAAGCCCTCTCTGGGTGGTAGGTGC'''
    >>> _reformat_fasta(fasta)
    'TCAGACCAGCCGGCTGGAGGGAGGGGCTCAGCAGGTCTGGCTTTGGCCCTGGGAGAGCAGGTGGAAGATCAGGCAGGCCATCGCTGCCACAGAACCCAGTGGATTGGCCTAGGTGGGATCTCTGAGCTCAACAAGCCCTCTCTGGGTGGTAGGTGC'
    """
    return ''.join(fasta.split('\n')[1:])


def gene_to_chr_loc(gene: str, genome='hg38') -> str:
    """
    # TODO (gdingle): we need this to return the length of the amplicon,
    not the whole gene. The amplicon should be less than 500 bp long,
    according to Jason Li.

    Takes the top result from USCS genome browser. See for example:
    https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=ATL2

    >>> gene_to_chr_loc('ATL2')
    'chr2:38294880-38377262'
    >>> gene_to_chr_loc('XXXX') is None
    True
    """
    url = 'https://genome.ucsc.edu/cgi-bin/hgTracks'
    response = requests.get(url, params={
        'db': genome,
        'position': gene,  # "position" is bad name
    })
    response.raise_for_status()
    match = re.search(CHR_REGEX, response.text)
    return match[0] if match else None


if __name__ == '__main__':
    doctest.testmod()
