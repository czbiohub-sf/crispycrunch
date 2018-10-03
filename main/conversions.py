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
    >>> chr_loc_to_seq('chr2:136114349-136116243', 'hg38')
    'TAAGTTTAACATGTACTTTTATTAACAACTTAATACAAGACTGTACACTGTAGGTGCTGAAATCAACCCACTCCTGAAAACTGAAAAACCAGCATTTCTATACCACTTTGGGCTTTGGTTATAAGTGCCATCTTCTACAGCAAAATCACGTCTTAAGAACAGGAAAAACGTTCCACGGGAATGGAGAGATTATCTATGCATAAACAGCTGGGGATCATTTCTAGCTTTACGTGATTCACTACACGCTCTGGAATGTTCAGTTCCCTTTTCTACAGTCCTACCACGAGACATACAGCAACTAAGAACTTGGCCACAGGTCCTGCCTAGACACACATCAATATGAAACAAAAAAAATTTATATAAATAAGTCAATTAAACTTCACAAAAACTAAAGAAACACAAGACAAAAATCCAACAAGCAATAAAAACTGTACAATATTGGTCAGTCTTTTATATCTGAAAAATGTGTAACTTAAAAAAAAGTTATTTATCGTATAAAAAAAAGTCTTTTACATCTGTGTTAGCTGGAGTGAAAACTTGAAGACTCAGACTCAGTGGAAACAGATGAATGTCCACCTCGCTTTCCTTTGGAGAGGATCTTGAGGCTGGACCCTCTGCTCACAGAGGTGAGTGCGTGCTGGGCAGAGGTTTTAAATTTGGCTCCAAGGAAAGCATAGAGGATGGGGTTCAGACAACAGTGGAAGAAAGCTAGGGCCTCGGTGATGGAAATCCACTTGTGCACAGTGTTCTCAAACTCACACCCTTGCTTGATGATTTCCAGGAGGATGAAGGAGTCGATGCTGATCCCAATGTAGTAAGGCAGCCAACAGGCGAAGAAAGCCAGGATGAGGATGACTGTGGTCTTGAGGGCCTTGCGCTTCTGGTGGCCCTTGGAGTGTGACAGCTTGGAGATGATAATGCAATAGCAGGACAGGATGACAATACCAGGCAGGATAAGGCCAACCATGATGTGCTGAAACTGGAACACAACCACCCACAAGTCATTGGGGTAGAAGCGGTCACAGATATATCTGTCATCTGCCTCACTGACGTTGGCAAAGATGAAGTCGGGAATAGTCAGCAGGAGGGCAGGGATCCAGACGCCAACATAGACCACCTTTTCAGCCAACAGCTTCCTTGGCCTCTGACTGTTGGTGGCGTGGACGATGGCCAGGTAGCGGTCCAGACTGATGAAGGCCAGGATGAGGACACTGCTGTAGAGGTTGACTGTGTAGATGACATGGACTGCCTTGCATAGGAAGTTCCCAAAGTACCAGTTTGCCACGGCATCAACTGCCCAGAAGGGAAGCGTGATGACAAAGAGGAGGTCGGCCACTGACAGGTGCAGCCTGTACTTGTCCGTCATGCTTCTCAGTTTCTTCTGGTAACCCATGACCAGGATGACCAATCCATTGCCCACAATGCCAGTTAAGAAGATGATGGAGTAGATGGTGGGCAGGAAGATTTTATTGAAATTAGCATTTTCTTCACGGAAACAGGGTTCCTTCATGGAGTCATAGTCCCCTGAGCCCATTTCCTCGGTGTAGTTATCTGAAGTGTATATCTGCAAAAGAGGCAAAGGAATGGACATTCACTTCCAATTCAGCAAGCATTAACCCAGTTAAAAAAAATTTTTAAAGCAATTTAAAAAACCAATTCAGGCTTGCTTTCTTCAGGAAATTCTGAAGTAGTGGGCTAAGGGCACAAGAGAATTAATGTAGAATCCTACAACTCTCCTCCCCATCTTTTCCCATAGTGACTTCATTATATCCTTCTTTGGTAGAACCAATTACAAAATTCTTTGTTTAGAACAAAAGGGCACTGAGACGCTGAGGGTTTCAAAGTCACATCTTGGCTAACTCCTCTGCCCCGCCCACTAGAGGGAAGAAAAAAAA'

    >>> chr_loc_to_seq('chr2:136114380-136114399', 'hg38')
    'AATACAAGACTGTACACTGT'
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


def gene_to_chr_loc(gene: str, genome: str ='hg38') -> str:
    """
    Takes the top result from USCS genome browser. See for example:
    https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=ATL2

    # TODO (gdingle): how does this deal with reverse strand?

    >>> gene_to_chr_loc('ATL2')
    'chr2:38294880-38377262'
    >>> gene_to_chr_loc('XXXX')
    Traceback (most recent call last):
    ...
    ValueError: No chr location for XXXX
    >>> gene_to_chr_loc('cxcr4')
    'chr2:136114349-136116243'

    Also handles ENST transcripts.
    # TODO (gdingle): how reliable is this?

    >>> gene_to_chr_loc('ENST00000398844')
    'chr5:134648789-134727823'
    """
    url = 'https://genome.ucsc.edu/cgi-bin/hgTracks'
    response = requests.get(url, params={
        'db': genome,
        'position': gene,  # "position" is bad name
    })
    response.raise_for_status()
    match = re.search(CHR_REGEX, response.text)
    if not match:
        raise ValueError('No chr location for {}'.format(gene))
    return match[0]


if __name__ == '__main__':
    doctest.testmod()
