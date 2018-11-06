"""
Convert between DNA sequence formats using using togows.org and crispycrunch
defaults.
"""
import doctest
import re

import requests
import requests_cache  # type: ignore

# See also CHR_REGEX in validators.py
CHR_REGEX = r'chr([0-9XY]+):([0-9,]+)-([0-9,]+[0-9])'

_cached_session = requests_cache.CachedSession(
    cache_name=__name__ + '_cache',
    # TODO (gdingle): what's the best timeout?
    expire_after=3600 * 24 * 14,
    allowable_methods=('GET', 'POST'),
)

# Avoid too many connections error. See:
# https://stackoverflow.com/questions/23632794/
adapter = requests.adapters.HTTPAdapter(pool_connections=96 * 4, pool_maxsize=96 * 4)
_cached_session.mount('http://', adapter)
_cached_session.mount('https://', adapter)


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
    response = _cached_session.get(url)
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


# TODO (gdingle): switch to togows?
# TODO (gdingle): need to filter better among all the chr results for a gene
# TODO (gdingle): how does this deal with reverse strand?
def gene_to_chr_loc(gene: str, genome: str ='hg38') -> str:
    """
    Takes the top result from USCS genome browser. See for example:
    https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=ATL2

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
    response = _cached_session.get(url, params={
        'db': genome,
        'position': gene,  # "position" is bad name
    })
    response.raise_for_status()
    match = re.search(CHR_REGEX, response.text)
    if not match:
        raise ValueError('No chr location for {}'.format(gene))
    return match[0]


def enst_to_gene(enst: str, genome: str = 'hg38') -> str:
    """
    Gets NCBI gene associated with with an ENST transcript ID.

    >>> enst_to_gene('ENST00000398844')
    'SEC24A'
    >>> enst_to_gene('ENST00000617316')
    Traceback (most recent call last):
    ...
    requests.exceptions.HTTPError: 404 Client Error: Not Found for url: http://togows.org/search/ncbi-gene/ENST00000617316/1,50.json
    """
    url = 'http://togows.org/search/ncbi-gene/{}/1,50.json'.format(enst)
    response = _cached_session.get(url)
    response.raise_for_status()
    res = response.json()
    assert len(res) == 1
    ncbi_id = res[0]

    url = 'http://togows.org/entry/ncbi-gene/{}/official_symbol.json'.format(ncbi_id)
    response = _cached_session.get(url)
    response.raise_for_status()
    res = response.json()
    assert len(res) == 1
    return res[0]


def enst_to_gene_or_unknown(enst: str, genome: str = 'hg38') ->str:
    """
    Suppresses not foundn exceptions.
    >>> enst_to_gene_or_unknown('ENST00000617316')
    'UNKNOWN'
    """
    try:
        return enst_to_gene(enst, genome)
    except IOError:
        return 'UNKNOWN'


def chr_loc_to_gene(chr_loc: str, genome: str = 'hg38', straddle: bool = True) -> str:
    """
    Gets the gene symbols associated with a chromosome range.

    By default, any part of the range may overlap the gene.

    Returns '' if no match found.

    See for example:
        http://togows.org/api/ucsc/hg38/refGene/exclusive/chr4:1,350,000-1,400,000.json
        http://togows.org/api/ucsc/hg38/refGene/exclusive/chr2:38294880-38377262.json

    >>> chr_loc_to_gene('chr2:136114349-136116243', straddle=True)
    'CXCR4'

    >>> chr_loc_to_gene('chr2:38294880-38377262', straddle=True)
    'ATL2'

    >>> chr_loc_to_gene('chr2:136114349-136116243', straddle=False)
    'CXCR4'

    >>> chr_loc_to_gene('chr2:38294880-38377262', straddle=False)
    ''
    """
    url = 'http://togows.org/api/ucsc/{}/refGene/{}/{}.json'.format(
        genome, 'inclusive' if straddle else 'exclusive', chr_loc)
    response = _cached_session.get(url)
    if response.status_code == 404:
        return ''
    else:
        response.raise_for_status()
    res = response.json()

    matches = set(r['name2'] for r in res)
    # TODO (gdingle): how to handle multiple gene matches?
    assert len(matches) == 1, matches
    return matches.pop()

# TODO (gdingle): share code with protospacex
# For gene_to_enst, see protospacex


if __name__ == '__main__':
    doctest.testmod()

    # print(chr_loc_to_seq('chr7:2251813-2251848'))
    # print(chr_loc_to_seq('chr7:2251848-2251813'))
