"""
Convert between DNA sequence formats using using togows.org and crispycrunch
defaults.
"""
import doctest
import logging
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

logger = logging.getLogger(__name__)


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

    >>> chr_loc_to_seq('chr3:128067063-128067085', 'hg38')
    'CCCTGGTGTCCAACCTTTATGTC'

    >>> chr_loc_to_seq('chr3:128067063-128067085:-', 'hg38')
    'GACATAAAGGTTGGACACCAGGG'

    Mouse.
    >>> chr_loc_to_seq('chr7:28179081-28179112:+', 'mm10')
    'CTGAGAAGCCAAAAGTGGTTACAACTCGACCC'
    """
    if chr_loc.endswith((':+', ':-')):
        chr_loc, strand = chr_loc[:-2], chr_loc[-1]
        if strand == '-':
            matches = re.search(CHR_REGEX, chr_loc)
            # flip for togows
            if matches:
                chr_loc = f'chr{matches[1]}:{matches[3]}-{matches[2]}'

    url = 'http://togows.org/api/ucsc/{}/{}.fasta'.format(
        genome, chr_loc)
    response = _cached_session.get(url)
    response.raise_for_status()
    return _reformat_fasta(response.text).upper()


def seq_to_chr_loc(seq: str, genome: str = 'hg38') -> str:
    """
    Looks up chr loc by exact sequence match.

    >>> seq_to_chr_loc('AAGGTGAAGAACTGAAGTTCAGCGCTGTCA')
    'chr19:39834523-39834552:-'

    Matched with Crispor: http://crispor.tefor.net/crispor.py?batchId=PpKGGNjVF5Bjt5u3G88P
    >>> seq_to_chr_loc('CACCTCGAGCTCTCGCACCAGGG')
    'chr11:134177387-134177409:+'

    No match.
    >>> seq_to_chr_loc('CACCTCGAGCTCTCGCACCAGGC')
    Traceback (most recent call last):
    ...
    ValueError: No match for sequence CACCTCGAGCTCTCGCACCAGGC in genome hg38

    Mouse.
    >>> seq_to_chr_loc('CTGAGAAGCCAAAAGTGGTTACAACTCGACCC', 'mm10')
    'chr7:28179081-28179112:+'
    """

    url = 'http://gggenome.dbcls.jp/{}/0/{}.json'.format(
        genome, seq)
    response = _cached_session.get(url)
    response.raise_for_status()
    data = response.json()

    if data['error'] != 'none':
        raise RuntimeError('gggenome error: "{}"'.format(data['error']))

    results = data['results']
    if len(results) == 0:
        raise ValueError('No match for sequence {} in genome {}'.format(
            seq, genome))
    if len(results) > 1:
        raise ValueError('More than one match for sequence {} in genome {}'.format(
            seq, genome))
    return '{}:{}-{}:{}'.format(
        results[0]['name'],
        results[0]['position'],
        results[0]['position_end'],
        results[0]['strand'],
    )


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
    'chr2:38295901-38377273'
    >>> gene_to_chr_loc('XXXX')
    Traceback (most recent call last):
    ...
    ValueError: No chr location for XXXX
    >>> gene_to_chr_loc('cxcr4')
    'chr2:136114349-136116243'

    Mouse
    >>> gene_to_chr_loc('CNTNAP1', 'mm10')
    'chr11:101176117-101190720'
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


def enst_to_gene(enst: str, genome: str = 'hg38', timeout=4.0) -> str:
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
    response = _cached_session.get(url, timeout=timeout / 2)
    response.raise_for_status()
    res = response.json()
    assert len(res) == 1
    ncbi_id = res[0]

    url = 'http://togows.org/entry/ncbi-gene/{}/official_symbol.json'.format(ncbi_id)
    response = _cached_session.get(url, timeout=timeout / 2)
    response.raise_for_status()
    res = response.json()
    assert len(res) == 1
    return res[0]


def enst_to_gene_or_unknown(enst: str, genome: str = 'hg38') ->str:
    """
    Suppresses not found exceptions.
    >>> enst_to_gene_or_unknown('ENST00000617316')
    'UNKNOWN'
    """
    try:
        return enst_to_gene(enst, genome, timeout=4)
    except requests.exceptions.Timeout:
        return 'TIMEOUT'
    except IOError as e:
        logger.warning(e)
        return 'UNKNOWN'


def chr_loc_to_gene(
    chr_loc: str,
    genome: str = 'hg38',
    straddle: bool = True,
    strand: str = None,
) -> str:
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
    >>> chr_loc_to_gene('chr3:128067063-128067085')
    'MULTIPLE'
    >>> chr_loc_to_gene('chr3:128067063-128067085', strand='+')
    'SEC61A1'
    >>> chr_loc_to_gene('chr3:128067063-128067085', strand='-')
    'RUVBL1'

    Mouse.
    >>> chr_loc_to_gene('chr11:101176117-101190720', genome='mm10')
    'Cntnap1'
    """
    if chr_loc.endswith((':+', ':-')):
        chr_loc, strand = chr_loc[:-2], chr_loc[-1]

    url = 'http://togows.org/api/ucsc/{}/refGene/{}/{}.json'.format(
        genome, 'inclusive' if straddle else 'exclusive', chr_loc)
    response = _cached_session.get(url)
    if response.status_code == 404:
        return ''
    else:
        response.raise_for_status()
    res = response.json()

    strands = [strand] if strand else ['+', '-']
    matches = set(r['name2'] for r in res if r['strand'] in strands)
    if not len(matches):
        return 'UNKNOWN'
    if len(matches) > 1:
        return 'MULTIPLE'
    return matches.pop()

# TODO (gdingle): share code with protospacex
# For gene_to_enst, see protospacex


if __name__ == '__main__':
    doctest.testmod()

    # print(enst_to_gene_or_unknown('ENST00000278840'))
    # print(chr_loc_to_seq('chr7:2251813-2251848'))
    # print(chr_loc_to_seq('chr7:2251848-2251813'))
