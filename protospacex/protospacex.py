"""
This module was extracted and adapted from https://github.com/czbiohub/protospacex,
which was designed for https://czi.quip.com/YbAhAbOV4aXi/.

Protospacex: automated guide design for Cas9 knock-in experiments.

The code here returns different regions of interest for HDR from a ENST transcript.
"""
import logging
import time

from functools import lru_cache
# lru_cache appears to reduce page time by 10s for 96 well plate!

import requests
import requests_cache  # type: ignore

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqFeature import FeatureLocation  # type: ignore
from Bio.SeqFeature import SeqFeature  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.WARN)


_cached_session = requests_cache.CachedSession(
    cache_name=__name__ + '_cache',
    # TODO (gdingle): what's the best timeout?
    expire_after=3600 * 24 * 14,
    allowable_methods=('GET', 'POST'),
)

# TODO (gdingle): refactor with conversions.py
# Avoid too many connections error. See:
# https://stackoverflow.com/questions/23632794/
adapter = requests.adapters.HTTPAdapter(pool_connections=96 * 4, pool_maxsize=96 * 4)
_cached_session.mount('http://', adapter)
_cached_session.mount('https://', adapter)


@lru_cache(maxsize=1024)
def fetch_ensembl_transcript(
        ensembl_transcript_id: str,
        is_retry: bool = False) -> SeqRecord:
    """Fetch the requested Ensembl transcript.

    Get the requested Ensembl transcript, together with exon and
    coding region (CDS) boundaries.

    Parameters
    ----------
    ensembl_transcript_id : str
      the ensembl transcript id, of the form ENST...

    Returns
    -------
    `Bio.SeqRecord`

      The requested transcript sequence, in 5' -> 3' order, together
      with exon and CDS features. The coordinates of exons and CDS
      features are relative to the sequence fragment.

    >>> fetch_ensembl_transcript('ENST00000398844').description
    'chromosome:GRCh38:5:134648789:134727823:1'

    >>> fetch_ensembl_transcript('ATL3').description
    'Reverse complement of chromosome:GRCh38:11:63624087:63671612:-1'

    >>> fs = fetch_ensembl_transcript('ENST00000398844').features
    >>> len([f for f in fs if f.type == 'exon'])
    23
    """
    base_url = "http://rest.ensembl.org"

    if not ensembl_transcript_id.startswith('ENS'):
        # could be a gene symbol
        ensembl_transcript_id = _gene_to_enst(ensembl_transcript_id)

    # First, fetch the transcript sequence
    url = base_url + f"/sequence/id/{ensembl_transcript_id}"

    log.debug(f"Querying Ensembl for sequence of {ensembl_transcript_id}")
    response = _cached_session.get(url, params={"type": "genomic",
                                                "content-type": "application/json"})
    log.debug('Request cached: {}'.format(getattr(response, 'from_cache', False)))
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError:
        log.error("Ensembl sequence REST query returned error "
                  "{}".format(response.text))
        raise ValueError(response.text)

    response_data = response.json()

    try:
        description = response_data['desc'].split(':')
        species = description[1]
        chromosome_number = description[2]  # may be X
        sequence_left = int(description[3])
        sequence_right = int(description[4])
        transcript_strand = int(description[5])

        if sequence_left > sequence_right:
            raise ValueError(f"Expected left sequence boundary {sequence_left} "
                             f"<= right sequence boundary {sequence_right}: did "
                             "the format of the Ensembl REST response change?")

        sequence_id = response_data['id']

        seq_str = response_data['seq']

        log.debug(f"Retrieved sequence {response_data['desc']} of length "
                  f"{sequence_right - sequence_left} for species {species} on "
                  f"strand {transcript_strand}")
    except (KeyError, ValueError) as e:
        log.error(e)
        log.error('Error parsing sequence metadata from Ensembl REST response - '
                  'did the format of the response change?')
        raise ValueError(e)

    if transcript_strand == -1:
        # If the transcript strand is -1, the sequence returned by
        # Ensembl is on the strand opposite the reference strand,
        # which is the strand of the Ensembl coordinates for
        # exons/coding regions. In this case, we initially store the
        # reverse complement of the sequence, and after fetching the
        # exon/coding regions, we'll return the reverse complement of
        # the `Bio.SeqRecord` object, which will properly re-index the
        # exon/coding regions.
        seq = Seq(seq_str, IUPACUnambiguousDNA()).reverse_complement()
    else:
        seq = Seq(seq_str, IUPACUnambiguousDNA())

    record = SeqRecord(seq, id=sequence_id,
                       description=":".join(description))

    url = base_url + f"/overlap/id/{ensembl_transcript_id}"

    log.debug(f"Querying Ensembl for overlaps of {ensembl_transcript_id}")
    response = _cached_session.get(url, params={"feature": ["cds", "exon"],
                                                "content-type": "application/json"})
    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError:
        log.error("Ensembl sequence REST query returned error "
                  "{}".format(response.text))
        # You have exceeded the limit of 15 requests per second; please reduce your concurrent connections
        if not is_retry and \
                response.text.strip().startswith('You have exceeded the limit of'):
            log.warn('Waiting 4 seconds then retrying fetch_ensembl_transcript')
            time.sleep(4)
            return fetch_ensembl_transcript(ensembl_transcript_id, is_retry=True)

        raise ValueError(response.text)

    response_data = response.json()

    try:
        # Handle the unlikely event of a single piece of information
        # overlapping a lonely transcript
        if not hasattr(response_data, '__iter__'):
            response_data = [response_data]

        for response_datum in response_data:
            if response_datum['Parent'] != ensembl_transcript_id:
                continue

            if response_datum['assembly_name'] != species:
                continue

            # We store feature locations 0-indexed from the left-most
            # sequence boundary
            record.features.append(SeqFeature(
                location=FeatureLocation(
                    int(response_datum['start']) - sequence_left,
                    int(response_datum['end']) - sequence_left + 1,
                    strand=int(response_datum['strand'])),
                type=response_datum['feature_type']))
        num_exon_boundaries = len([f for f in record.features
                                   if f.type == 'exon'])

        num_cds_boundaries = len([f for f in record.features
                                  if f.type == 'cds'])

        log.debug(f"Retrieved {num_exon_boundaries} exons and "
                  f"{num_cds_boundaries} coding regions for transcript "
                  f"{ensembl_transcript_id}")
    except (KeyError, ValueError) as e:
        log.error(e)
        log.error('Error parsing overlap metadata from Ensembl REST response - '
                  'did the format of the response change?')
        raise ValueError(e)

    if transcript_strand == -1:
        # By default `reverse_complement` doesn't preserve
        # description, so force it...
        record = record.reverse_complement(description=True)

        # ...but update the description to make clear the sequence
        # we're storing is the reverse complement of the sequence
        # described by the metadata in the description
        record.description = "Reverse complement of " + record.description

    record.annotations['reference_species'] = species
    record.annotations['reference_chromosome_number'] = chromosome_number
    record.annotations['reference_left_index'] = sequence_left
    record.annotations['reference_right_index'] = sequence_right
    record.annotations['transcript_strand'] = transcript_strand

    # Finally, sort features by their start locations
    record.features.sort(key=lambda f: f.location.start)

    return record


def get_cds_seq(
        ensembl_transcript_id: str,
        cds_index: int = 0,
        length: int = 36) -> str:
    """
    Return base pair sequence surrounding codon.

    The codon may the start or stop codon or other depending on cds_index.

    See https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000113615;r=5:134648789-134727823;t=ENST00000398844
    See https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS43363

    For param docs, see get_cds_chr_loc

    Start codon.

    >>> get_cds_seq('ENST00000398844')
    'CCTTCCAACCCAGTCATCATGTCCCAGCCGGGAATA'

    >>> get_cds_seq('ENST00000411809')
    'CGTCCCGCGCGGGGCACGATGTTGAACATGTGGAAG'

    Stop codon.

    >>> get_cds_seq('ENST00000398844', -1)
    'CAACAAGTGAATAAATGAATGAATGAAGAAATTTGA'

    >>> get_cds_seq('ENST00000411809', -1)
    'GCCAATTTTAGCAAATAAGAGATTGTAAAAGAAGCA'

    Length 30.
    >>> len(get_cds_seq('ENST00000221801', -1, length=30))
    30

    No length change. Intron/exon junction.

    >>> get_cds_seq('ENST00000221801', -1, length=-1)
    'gccaCCCCCCAAGGTGAAGAACTGA'

    Too short.
    >>> get_cds_seq('ENST00000262033', length=-1)
    ''
    >>> get_cds_seq('ENST00000295682', -1, length=72)
    'GCCAAGGTCACAGGCAAGAGCAAGAAGAGAAACTGACCCTGAATGTTCAATAAAGTTGATTCTTTGTAGCTC'
    """
    return _get_cds_seq_and_codon_at(
        ensembl_transcript_id,
        cds_index,
        length)[0]


def get_cds_codon_at(
        ensembl_transcript_id: str,
        cds_index: int = 0,
        length: int = 36) -> tuple:
    return _get_cds_seq_and_codon_at(
        ensembl_transcript_id,
        cds_index,
        length)[1]


def _get_cds_seq_and_codon_at(
        ensembl_transcript_id: str,
        cds_index: int = 0,
        length: int = 36) -> tuple:
    record = fetch_ensembl_transcript(ensembl_transcript_id)
    cds = [f for f in record.features if f.type == 'cds']
    assert len(cds)

    cds_seq = _get_cds_seq(record, cds, cds_index)

    # enforce length
    _validate_length(length)
    start, end, codon_at = _get_start_end(
        cds[cds_index].location,
        length,
        cds_index
    )

    # Edge case
    if length == -1 and end - start < 3:
        return '', -1

    if end > len(record.seq) or start < 0:
        cds_seq = _get_seq_from_surrounding(record, start, end)
    else:
        cds_seq = record.seq[start:end]
        # TODO (gdingle): do for above also somehow... how to deal with out of range?
        cds_seq = _lowercase_exon_boundaries(cds_seq, record, start)

    if cds_index == 0:
        target_codons = ['ATG']
    elif cds_index == -1:
        target_codons = ['TAG', 'TGA', 'TAA']
    assert cds_seq[codon_at:codon_at + 3].upper() in target_codons, (
        str(cds_seq), codon_at)

    if length != -1:
        assert len(cds_seq) == length, (len(cds_seq), start, end, length)
        assert len(cds_seq) % 3 == 0, 'must be codon aligned'

    if len(cds_seq) < 3:
        log.warning('CCDS of {} is too short "{}", so returning empty string'.format(
            ensembl_transcript_id, cds_seq))
        return '', 0

    return str(cds_seq), codon_at


def _lowercase_exon_boundaries(
        cds_seq: str,
        record: SeqRecord,
        start: int = 0,
        bound_len: int = 4) -> str:
    """
    >>> seq = Seq('AAGGTGAAGAACTGAAGTTCAGCGCTGTCA', IUPACUnambiguousDNA())
    >>> features = [SeqFeature(location=FeatureLocation(15, 15), type='exon')]
    >>> record = SeqRecord(seq, features=features)
    >>> _lowercase_exon_boundaries(str(seq), record)
    'AAGGTGAAGAActgaagttCAGCGCTGTCA'
    >>> features[0].location = FeatureLocation(0, 0)
    >>> _lowercase_exon_boundaries(str(seq), record)
    'aaggTGAAGAACTGAAGTTCAGCGCTGTCA'
    >>> features[0].location = FeatureLocation(0, 32)
    >>> _lowercase_exon_boundaries(str(seq), record)
    'aaggTGAAGAACTGAAGTTCAGCGCTGTca'
    """
    exons = [f.location for f in record.features if f.type == 'exon']
    for exon in exons:
        for exon_boundary in [exon.start - start, exon.end - start]:
            left = max(0, exon_boundary - bound_len)
            right = min(len(cds_seq), exon_boundary + bound_len)
            cds_seq = cds_seq[:left] \
                + cds_seq[left:right].lower() \
                + cds_seq[right:]
    return cds_seq


def get_cds_chr_loc(
        ensembl_transcript_id: str,
        cds_index: int = 0,
        length: int = 36) -> str:
    """
    Return genome region surrounding codon.

    The codon may the start or stop codon or other depending on cds_index.

    length specifies the region around the codon
        If start codon, within length *after* codon.
        If stop codon, within length *centered* on codon.

    If stop codon, the start is adjusted, if start codon, the end is adjusted.

    length must be divisible by 2 and 3 to maintain codon frame and symmetry.

    length -1 will return the cds location unchanged.

    See https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000113615;r=5:134648789-134727823;t=ENST00000398844
    See https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS43363

    >>> get_cds_chr_loc('ENST00000398844', length=990)
    'chr5:134648582-134649571:+'

    >>> get_cds_chr_loc('ENST00000411809', length=990)
    'chr5:157858476-157859465:-'

    Length.

    >>> get_cds_chr_loc('ENST00000221801', length=990)
    'chr19:39845806-39846795:-'

    >>> get_cds_chr_loc('ENST00000398844', -1, 90)
    'chr5:134725050-134725139:+'

    Get last codon.

    >>> get_cds_chr_loc('ENST00000398844', -1, 990)
    'chr5:134724600-134725589:+'

    >>> get_cds_chr_loc('ENST00000411809', -1, 990)
    'chr5:157787151-157788140:-'

    Length, last codon.

    >>> get_cds_chr_loc('ENST00000221801', -1, length=30)
    'chr19:39834523-39834552:-'

    Gggenome says:
    http://gggenome.dbcls.jp/hg38/AAGGTGAAGAACTGAAGTTCAGCGCTGTCA
    chr19:39834523-39834552

    Length.

    >>> get_cds_chr_loc('ENST00000398844')
    'chr5:134649059-134649094:+'

    >>> get_cds_chr_loc('ENST00000398844', -1)
    'chr5:134725077-134725112:+'

    No length change.

    >>> get_cds_chr_loc('ENST00000221801', -1, length=-1)
    'chr19:39834538-39834562:-'

    Weird case.
    >>> get_cds_chr_loc('ENST00000638706', -1)
    'chrCHR_HG30_PATCH:179729500-179729535:+'

    Across splice boundaries.
    >>> get_cds_chr_loc('ENST00000262033', length=36)
    'chr12:56687982-56688017:-'
    >>> get_cds_chr_loc('ENST00000054666', length=36)
    'chr1:7771366-7771401:+'
    """
    record = fetch_ensembl_transcript(ensembl_transcript_id)
    cds = [f for f in record.features if f.type == 'cds']
    assert len(cds)
    cds_location = cds[cds_index].location

    species = record.annotations['reference_species']
    chromosome_number = record.annotations['reference_chromosome_number']
    sequence_left = record.annotations['reference_left_index']
    sequence_right = record.annotations['reference_right_index']
    transcript_strand = record.annotations['transcript_strand']

    seq_len = sequence_right - sequence_left + 1
    assert seq_len == len(record.seq)

    # Transcript of reverse strand is translated above.
    assert cds_location.strand == 1
    assert species == 'GRCh38'

    # enforce length
    _validate_length(length)
    start, end, codon_at = _get_start_end(
        cds_location,
        length,
        cds_index
    )

    if transcript_strand == 1:
        start = sequence_left + start
        end = sequence_left + end
    else:
        start, end = sequence_right - end, sequence_right - start
        # Not sure why, but this makes locations agree with gggenome and crispor
        start, end = start + 1, end + 1

    return 'chr{}:{}-{}:{}'.format(
        chromosome_number,
        start,
        end - 1,  # change to inclusive range
        '+' if transcript_strand == 1 else '-',
    )


def get_ultramer_seq(
        ensembl_transcript_id: str,
        cds_index: int = 0,
        length: int = 110) -> tuple:
    """
    Function to get the precise sequence centered around the target codon
    insert position needed for IDT ultramer ordering (donor DNA template).

    Length default of 110 is determined by need for max donor length of 200 bp.
    See https://www.idtdna.com/pages/education/decoded/article/crispr-cas9-mediated-hdr-tips-for-successful-experimental-design .

    >>> get_ultramer_seq('ENST00000398844')[0]
    'TCTTCTTGTGCGCTGTTGTCGACCCCGACCAGCCCCTTCCAACCCAGTCATCATGTCCCAGCCGGGAATACCGGCCTCCGGCGGCGCCCCAGCCAGCCTCCAGGCCCAGA'

    >>> get_ultramer_seq('ENST00000411809')[0]
    'GGTCCGTGGGGAGCAGGAGAGGGAGGCGGCGGACCGTCCCGCGCGGGGCACGATGTTGAACATGTGGAAGGTGCGCGAGCTGGTGGACAAAGCGTGAGTATCGGGGGGCA'

    Stop codon.

    >>> get_ultramer_seq('ENST00000398844', -1)[0]
    'ATCTGCATTATCATATTATGAATTCCTGTTGCATATACAGCAACAAGTGAATAAATGAATGAATGAAGAAATTTGACTTATTTTTAAGGAATGTCACGATAGTGCAGAAT'

    >>> get_ultramer_seq('ENST00000411809', -1)[0]
    'TGGAACTGTGCAACCCAAGCAAGATGCCTTTGCAAATTTCGCCAATTTTAGCAAATAAGAGATTGTAAAAGAAGCAGATTGAATGAAGAATTTTTAGCTGTGCAGATAGG'

    >>> get_ultramer_seq('ENST00000221801', -1)[0]
    'GACCCTCCTTCATCACCTATCTTCCTCTCACAGGCCACCCCCCAAGGTGAAGAACTGAAGTTCAGCGCTGTCAGGATTGCGAGAGATGTGTGTTGATACTGTTGCACGTG'

    Not enough in the transcript for the ultramer.
    >>> get_ultramer_seq('ENST00000258648')[0]
    'CACCTGCGTCAGCTCGCTCTGCGCGTGCGCCGGTGGCGGGACTCTGGGGAAAATGGCTGCGTCTTCGAGTGGTGAGAAGGAGAAGGAGCGGCTGGGAGGCGGTTTGGGAG'

    >>> get_ultramer_seq('ENST00000267113')[0]
    'GGCAACCTCCAGCCAGTCCCTGGGTCGGGCGGATCCTCCCAGAGGTGGCACAATGGAGCGATCTCCAGGAGAGGGCCCCAGCCCCAGCCCCATGGACCAGCCCTCTGCTC'
    """

    record = fetch_ensembl_transcript(ensembl_transcript_id)
    cds = [f for f in record.features if f.type == 'cds']
    assert len(cds)

    location = cds[cds_index].location
    # see also _get_start_end
    start = location.start
    end = location.end
    if cds_index == -1:  # stop codon
        codon_at = end - 3
        start = codon_at - length // 2
        end = codon_at + length // 2
    elif cds_index == 0:  # start codon
        codon_at = start
        insert_at = codon_at + 3
        end = insert_at + length // 2
        start = insert_at - length // 2
    assert end - start == length

    ult_seq = record.seq[start:end]

    if start < 0 or end > len(record.seq):
        log.warning('Transcript {} is not large enough: {}. Fetching sequence.'.format(
            ensembl_transcript_id, len(ult_seq)))
        ult_seq = _get_seq_from_surrounding(record, start, end)
        assert len(ult_seq) == length, (len(ult_seq), length)

    codon_at = codon_at - start  # make relative to
    if cds_index == 0:
        assert ult_seq[codon_at:codon_at + 3] == 'ATG', (start, end, codon_at, ult_seq)
    elif cds_index == -1:
        assert ult_seq[codon_at:codon_at + 3] in ['TAG', 'TGA', 'TAA']

    return str(ult_seq), codon_at


def _get_seq_from_surrounding(record, start: int, end: int) -> str:
    species = record.annotations['reference_species']
    chr_num = record.annotations['reference_chromosome_number']
    transcript_strand = record.annotations['transcript_strand']
    if transcript_strand == -1:
        seq_end = record.annotations['reference_right_index']
        seq = _fetch_seq(
            species,
            chr_num,
            seq_end - end + 1,
            # Note: ensembl is inclusive range and biopython is exclusive
            seq_end - start
        )
        seq = Seq(seq, IUPACUnambiguousDNA()).reverse_complement()
    else:
        seq_start = record.annotations['reference_left_index']
        seq = _fetch_seq(
            species,
            chr_num,
            seq_start + start,
            # Note: ensembl is inclusive range and biopython is exclusive
            seq_start + end - 1
        )
    return seq


def _get_start_end(
        location: FeatureLocation,
        length: int,
        cds_index: int) -> tuple:
    start = location.start
    end = location.end

    # No-op
    if length == -1:
        return start, end, 0 if cds_index == 0 else end - start - 3

    assert length > 0

    if cds_index == -1:  # stop_codon
        start = end - length // 2
        codon_at = end - 3
        end = end + length // 2
    elif cds_index == 0:  # start_codon2
        codon_at = start
        start = start - length // 2
        end = codon_at + length // 2
    # TODO (gdingle): need both start_codon and start_codon2?
    # elif cds_index == 0:  # legacy start_codon
    #     end = start + length
    #     codon_at = start

    assert end - start == length
    return start, end, codon_at - start


def _validate_length(length: int) -> None:
    if length == -1:
        return None

    def divisible(i: int) -> bool:
        return i % 3 == 0 and i % 2 == 0

    if not divisible(length):
        raise ValueError(
            f'length {length} must be divisible by 3 and 2 to ensure codon frame and symmetry')


# TODO (gdingle): share code with conversions.py
def _gene_to_enst(gene: str, genome: str = 'hg38') -> str:
    """
    >>> _gene_to_enst('ATL3')
    'ENST00000398868'

    >>> _gene_to_enst('ATL3', 'wooky')
    Traceback (most recent call last):
    ...
    ValueError: Unsupported genome: "wooky"

    >>> _gene_to_enst('ATL3', 'mm10')
    'ENSMUST00000025668'

    >>> _gene_to_enst('ATL9')
    Traceback (most recent call last):
    ValueError: No gene of name "ATL9" found
    """
    if genome.startswith('hg'):
        species = 'human'
    elif genome.startswith('mm'):
        species = 'mouse'
    else:
        raise ValueError('Unsupported genome: "{}"'.format(genome))

    url = 'http://rest.ensembl.org/lookup/symbol/{}/{}?expand=1;content-type=application/json'.format(
        species, gene)
    response = _cached_session.get(url)
    if response.status_code == 400:
        raise ValueError('No gene of name "{}" found'.format(gene))
    else:
        response.raise_for_status()

    res = response.json()
    transcripts = [t for t in res['Transcript'] if t['is_canonical']]
    assert len(transcripts) == 1
    transcript = transcripts[0]['id']
    assert transcript.startswith('ENS')
    return transcript


def _fetch_seq(species: str, chr_num: str, start: int, end: int) -> str:
    """
    >>> _fetch_seq('human', 'X', 1000000, 1000100)
    'CTGTAGAAACATTAGCCTGGCTAACAAGGTGAAACCCCATCTCTACTAACAATACAAAATATTGGTTGGGCGTGGTGGCGGGTGCTTGTAATCCCAGCTAC'
    """
    # TODO (gdingle): refactor with conversions.chr_loc_to_seq?
    if species.startswith(('GRCh', 'human')):
        species = 'human'
    elif species.startswith(('GRCm', 'mouse')):
        species = 'mouse'
    else:
        assert False, f'Species "{species}" not supported yet'

    url = f'/sequence/region/{species}/{chr_num}:{start}..{end}:1'
    response = _cached_session.get(
        'http://rest.ensembl.org' + url,
        params={'content-type': 'application/json'})
    json = response.json()
    if 'error' in json:
        raise ValueError(json['error'])
    response.raise_for_status()
    return json['seq']


def _get_cds_seq(record, cds: list, cds_index: int) -> str:
    cds_seq = cds[cds_index].location.extract(record).seq
    try:
        if cds_index == 0:
            # start codon
            assert cds_seq[0:3] == 'ATG', cds_seq
        elif cds_index == -1:
            # stop codon
            assert cds_seq[-3:] in ('TAG', 'TGA', 'TAA'), cds_seq
    except AssertionError:
        log.warning('Transcript does not have codon in expected location in CCDS "{}"'.format(cds_seq))
    return cds_seq


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    # fetch_ensembl_transcript('ENST00000398844').description
