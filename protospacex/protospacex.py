"""
This module was extracted from https://github.com/czbiohub/protospacex,
which was designed for https://czi.quip.com/YbAhAbOV4aXi/.

protospacex: automated guide design for Cas9 knock-in experiments

Given an Ensembl transcript id, return target regions for Crispr guides.
"""
import logging
import requests

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqFeature import FeatureLocation  # type: ignore
from Bio.SeqFeature import SeqFeature  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def fetch_ensembl_transcript(ensembl_transcript_id: str) -> SeqRecord:
    """Fetch the requested Ensembl transcript.

    # TODO (gdingle): do we actually need exon features? or just cds? #"exon"
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
    # TODO (gdingle): do we actually need exon features? or just cds? #"exon"
      with exon and CDS features. The coordinates of exons and CDS
      features are relative to the sequence fragment.

    >>> fetch_ensembl_transcript('ENST00000398844').description
    'chromosome:GRCh38:5:134648789:134727823:1'
    """

    # TODO: Validate ensembl_transcript_id is a valid transcript id

    base_url = "http://rest.ensembl.org"

    # First, fetch the transcript sequence
    url = base_url + f"/sequence/id/{ensembl_transcript_id}"

    log.info(f"Querying Ensembl for sequence of {ensembl_transcript_id}")
    response = requests.get(url, {"type": "genomic",
                                  "content-type": "application/json"})

    try:
        response.raise_for_status()
    except requests.HTTPError:
        log.error("Ensembl sequence REST query returned error "
                  "{}".format(response.text))
        raise ValueError(response.text)

    response_data = response.json()

    try:
        description = response_data['desc'].split(':')
        species = description[1]
        chromosome_number = description[2] # may be X
        sequence_left = int(description[3])
        sequence_right = int(description[4])
        transcript_strand = int(description[5])

        if sequence_left > sequence_right:
            raise ValueError(f"Expected left sequence boundary {sequence_left} "
                             f"<= right sequence boundary {sequence_right}: did "
                             "the format of the Ensembl REST response change?")

        sequence_id = response_data['id']

        seq_str = response_data['seq']

        log.info(f"Retrieved sequence {response_data['desc']} of length "
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

    log.info(f"Querying Ensembl for overlaps of {ensembl_transcript_id}")
    response = requests.get(url, {"feature": ["cds", "exon"],
                                  "content-type": "application/json"})
    try:
        response.raise_for_status()
    except requests.HTTPError:
        log.error("Ensembl sequence REST query returned error "
                  "{}".format(response.text))
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

        log.info(f"Retrieved {num_exon_boundaries} exons and "
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
    'ATGTCCCAGCCGGGAATACCGGCCTCCGGCGGCGCC'

    >>> get_cds_seq('ENST00000411809')
    'ATGTTGAACATGTGGAAGGTGCGCGAGCTGGTGGAC'

    Stop codon.

    >>> get_cds_seq('ENST00000398844', -1)
    'CAACAAGTGAATAAATGAATGAATGAAGAAATTTGA'

    >>> get_cds_seq('ENST00000411809', -1)
    'GCCAATTTTAGCAAATAAGAGATTGTAAAAGAAGCA'

    Length 30.
    >>> len(get_cds_seq('ENST00000221801', -1, length=30))
    30
    """
    record = fetch_ensembl_transcript(ensembl_transcript_id)
    cds = [f for f in record.features if f.type == 'cds']
    assert len(cds)
    location = cds[cds_index].location

    # enforce length
    start, end, codon_at = _get_start_end(
        location,
        length,
        cds_index
    )

    cds_seq = record.seq[start:end]

    assert len(cds_seq) == length
    assert len(cds_seq) % 3 == 0, 'must be codon aligned'

    if cds_index == 0:
        assert cds_seq[codon_at:codon_at + 3] == 'ATG'
    elif cds_index == -1:
        assert cds_seq[codon_at:codon_at + 3] in ['TAG', 'TGA', 'TAA']

    return str(cds_seq)


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

    See https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000113615;r=5:134648789-134727823;t=ENST00000398844
    See https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS43363

    >>> get_cds_chr_loc('ENST00000398844', length=990)
    'chr5:134649077-134650066'

    >>> get_cds_chr_loc('ENST00000411809', length=990)
    'chr5:157786494-157787483'

    Length.

    >>> get_cds_chr_loc('ENST00000221801', length=990)
    'chr19:39834572-39835561'

    >>> get_cds_chr_loc('ENST00000398844', -1, 90)
    'chr5:134725050-134725139'

    Get last codon.

    >>> get_cds_chr_loc('ENST00000398844', -1, 990)
    'chr5:134724600-134725589'

    >>> get_cds_chr_loc('ENST00000411809', -1, 990)
    'chr5:157857324-157858313'

    Length, last codon.

    >>> get_cds_chr_loc('ENST00000221801', -1, length=30)
    'chr19:39846320-39846349'

    Length.

    >>> get_cds_chr_loc('ENST00000398844')
    'chr5:134649077-134649112'

    >>> get_cds_chr_loc('ENST00000398844', -1)
    'chr5:134725077-134725112'
    """
    record = fetch_ensembl_transcript(ensembl_transcript_id)
    cds = [f for f in record.features if f.type == 'cds']
    assert len(cds)
    cds_location = cds[cds_index].location
    cds_seq = cds[cds_index].location.extract(record).seq
    # TODO (gdingle): why isn't cds_seq always divisible by 3?
    # assert len(cds_seq) % 3 == 0, len(cds_seq)
    if cds_index == 0:
        # start codon
        assert cds_seq[0:3] == 'ATG', cds_seq
    elif cds_index == -1:
        # stop codon
        assert cds_seq[-3:] in ('TAG', 'TGA', 'TAA'), cds_seq

    species = record.annotations['reference_species']
    chromosome_number = record.annotations['reference_chromosome_number']
    sequence_left = record.annotations['reference_left_index']

    # Transcript of reverse strand is translated above.
    assert cds_location.strand == 1
    assert species == 'GRCh38'

    # TODO (gdingle): exon junctions
    # * Have a filter/flag for intron/exon junctions: do not cut or mutate less than 3 nt away

    start = sequence_left + cds_location.start
    end = sequence_left + cds_location.end
    assert end - start == len(cds_seq)

    # enforce length
    start, end, codon_at = _get_start_end(
        FeatureLocation(start, end),
        length,
        cds_index
    )

    return 'chr{}:{}-{}'.format(
        chromosome_number,
        start,
        end - 1,  # change to inclusive range
    )


def _get_start_end(
        location: FeatureLocation,
        length: int,
        cds_index: int) -> tuple:
    assert length > 0

    def divisible(i: int) -> bool:
        return i % 3 == 0 and i % 2 == 0

    if not divisible(length):
        raise ValueError(
            f'length {length} must be divisible by 3 and 2 to ensure codon frame and symmetry')

    start = location.start
    end = location.end
    if cds_index == -1:  # stop codon
        start = end - length // 2
        codon_at = end - 3
        end = end + length // 2
    elif cds_index == 0:  # start codon
        end = start + length
        codon_at = start
    assert end - start == length
    return start, end, codon_at - start


if __name__ == '__main__':
    import requests_cache  # type: ignore
    requests_cache.install_cache()

    import doctest
    doctest.testmod()

    # s = get_cds_chr_loc('ENST00000221801', -1, 990, 0)
    # print(s)
    # s = get_cds_seq('ENST00000221801', -1)
    # print(s, len(s))
