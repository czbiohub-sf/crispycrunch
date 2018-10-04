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
        sequence_left = int(description[3])
        sequence_right = int(description[4])
        transcript_strand = int(description[5])

        if sequence_left > sequence_right:
            raise ValueError(f"Expected left sequence boundary {sequence_left} "
                             f"<= right sequence boundary {sequence_right}: did "
                             "the format of the Ensembl REST response change?")

        sequence_id = response_data['id']

        seq_str = response_data['seq']

        log.info(f"Retrieved sequence of length {sequence_right - sequence_left} "
                 f"for species {species} on strand {transcript_strand}")
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
    # TODO (gdingle): do we actually need exon features? or just cds? #"exon"
    response = requests.get(url, {"feature": ["cds"],
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
            # log.info(response_datum)
            record.features.append(SeqFeature(
                location=FeatureLocation(
                    int(response_datum['start']) - sequence_left,
                    int(response_datum['end']) - sequence_left + 1,
                    strand=int(response_datum['strand'])),
                type=response_datum['feature_type']))
            # log.info(record.features[-1])
            # if response_datum['feature_type'] == 'exon' and response_datum['rank'] == 1:
            #    log.info(seq[response_datum['start'] - sequence_left:response_datum['end'] - sequence_left])
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
    # log.info(sequence_left)
    # log.info(sequence_right)
    # log.info(len(seq))

    if transcript_strand == -1:
        record = record.reverse_complement()

    # Finally, sort features by their start locations
    record.features.sort(key=lambda f: f.location.start)

    return record


def start_codon_seq(ensembl_transcript_id: str) -> str:
    """
    Convience function to return string sequence of start codon.

    See https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000113615;r=5:134648789-134727823;t=ENST00000398844
    See https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS43363

    >>> start_codon_seq('ENST00000398844')
    'ATGTCCCAGCCGGGAATACCGGCCTCCGGCGGCGCCCCAGCCAGCCTCCAGGCCCAGAACGGAGCCGCCTTGGCCTCGGGGTCTCCCTACACCAACG'
    """
    record = fetch_ensembl_transcript(ensembl_transcript_id)
    cds = [f for f in record.features if f.type == 'cds']
    assert len(cds)
    start_codon_seq = cds[0].location.extract(record).seq
    assert len(start_codon_seq)
    return str(start_codon_seq)


def start_codon_chr_loc(ensembl_transcript_id: str) -> str:
    """
    Convience function to return chromosome location of start codon.

    See https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?g=ENSG00000113615;r=5:134648789-134727823;t=ENST00000398844
    See https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS43363

    >>> start_codon_chr_loc('ENST00000398844')
    'chr5:134649077-134649174'

    # TODO (gdingle): upgrade when able
    >> start_codon_chr_loc('ENST00000398844')
    'GRCh38:chr5:134649077-134649174:+'
    """
    record = fetch_ensembl_transcript(ensembl_transcript_id)
    cds = [f for f in record.features if f.type == 'cds']
    assert len(cds)
    codon_location = cds[0].location
    codon_seq = cds[0].location.extract(record).seq
    assert codon_seq[0:3] == 'ATG', codon_seq

    description = record.description.split(':')
    assert len(description) == 6
    # TODO (gdingle): deal with other strands
    assert codon_location.strand == 1
    # TODO (gdingle): deal with other genomes
    assert description[1] == 'GRCh38'

    start = int(description[3]) + codon_location.start
    end = int(description[3]) + codon_location.end
    assert end - start == len(codon_seq)
    return 'chr{}:{}-{}'.format(
        description[2],
        start,
        end,
    )

    # return '{}:chr{}:{}-{}:{}'.format(
    #     description[1],
    #     description[2],
    #     int(description[3]) + codon_location.start,
    #     int(description[3]) + codon_location.end,
    #     '+' if codon_location.strand == 1 else '-',
    # )


if __name__ == '__main__':
    import doctest
    doctest.testmod()
