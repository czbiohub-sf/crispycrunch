import os
import pickle
import pytest
import re

from httmock import HTTMock  # type: ignore
from httmock import all_requests  # type: ignore

from protospacex import fetch_ensembl_transcript


def _load_testing_data(ensembl_transcript_id: str):
    path = os.path.join(
        os.path.dirname(__file__),
        f'data/{ensembl_transcript_id}.pkl')
    with open(path, 'rb') as testing_file:
        testing_data = pickle.load(testing_file)
    return testing_data


@all_requests
def ensembl_mock(url, request):
    match = re.match(r'.*/(.*?)$', url.path)

    if match is None:
        raise ValueError(f"URL {url} doesn't end in an Ensembl ID. Is this a "
                         "valid Ensembl REST URL?")

    ensembl_transcript_id = match.group(1)

    try:
        testing_data = _load_testing_data(ensembl_transcript_id)
    except FileNotFoundError:
        raise ValueError(f"Testing file for transcript {ensembl_transcript_id} "
                         "not found")

    assert testing_data['transcript_name'] == ensembl_transcript_id

    print(f"Loaded test data for {ensembl_transcript_id}")
    print(f"Test data generated on {testing_data['time_generated']}")

    if 'sequence' in url.path:
        return testing_data['sequence_response']
    elif 'overlap' in url.path:
        return testing_data['overlap_response']
    else:
        raise NotImplementedError(f"Endpoint {url.path} not implemented for "
                                  "mocking")


@pytest.mark.parametrize("ensembl_transcript_id", [
    "ENST00000411809",
    "ENST00000398844"
])
def test_fetch_ensembl_transcript(ensembl_transcript_id):
    with HTTMock(ensembl_mock):
        testing_data = _load_testing_data(ensembl_transcript_id)

        record = fetch_ensembl_transcript(ensembl_transcript_id)

        exons = [f for f in record.features if f.type == "exon"]
        cdss = [f for f in record.features if f.type == "cds"]

        # We don't care about lower case vs. upper case bases
        UTR_5 = record.seq[exons[0].location.start:cdss[0].location.start]
        assert UTR_5.upper() == testing_data['UTR_5'].upper()

        UTR_3 = record.seq[cdss[-1].location.end:]
        assert UTR_3.upper() == testing_data['UTR_3'].upper()

        coding_seqs = [record.seq[cds.location.start:cds.location.end] for cds in cdss]

        assert coding_seqs[0].upper() == testing_data['first_coding_region'].upper()
        assert coding_seqs[-1].upper() == testing_data['last_coding_region'].upper()
        assert testing_data['another_coding_region'].upper() in coding_seqs

        exon_seqs = [record.seq[exon.location.start:exon.location.end] for exon in exons]

        assert exon_seqs[0].upper() == testing_data['first_exon'].upper()
        assert exon_seqs[-1].upper() == testing_data['last_exon'].upper()
        assert testing_data['another_exon'].upper() in exon_seqs

        random_intron = record.seq[testing_data['random_intron'].features[0]
                                   .location.start:testing_data['random_intron'].features[0].location.end + 1]
        assert random_intron.upper() == testing_data['random_intron'].seq.upper()
