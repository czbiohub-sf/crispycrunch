
"""
A collection of web clients that make requests to various websites. The clients
return data extracted from HTML. The clients may make multiple
dependent requests to get results. They may also retry in case of failure.

Server responses are cached by default using requests_cache.

Doctests will run slow on the first run before the cahce is warm.

See also SampleSheetTestCase for sample return data.
"""

import logging
import time
import urllib.parse

from abc import abstractmethod
# TODO (gdingle): OrderedDict no longer needed in python3.7
from collections import OrderedDict
from typing import Any, Dict

import requests
import requests_cache  # type: ignore
import urllib3

from bs4 import BeautifulSoup

logger = logging.getLogger(__name__)
_cached_session = requests_cache.CachedSession(
    cache_name=__name__ + '_cache',
    # TODO (gdingle): what's the best timeout?
    expire_after=3600 * 24 * 14,
    allowable_methods=('GET', 'POST'),
)
_cache = _cached_session.cache
# _cache.clear()

# NOTE: This monkey-patch is needed for a stable cache key for file uploads.
urllib3.filepost.choose_boundary = lambda: 'crispycrunch_super_special_form_boundary'

# CRISPOR_BASE_URL = 'http://crispor.tefor.net/crispor.py'
CRISPOR_BASE_URL = 'http://ec2-34-219-237-20.us-west-2.compute.amazonaws.com/crispor.py'


class AbstractScrapeRequest:

    def __repr__(self):
        return '{}({}, {})'.format(
            self.__class__, self.endpoint, self.__dict__.get('data'))

    def __str__(self):
        return self.__repr__()

    def in_cache(self) -> bool:
        return _cache.has_key(self.cache_key)

    @property
    def cache_key(self):
        return _cache.create_key(self.request)

    @abstractmethod
    def run(self) -> Dict[str, Any]:
        """Requests self.endpoint and extracts relevant data from the HTML response"""


class CrispressoRequest(AbstractScrapeRequest):
    """
    Requests a Crispresso2 analysis and returns the resulting report.
    Typically takes 90 seconds.

    >>> amplicon = 'cgaggagatacaggcggagggcgaggagatacaggcggagggcgaggagatacaggcggagagcgGCGCTAGGACCCGCCGGCCACCCCGCCGGCTCCCGGGAGGTTGATAAAGCGGCGGCGGCGTTTGACGTCAGTGGGGAGTTAATTTTAAATCGGTACAAGATGGCGGAGGGGGACGAGGCAGCGCGAGGGCAGCAACCGCACCAGGGGCTGTGGCGCCGGCGACGGACCAGCGACCCAAGCGCCGCGGTTAACCACGTCTCGTCCAC'
    >>> sgRNA = 'AATCGGTACAAGATGGCGGA'
    >>> fastq_r1 = '../crispresso/fastqs/A1-ATL2-N-sorted-180212_S1_L001_R1_001.fastq.gz'
    >>> fastq_r2 = '../crispresso/fastqs/A1-ATL2-N-sorted-180212_S1_L001_R2_001.fastq.gz'
    >>> req = CrispressoRequest(amplicon, sgRNA, fastq_r1, fastq_r2)
    >>> response = req.run()

    >>> len(response['report_files']) > 0
    True

    >>> len(response['report_stats']) > 0
    True

    >>> req.in_cache()
    True
    """

    # TODO (gdingle): parameterize
    base_url = 'http://ec2-52-12-22-81.us-west-2.compute.amazonaws.com'
    # TODO (gdingle): why is this broken currently?
    # base_url = 'http://crispresso.pinellolab.partners.org'

    def __init__(self,
                 amplicon: str,
                 sgRNA: str,
                 fastq_r1: str,
                 fastq_r2: str,
                 hdr_seq: str='',
                 optional_name: str='') -> None:
        self.endpoint = self.base_url + '/submit'
        self.data = {
            # NOTE: all post vars are required, even if empty
            'amplicon': amplicon,
            'amplicon_names': '',
            'be_from': 'C',
            'be_to': 'T',
            'demo_used': '',
            'email': 'gdingle@chanzuckerberg.com',
            'exons': '',
            'fastq_se': '',
            # TODO (gdingle): - [ ] Hook up crispresso2 to HDR
            # 'hdr_seq': 'cgaggagatacaggcggagggcgaggagatacaggcggagggcgaggagatacaggcggagagcgGCGCTAGGACCCGCCGGCCACCCCGCCGGCTCCCGGGAGGTTGATAAAGCGGCGGCGGCGTTTGACGTCAGTGGGGAGTTAATTTTAAATCGGTACAAGATGCGTGACCACATGGTCCTTCATGAGTATGTAAATGCTGCTGGGATTACAGGTGGCGGAttggaagttttgtttcaaggtccaggaagtggtGCGGAGGGGGACGAGGCAGCGCGAGGGCAGCAACCGCACCAGGGGCTGTGGCGCCGGCGACGGACCAGCGACCCAAGCGCCGCGGTTAACCACGTCTCGTCCAC',
            'hdr_seq': '',
            'optional_name': '',
            'optradio_exc_l': '15',
            'optradio_exc_r': '15',
            'optradio_hs': '60',
            'optradio_qc': '0',
            'optradio_qn': '0',
            'optradio_qs': '0',
            'optradio_trim': '',
            'optradio_wc': '-3',
            'optradio_ws': '1',
            'seq_design': 'paired',
            'sgRNA': sgRNA,
        }
        self.files = {
            'fastq_r1': open(fastq_r1, 'rb'),
            'fastq_r2': open(fastq_r2, 'rb'),
        }
        self.request = requests.Request(  # type: ignore
            'POST',
            self.endpoint,
            data=self.data,
            files=self.files,
        ).prepare()

    def run(self) -> Dict[str, Any]:
        logger.info('POST request to: {}'.format(self.endpoint))
        response = _cached_session.send(self.request)  # type: ignore
        response.raise_for_status()
        # for example: http://crispresso.pinellolab.partners.org/check_progress/P2S84K
        report_id = response.url.split('/')[-1]

        self._wait_for_success(report_id)

        report_data_url = self.base_url + '/reports_data/CRISPRessoRun{}'.format(report_id)
        report_files_url = '{}/CRISPResso_on_{}/'.format(report_data_url, report_id)
        report_zip = '{}/CRISPResso_Report_{}.zip'.format(report_data_url, report_id)
        report_url = self.base_url + '/view_report/' + report_id
        stats_url = report_files_url + 'CRISPResso_quantification_of_editing_frequency.txt'

        return {
            'report_url': report_url,
            'report_zip': report_zip,
            'log_params': self._get_log_params(report_url),
            'report_files': [report_files_url + file for file in self.report_files],
            'report_stats': self._get_stats(stats_url),
            'input_data': self.data,
            'input_files': [f.name for f in self.files.values()],
        }

    def _wait_for_success(self, report_id: str, retries: int = 3 * 96) -> None:
        """
        Poll for SUCCESS. Typically took 90 secs in testing.

        Crispresso2 appears to support running only 1 report in parallel, so in
        the worst case, 96 reports will take approx 3 hours.
        """
        total = 0
        while not self._check_report_status(report_id) and retries >= 0:
            amount = 30
            # amount = 90 // (retries + 1)  # backoff
            time.sleep(amount)
            retries -= 1
            total += amount

        if retries < 0:
            raise TimeoutError('Crispresso on {}: Retries exhausted after {}s.'.format(
                report_id, total))

    def _get_stats(self, stats_url: str) -> dict:
        logger.info('GET request to: {}'.format(stats_url))
        stats_response = _cached_session.get(stats_url)
        return self._parse_tsv(stats_response.text)

    @staticmethod
    def _parse_tsv(tsv: str) -> dict:
        r"""
        >>> tsv = '''Reference\tTotal\tUnmodified\tModified\tDiscarded\tInsertions\tDeletions\tSubstitutions\tOnly Insertions\tOnly Deletions\tOnly Substitutions\tInsertions and Deletions\tInsertions and Substitutions\tDeletions and Substitutions\tInsertions Deletions and Substitutions
        ... Reference\t14470\t12930\t1540\t0\t1537\t0\t6\t1534\t0\t3\t0\t3\t0\t0'''
        >>> CrispressoRequest._parse_tsv(tsv)
        OrderedDict([('Total', 14470), ('Unmodified', 12930), ('Modified', 1540), ('Discarded', 0), ('Insertions', 1537), ('Deletions', 0), ('Substitutions', 6), ('Only Insertions', 1534), ('Only Deletions', 0), ('Only Substitutions', 3), ('Insertions and Deletions', 0), ('Insertions and Substitutions', 3), ('Deletions and Substitutions', 0), ('Insertions Deletions and Substitutions', 0)])
        """
        lines = [line.split('\t')[1:] for line in tsv.split('\n')]
        headers = [h.strip() for h in lines[0]]
        values = [int(v.strip()) for v in lines[1]]
        return OrderedDict(zip(headers, values))

    def _get_log_params(self, report_url: str) -> str:
        logger.info('GET request to: {}'.format(report_url))
        report_response = _cached_session.get(report_url)
        soup = BeautifulSoup(report_response.text, 'html.parser')
        return soup.find(id='log_params').get_text()

    def _check_report_status(self, report_id: str) -> bool:
        status_endpoint = self.base_url + '/status/'
        status_url = status_endpoint + report_id
        logger.info('GET request to: {}'.format(status_url))
        # no cache here
        report_status = requests.get(status_url).json()
        if report_status['state'] == 'FAILURE':
            raise RuntimeError('Crispresso on {}: {}'.format(report_id, report_status['message']))
        elif report_status['state'] == 'SUCCESS':
            return True
        elif report_status['state'] == 'PENDING':
            # Sometimes pending status is never set to success though the report exists.
            # The bug was reported to Luca Pinello.
            report_url = self.base_url + '/view_report/' + report_id
            response = _cached_session.get(report_url)
            if response.status_code == 200:
                return True
            else:
                return False
        else:
            return False

    @property
    def report_files(self):
        return (
            'CRISPResso_RUNNING_LOG.txt',
            # TODO (gdingle): read pickle file and avoid "invalid start byte" error
            'CRISPResso2_info.pickle',  # contains figure captions
            'Alleles_frequency_table.txt',
            'CRISPResso_mapping_statistics.txt',
            'CRISPResso_quantification_of_editing_frequency.txt',
            'Mapping_statistics.txt',
            'Quantification_of_editing_frequency.txt',
            'Reference.Alleles_frequency_table_around_cut_site_for_AATCGGTACAAGATGGCGGA.txt',
            'Reference.deletion_histogram.txt',
            'Reference.effect_vector_combined.txt',
            'Reference.effect_vector_deletion.txt',
            'Reference.effect_vector_insertion.txt',
            'Reference.effect_vector_substitution.txt',
            'Reference.indel_histogram.txt',
            'Reference.insertion_histogram.txt',
            'Reference.modification_count_vectors.txt',
            'Reference.nucleotide_frequency_table.txt',
            'Reference.nucleotide_percentage_table.txt',
            'Reference.quantification_window_modification_count_vectors.txt',
            'Reference.quantification_window_nucleotide_frequency_table.txt',
            'Reference.quantification_window_nucleotide_percentage_table.txt',
            'Reference.quantification_window_substitution_frequency_table.txt',
            'Reference.substitution_frequency_table.txt',
            'Reference.substitution_histogram.txt',

            '1a.Read_Barplot.pdf',
            '1a.Read_Barplot.png',
            '1b.Alignment_Pie_Chart.pdf',
            '1b.Alignment_Pie_Chart.png',
            '1c.Alignment_Barplot.pdf',
            '1c.Alignment_Barplot.png',
            '2a.Reference.Nucleotide_Percentage_Quilt.pdf',
            '2a.Reference.Nucleotide_Percentage_Quilt.png',
            '2b.Reference.Nucleotide_Percentage_Quilt_For_AATCGGTACAAGATGGCGGA.pdf',
            '2b.Reference.Nucleotide_Percentage_Quilt_For_AATCGGTACAAGATGGCGGA.png',
            '3a.Reference.Indel_Size_Distribution.pdf',
            '3a.Reference.Indel_Size_Distribution.png',
            '3b.Reference.Insertion_Deletion_Substitutions_Size_Hist.pdf',
            '3b.Reference.Insertion_Deletion_Substitutions_Size_Hist.png',
            '4a.Reference.Combined_Insertion_Deletion_Substitution_Locations.pdf',
            '4a.Reference.Combined_Insertion_Deletion_Substitution_Locations.png',
            '4b.Reference.Insertion_Deletion_Substitution_Locations.pdf',
            '4b.Reference.Insertion_Deletion_Substitution_Locations.png',
            '4c.Reference.Quantification_Window_Insertion_Deletion_Substitution_Locations.pdf',
            '4c.Reference.Quantification_Window_Insertion_Deletion_Substitution_Locations.png',
            '4d.Reference.Position_Dependent_Average_Indel_Size.pdf',
            '4d.Reference.Position_Dependent_Average_Indel_Size.png',
            '9.Reference.Alleles_Frequency_Table_Around_Cut_Site_For_AATCGGTACAAGATGGCGGA.pdf',
            '9.Reference.Alleles_Frequency_Table_Around_Cut_Site_For_AATCGGTACAAGATGGCGGA.png',
        )


class CrisporGuideRequest(AbstractScrapeRequest):
    """
    Given a sequence or chromosome location, gets candidate Crispr guides from
    Crispor service.

    This module calls out to crispor.tefor.net endpoints. Crispor does not have
    an official public web API. You've been warned!

    >>> seq = 'chr1:11,130,540-11,130,751'
    >>> req = CrisporGuideRequest(name='test-crispr-guides', seq=seq)
    >>> data = req.run()

    >>> len(data['primer_urls']) > 3
    True

    >>> len(data['scores']) > 3
    True

    Responses are cached. Use in_cache to check status when many requests
    are in flight.

    >>> req.in_cache()
    True
    """

    def __init__(
            self,
            seq: str,
            name: str = '',
            org: str = 'hg38',
            pam: str = 'NGG',
            target: str = '',
            pre_filter: int = 2) -> None:

        self.data = {
            'name': name,
            'seq': seq,
            'org': org,
            'pam': pam,
            # TODO (gdingle): this is failing in the case of all low specificity guides
            # 'sortBy': 'offCount',
            # sort by number of off-targets
            'submit': 'SUBMIT',
        }
        self.endpoint = CRISPOR_BASE_URL
        self.request = requests.Request(  # type: ignore
            'POST', self.endpoint, data=self.data).prepare()
        self.target = target or seq
        self.pre_filter = pre_filter

    def run(self, retries: int=3) -> Dict[str, Any]:
        try:
            logger.info('POST request to: {}'.format(self.endpoint))
            response = _cached_session.send(self.request)  # type: ignore
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            return self._extract_data(soup, response.url)
        except TimeoutError as e:
            logger.warning(str(e))
            # IMPORTANT: Delete cache of "waiting" page
            _cache.delete(self.cache_key)
            if retries:
                time.sleep(60 // (retries + 1))  # backoff
                return self.run(retries - 1)
            else:
                raise
        except RuntimeError as e:
            logger.warning(str(e))
            # IMPORTANT: Delete cache of unexpected output
            _cache.delete(self.cache_key)
            raise
        raise RuntimeError('unknown error')

    def _extract_data(self, soup: BeautifulSoup, url: str) -> Dict[str, Any]:
        title = soup.find(class_='title')
        if title and 'not present in the selected genome' in title.get_text():
            raise ValueError('Crispor on {}: {}'.format(
                self.target, title.get_text()))

        if 'Input sequence range too long' in soup.get_text() or \
                'cannot handle sequences longer than' in soup.get_text():

            # This is a hack to re-route to the dev version of Crispor which has
            # a higher bp limit.
            # TODO (gdingle): unfortunately dev version does not have hg38!!!
            # TODO (gdingle): split up input into mulitple requests
            size = len(self.data['seq'])
            # if size <= 10000:
            #     self.endpoint = 'http://crispor-max.tefor.net/crispor.py'
            #     self.request = requests.Request('POST', self.endpoint, data=self.data).prepare()  # type: ignore
            #     raise TimeoutError('Large seq size {}. Retry on {}'.format(
            #         size, self.endpoint))
            # else:
            raise ValueError('Crispor on {}: Bad sequence size: {}'.format(
                self.target, size))

        if 'This page will refresh every 10 seconds' in soup.get_text():
            raise TimeoutError('Crispor on {}: Stuck in job queue. Please retry.'.format(
                self.target))

        output_table = soup.find('table', {'id': 'otTable'})
        if not output_table:
            # TODO (gdingle): convert these to exceptions to be handled on progress page?
            if 'Found no possible guide sequence' in soup.get_text():
                return dict(
                    target=self.target,
                    guide_seqs={
                        'not found': 'not found',
                        # TODO (gdingle): make url on error work somehow
                        # 'url': url,
                    },
                )
            if 'Server error: could not run command' in soup.get_text():
                # Delete because intermittent
                _cache.delete(self.cache_key)
                return dict(
                    target=self.target,
                    guide_seqs={
                        'server error': 'server error',
                        # TODO (gdingle): make url on error work somehow
                        # 'url': url,
                    },
                )
            if 'are not valid in the genome' in soup.get_text():
                return dict(
                    target=self.target,
                    guide_seqs={
                        'invalid chromosome range': 'invalid chromosome range',
                        # TODO (gdingle): make url on error work somehow
                        # 'url': url,
                    },
                )
            if 'An error occured during processing' in soup.get_text():
                raise RuntimeError(
                    'Crispor: An error occured during processing: {}'.format(self.data['seq']))

            raise RuntimeError('Crispor on {}: No output rows. "{}"'.format(
                self.target, soup.find('body')))

        batch_id = soup.find('input', {'name': 'batchId'})['value']
        url = self.endpoint + '?batchId=' + batch_id
        primers_url = self.endpoint.split('?')[0] + '?batchId={}&pamId={}&pam=NGG'

        rows = [
            [t['id']] +
            [cell for cell in t.find_all('td')[1:8]] +
            [primers_url.format(batch_id, urllib.parse.quote(t['id']))]
            for t in output_table.find_all(class_='guideRow')]

        # Filter for rows that have possible primers and better than low specificity
        # See http://crispor.tefor.net/manual/
        rows = [r for r in rows
                if 'primers' in r[1].get_text() and
                r[2].get_text().strip().isdigit() and
                int(r[2].get_text().strip()) > 20]

        # Filter for rows that have actual primers (by http request),
        # but only first `pre_filter` to save time.
        # TODO (gdingle): pool.map parallelize? Or are we already parallel enough at this point?
        if self.pre_filter:
            rows = [r for i, r in enumerate(rows)
                    if i > self.pre_filter or
                    'Warning: No primers were found'
                    not in _cached_session.get(r[-1]).text]
            if not rows:
                return dict(
                    target=self.target,
                    guide_seqs={'not found': 'not found'},
                )

        # TODO (gdingle): refactor to simple lists... see GuideDesign.to_df
        scores = OrderedDict((t[0],
                              [c.get_text().strip() for c in t[2:5]])
                             for t in rows)
        scores = OrderedDict((t[0],
                              [c.get_text().strip() for c in t[2:5]])
                             for t in rows)
        guide_seqs = OrderedDict((t[0], t[1].find_next('tt').get_text())
                                 for t in rows)
        primer_urls = OrderedDict((t[0],
                                   primers_url.format(batch_id, urllib.parse.quote(t[0])))
                                  for t in rows)

        return dict(
            target=self.target,
            url=url,
            batch_id=batch_id,
            guide_seqs=guide_seqs,
            scores=scores,
            primer_urls=primer_urls,
            # TODO (gdingle): are these links ever needed?
            fasta_url=self.endpoint + '?batchId={}&download=fasta'.format(batch_id),
            benchling_url=self.endpoint + '?batchId={}&download=benchling'.format(batch_id),
            guides_url=self.endpoint + '?batchId={}&download=guides&format=tsv'.format(batch_id),
            offtargets_url=self.endpoint + \
            '?batchId={}&download=offtargets&format=tsv'.format(batch_id),
        )


class CrisporGuideRequestByBatchId(CrisporGuideRequest):
    """
    For debugging Crispor errors.
    """

    def __init__(self, batch_id: str, pre_filter: int = 5) -> None:
        self.endpoint = CRISPOR_BASE_URL + '?batchId=' + batch_id
        self.request = requests.Request('GET', self.endpoint).prepare()  # type: ignore
        self.target = batch_id
        self.pre_filter = pre_filter


class CrisporPrimerRequest(AbstractScrapeRequest):
    """
    Gets primers for a guide generated by Crispor.

    NOTE: Crispor uses Primer3 under the covers.

    NOTE: These test came from
    CRISPOR_BASE_URL = 'http://ec2-34-219-237-20.us-west-2.compute.amazonaws.com/crispor.py'

    >>> req = CrisporPrimerRequest('gYvicTzp9e5VPFC9YwLR', 's45-')
    >>> data = req.run()
    >>> len(data['ontarget_primers']) == 2
    True

    >>> req.in_cache()
    True

    With HDR.

    >>> req = CrisporPrimerRequest('gYvicTzp9e5VPFC9YwLR', 's45-', hdr_dist=0)
    >>> data = req.run()
    >>> data['ontarget_primers']
    {}

    Test guide #2 has HDR primers.

    >>> req = CrisporPrimerRequest('gYvicTzp9e5VPFC9YwLR', 's11-', hdr_dist=0)
    >>> data = req.run()
    >>> len(data['ontarget_primers']) == 2
    True
    """

    def __init__(
            self,
            batch_id: str,
            pam_id: str,
            amp_len: int=250,
            tm: int=60,
            pam: str='NGG',
            target: str='',
            hdr_dist: int = None) -> None:

        self.pam_id = pam_id
        quoted_pam_id = urllib.parse.quote(pam_id)  # percent encode the '+' symbol

        self.endpoint = CRISPOR_BASE_URL
        self.endpoint += f'?ampLen={amp_len}&tm={tm}&batchId={batch_id}&pamId={quoted_pam_id}&pam={pam}'
        if hdr_dist is not None:
            self.endpoint += f'&hdrDist={hdr_dist}'

        self.target = target  # just for metadata
        self.request = requests.Request('GET', self.endpoint).prepare()  # type: ignore

    def __repr__(self):
        return 'CrisporPrimerRequest({})'.format(self.endpoint)

    def run(self,
            retries: int=1) -> dict:
        # TODO (gdingle): temp for working on crispor
        # _cache.delete(self.cache_key)
        try:
            logger.info('GET request to: {}'.format(self.endpoint))
            response = _cached_session.send(self.request)  # type: ignore
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            return self._extract_data(soup)
        except RuntimeError as e:
            logger.warning(str(e))
            _cache.delete(self.cache_key)
            if retries:
                return self.run(retries - 1)
            else:
                raise

    def _extract_data(self, soup: BeautifulSoup) -> dict:
        if soup is None:
            raise RuntimeError('Cannot parse HTML {}'.format(soup))

        if 'exceptions.ValueError' in soup.get_text():
            raise RuntimeError('Crispor exceptions.ValueError')

        if 'Error:' in soup.get_text():
            raise RuntimeError('Crispor error: {}'.format(
                soup.get_text().split('Error:')[1].strip()))

        return dict(
            pam_id=self.pam_id,
            target=self.target,
            url=self.endpoint,
            # TODO (gdingle): this is broken sometimes
            # amplicon_length=soup
            # .find('select', {'name': 'ampLen'})
            # .find('option', {'selected': 'selected'})['value'],
            # TODO (gdingle): this is broken sometimes
            # primer_temp=soup
            # .find('select', {'name': 'tm'})
            # .find('option', {'selected': 'selected'})['value'],
            # TODO (gdingle): may not need primer_seqs
            # primer_tables=primer_tables,
            # primer_seqs=primer_seqs,
            ontarget_primers=self._extract_ontarget_primers(soup),
        )

    def _extract_ontarget_primers(self, soup: BeautifulSoup) -> Dict[str, tuple]:
        ontargetPcr = soup.find(id='ontargetPcr')
        table = ontargetPcr.find_next(class_='primerTable')
        message = ontargetPcr.find_next('strong')
        if message and 'Warning' in message.get_text():
            return {}  # will be interpreted as 'not found'
            # TODO (gdingle): better to raise exception?
            # raise ValueError('Cripor at {}: "{}"'.format(
            #     self.endpoint,
            #     message.get_text()))

        if table is None:
            # TODO (gdingle): catch crispor execptions here?
            text = ontargetPcr.find_next('div').get_text()
            if 'No perfect match found' in text:
                raise ValueError('Cripor at {}: "{}"'.format(
                    self.endpoint,
                    'No perfect match found for guide sequence in the genome. Cannot design primers for a non-matching guide sequence.'))
            else:
                raise ValueError('Cripor at {}: "{}"'.format(self.endpoint, text))

        rows = (row.find_all('td') for row in table.find_all('tr'))  # get primers
        tts = ontargetPcr.find_next('div').find_all('tt')  # get products
        # TODO (gdingle): figure out why both products are sometimes not grabbed
        # see http://crispor.tefor.net/crispor.py?ampLen=400&tm=60&batchId=oI0f65TEwlZdrH9NMAat&pamId=s97-&pam=NGG
        # Temp fix because we know both are always the same
        tts = (list(tts)[0], list(tts)[0])
        return dict(
            (row[0].get_text().split('_')[-1],
                (row[1].get_text(), tt.get_text())
             )
            for row, tt in zip(rows, tts))


if __name__ == '__main__':
    # import doctest  # noqa
    # doctest.testmod()

    # seq = 'chrX:49250226-49250297'
    seq = 'AACGCTGGCCTCCTAGGCCGGGTCCGGGAGGCCTGAGGAGCCGCCGGCAGAGGTCTCTCCCCAGCCTCAGGC'
    req = CrisporGuideRequest(name='test-crispr-guides', seq=seq)
    data = req.run()
    print(data)

    # req = CrisporGuideRequestByBatchId('tZgMsg3spbVL3Irgvhvl', pre_filter=5)
    # data = req.run()
    # print(data['guide_seqs'])

    # req = CrisporPrimerRequest('oI0f65TEwlZdrH9NMAat', 's97-')
    # data = req.run()
    # print(data)
