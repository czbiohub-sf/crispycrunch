
"""
A collection of web clients that make requests to various websites. The clients
return data extracted from the websidote HTML. The clients may make multiple
dependent requests to get results. They may also retry in case of failure.

Server responses are cached by default using requests_cache.

Doctests will run slow on the first run before the cahce is warm.

See also SampleSheetTestCase for sample return data.
"""
import json
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

CACHE_ARGS = dict(
    cache_name=__name__ + '_cache',
    # TODO (gdingle): what's the best timeout?
    expire_after=3600 * 24 * 14,
    allowable_methods=('GET', 'POST'),
)
# TODO (gdingle): IMPORTANT! Things seem to break if we install the cache simulatenously
# in other modules.
# TODO (gdingle): use one cache session per module
requests_cache.install_cache(**CACHE_ARGS)
CACHE = requests_cache.core.get_cache()
# CACHE.clear()

# NOTE: This monkey-patch is needed for a stable cache key for file uploads.
urllib3.filepost.choose_boundary = lambda: 'crispycrunch_super_special_form_boundary'


class AbstractScrapeRequest:

    def __repr__(self):
        return '{}({}, {})'.format(
            self.__class__, self.endpoint, self.__dict__.get('data'))

    def __str__(self):
        return self.__repr__()

    def in_cache(self) -> bool:
        return CACHE.has_key(self.cache_key)

    @property
    def cache_key(self):
        return CACHE.create_key(self.request)

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

    base_url = 'http://crispresso.pinellolab.partners.org'

    def __init__(self,
                 amplicon: str,
                 sgRNA: str,
                 fastq_r1: str,
                 fastq_r2: str,
                 hdr_seq: str='',
                 optional_name: str='') -> None:
        self.endpoint = self.base_url + '/submit'
        # TODO (gdingle): compare crispresso2 values to chosen crispresso1 values
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
            # TODO (gdingle): why is this not working?
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
        # This is necessary because of multiple threads and .disabled() call below
        # TODO (gdingle): is this thread-safe enough?
        with requests_cache.enabled(**CACHE_ARGS):
            response = requests.Session().send(self.request)  # type: ignore
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
        # TODO (gdingle): upgrade crispresso somehow
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
        stats_response = requests.get(stats_url)
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
        report_response = requests.get(report_url)
        soup = BeautifulSoup(report_response.text, 'html.parser')
        return soup.find(id='log_params').get_text()

    def _check_report_status(self, report_id: str) -> bool:
        status_endpoint = self.base_url + '/status/'
        status_url = status_endpoint + report_id
        with requests_cache.disabled():
            logger.info('GET request to: {}'.format(status_url))
            report_status = requests.get(status_url).json()
        if report_status['state'] == 'FAILURE':
            raise RuntimeError('Crispresso on {}: {}'.format(report_id, report_status['message']))
        elif report_status['state'] == 'SUCCESS':
            return True
        elif report_status['state'] == 'PENDING':
            # TODO (gdingle): report bug to Luca Pinello
            # Sometimes pending status is never set to success though the report exists.
            report_url = self.base_url + '/view_report/' + report_id
            response = requests.get(report_url)
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
            pre_filter: int = 20) -> None:

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
        self.endpoint = 'http://crispor.tefor.net/crispor.py'
        self.request = requests.Request('POST', self.endpoint, data=self.data).prepare()  # type: ignore
        self.target = target or seq
        self.pre_filter = pre_filter

    def run(self, retries: int=3) -> Dict[str, Any]:
        try:
            logger.info('POST request to: {}'.format(self.endpoint))
            response = requests.Session().send(self.request)  # type: ignore
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            return self._extract_data(soup, response.url)
        except TimeoutError as e:
            logger.warning(str(e))
            # IMPORTANT: Delete cache of "waiting" page
            CACHE.delete(self.cache_key)
            if retries:
                time.sleep(60 // (retries + 1))  # backoff
                return self.run(retries - 1)
            else:
                raise
        except RuntimeError as e:
            logger.warning(str(e))
            # IMPORTANT: Delete cache of unexpected output
            CACHE.delete(self.cache_key)
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
                raise RuntimeError('Crispor: An error occured during processing.')

            raise RuntimeError('Crispor on {}: No output rows. "{}"'.format(
                self.target, soup.find('body').get_text().strip()))

        batch_id = soup.find('input', {'name': 'batchId'})['value']
        url = self.endpoint + '?batchId=' + batch_id
        primers_url = self.endpoint.split('?')[0] + '?batchId={}&pamId={}&pam=NGG'

        # TODO (gdingle): prefilter scores, guides, primer_urls
        # by existence of primer:
        # 'Warning: No primers were found'
        # Use new threadpool? limit to first 5?
        # TODO (gdingle): check for iteration strat of yunfang

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
        # TODO (gdingle): pool.map parallelize
        if self.pre_filter:
            rows = [r for i, r in enumerate(rows)
                    if i > self.pre_filter or
                    'Warning: No primers were found'
                    not in requests.get(r[-1]).text]
            if not rows:
                return dict(
                    target=self.target,
                    guide_seqs={'not found': 'not found'},
                )

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
            # TODO (gdingle): are these links ever needed?
            primer_urls=primer_urls,
            fasta_url=self.endpoint + '?batchId={}&download=fasta'.format(batch_id),
            benchling_url=self.endpoint + '?batchId={}&download=benchling'.format(batch_id),
            guides_url=self.endpoint + '?batchId={}&download=guides&format=tsv'.format(batch_id),
            offtargets_url=self.endpoint + '?batchId={}&download=offtargets&format=tsv'.format(batch_id),
        )


class CrisporGuideRequestByBatchId(CrisporGuideRequest):
    """
    For debugging Crispor errors.
    """

    def __init__(self, batch_id: str, pre_filter: int = 5) -> None:
        self.endpoint = 'http://crispor.tefor.net/crispor.py?batchId=' + batch_id
        self.request = requests.Request('GET', self.endpoint).prepare()  # type: ignore
        self.target = batch_id
        self.pre_filter = pre_filter


class CrisporPrimerRequest(AbstractScrapeRequest):
    """
    Gets primers for a guide generated by Crispor.

    NOTE: Crispor uses Primer3 under the covers.

    >>> req = CrisporPrimerRequest('9cJNEsbfWiSKa8wlaJMZ', 's185+')
    >>> data = req.run()
    >>> len(data['ontarget_primers']) == 2
    True

    >>> req.in_cache()
    True
    """

    def __init__(
            self,
            batch_id: str,
            pam_id: str,
            amp_len: int=250,
            tm: int=60,
            pam: str='NGG',
            target: str='') -> None:

        self.pam_id = pam_id
        quoted_pam_id = urllib.parse.quote(pam_id)  # percent encode the '+' symbol
        self.endpoint = 'http://crispor.tefor.net/crispor.py' + \
            '?ampLen={amp_len}&tm={tm}&batchId={batch_id}&pamId={quoted_pam_id}&pam={pam}'.format(**locals())
        self.target = target  # just for metadata
        self.request = requests.Request('GET', self.endpoint).prepare()  # type: ignore

    def __repr__(self):
        return 'CrisporPrimerRequest({})'.format(self.endpoint)

    def run(self,
            retries: int=1) -> dict:
        try:
            logger.info('GET request to: {}'.format(self.endpoint))
            response = requests.Session().send(self.request)  # type: ignore
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            return self._extract_data(soup)
        except RuntimeError as e:
            logger.warning(str(e))
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
        if message and 'Warning: No primers were found' in message.get_text():
            return {}  # will be interpreted as 'not found'
            # TODO (gdingle): better to raise exception?
            # raise ValueError('Cripor at {}: "{}"'.format(
            #     self.endpoint,
            #     message.get_text()))

        if table is None:
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

# TODO (gdingle): remove me when protospacex is ready


class TagInRequest(AbstractScrapeRequest):
    """
    Given an an Ensembl Transcript or a custom sequence, gets candidate sgRNAs
    and ssDNAs from tagin.stembio.org service for HDR experiments.

    To bypass CSRF protection, we make a initial request to extract the CSRF
    token out of the form, so we can submit it with the main request.

    The sessionid is also needed, for unknown reasons.

    TODO: re-enable doctest
    data = TagInRequest('ENST00000330949').run()
    > len(data['guide_seqs']) >= 1
    True
    > len(data['donor_seqs']) >= 1
    True
    > all(s in data['guide_seqs'].values() for s in data['donor_seqs'].keys())
    True
    """

    # TODO (gdingle): convert __init__ from acc_number to seq and compute stop codon position
    # https://www.genscript.com/tools/codon-frequency-table
    # TAG
    # TAA
    # TGA

    def __init__(self, acc_number: str, tag: str='FLAG', species: str='GRCh38') -> None:
        self.data = {
            'acc_number': acc_number,
            'tag': tag,
            'species': species,
        }
        self.endpoint = 'http://tagin.stembio.org/submit/'

    def run(self) -> Any:
        self._init_session()
        url = self.endpoint + self._get_design_id()

        cookies = {
            # Both of these are necessary!
            'csrftoken': self.csrftoken,
            'sessionid': self.sessionid,
        }
        logger.info('GET request to: {}'.format(url))
        response = requests.get(url, cookies=cookies)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        return self._extract_data(soup, url)

    def _init_session(self) -> None:
        logger.info('GET request to: {}'.format(self.endpoint))
        initial = requests.get(self.endpoint)
        initial.raise_for_status()
        self.csrftoken = initial.cookies['csrftoken']

    def _get_design_id(self) -> str:
        self.data['csrfmiddlewaretoken'] = self.csrftoken
        headers = {
            'X-Requested-With': 'XMLHttpRequest',
            'Cookie': 'csrftoken={}'.format(self.csrftoken),
        }

        logger.info('POST request to: {}'.format(self.endpoint))
        json_response = requests.post(self.endpoint, data=self.data, headers=headers)
        self.sessionid = json_response.cookies['sessionid']

        json_response.raise_for_status()
        assert json_response.json()['Success']

        return json_response.json()['Id']

    def _extract_data(self, soup: BeautifulSoup, url: str) -> Dict[str, Any]:
        # TODO (gdingle): are these useful?
        # sgRNA = soup.find(id='table_id1').get_text()
        # off_targets = soup.find(id='table_id2').get_text()

        data_user = [json.loads(s.attrs['data-user'])
                     for s in soup.findAll('div')
                     if s.has_attr('data-user')]
        top_guides = sorted(data_user[1], key=lambda d: d['sgRNA_score'], reverse=True)
        # TODO (gdingle): remap to offset position key used by Crispor
        guide_seqs = OrderedDict((round(d['sgRNA_score'], 2), d['sgRNA']) for d in top_guides)

        rows = [row.findAll('td') for row in soup.find(id='table_id3').findAll('tr')]

        def valid_row(row) -> bool:
            return bool(len(row) and row[0] and row[1] and
                        row[1].get_text().strip() != 'sgRNA too far from stop codon')

        donor_seqs = OrderedDict(
            (row[0].get_text(), row[1].get_text())
            for row in rows
            if valid_row(row) and row[0].get_text() in guide_seqs.values()
        )

        metadata = data_user[0][0]
        metadata['chr_loc'] = 'chr{}:{}-{}'.format(metadata['chrm'], metadata['tx_start'], metadata['tx_stop'])
        # This is the narrower range that contains all sgRNA guides.
        metadata['guide_chr_range'] = 'chr{}:{}-{}'.format(
            metadata['chrm'],
            min(g['sgRNA_start'] for g in data_user[1]),
            max(g['sgRNA_stop'] for g in data_user[1]))

        return {
            'guide_seqs': guide_seqs,
            'donor_seqs': donor_seqs,
            'metadata': metadata,
            'url': url,
        }


if __name__ == '__main__':
    import doctest  # noqa
    doctest.testmod()

    # seq = 'chr2:150500625-150500725'
    # req = CrisporGuideRequest(name='test-crispr-guides', seq=seq)
    # data = req.run()
    # print(data)

    # req = CrisporGuideRequestByBatchId('tZgMsg3spbVL3Irgvhvl', pre_filter=5)
    # data = req.run()
    # print(data['guide_seqs'])

    # req = CrisporPrimerRequest('oI0f65TEwlZdrH9NMAat', 's97-')
    # data = req.run()
    # print(data)
