"""
Crispor Client.

This module calls out to crispor.tefor.net endpoints. Crispor does not have
an official public web API. You've been warned!
"""
import requests
import requests_cache

import json
import time

from abc import abstractmethod
from bs4 import BeautifulSoup
from collections import OrderedDict
from typing import Any, Dict, Tuple
from urllib.parse import quote

requests_cache.install_cache(
    # TODO (gdingle): what's the best timeout?
    expire_after=3600,
    allowable_methods=('GET', 'POST'))
CACHE = requests_cache.core.get_cache()
# CACHE.clear()


class AbstractCrisporRequest:

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__, self.endpoint, self.__dict__.get('data'))

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


class CrisporGuideRequest(AbstractCrisporRequest):
    """
    Given a sequence or chromosome location, gets candidate Crispr guides from
    Crispor service.

    >>> seq = 'chr1:11,130,540-11,130,751'
    >>> req = CrisporGuideRequest(name='test-crispr-guides', seq=seq)
    >>> data = req.run()
    >>> len(data['primer_urls']) > 3
    True

    Responses are cached. Use in_cache to check status when many requests
    are in flight.

    >>> req.in_cache()
    True
    """

    def __init__(self, seq: str, name: str = '', org: str = 'hg19', pam: str = 'NGG') -> None:
        self.data = {
            'name': name,
            'seq': seq,
            'org': org,
            'pam': pam,
            'sortBy': 'spec',
            'submit': 'SUBMIT',
        }
        self.endpoint = 'http://crispor.tefor.net/crispor.py'
        self.request = requests.Request('POST', self.endpoint, data=self.data).prepare()

    def run(self, retries: int=3) -> Dict[str, Any]:
        try:
            response = requests.Session().send(self.request)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            return self._extract_data(soup)
        except TimeoutError:
            # IMPORTANT: Delete cache of "waiting" page
            CACHE.delete(self.cache_key)
            if retries:
                time.sleep(60 // (retries + 1))  # backoff
                return self.run(retries - 1)
            else:
                raise
        except RuntimeError:
            # IMPORTANT: Delete cache of unexpected output
            CACHE.delete(self.cache_key)

    # TODO (gdingle): need to also handle case of queued request
    def _extract_data(self, soup: BeautifulSoup) -> Dict[str, Any]:
        title = soup.find(class_='title')
        if title and 'not present in the selected genome' in title.get_text():
            raise ValueError('Crispor on {}: {}'.format(
                self.data['seq'], title.get_text()))

        if 'retry with a sequence range shorter than 2000 bp' in soup.find(class_='contentcentral').get_text():
            raise ValueError('Crispor on {}: retry with a sequence range shorter than 2000 bp'.format(
                self.data['seq']))

        if 'This page will refresh every 10 seconds' in soup.find(class_='contentcentral').get_text():
            raise TimeoutError('Crispor on {}: Stuck in job queue. Please retry.'.format(
                self.data['seq']))

        output_table = soup.find('table', {'id': 'otTable'})
        if not output_table:
            if 'Found no possible guide sequence' in soup.get_text():
                return dict(
                    seq=self.data['seq'],
                    guide_seqs={'not found': 'not found'},
                )
            if 'Server error: could not run command' in soup.get_text():
                return dict(
                    seq=self.data['seq'],
                    guide_seqs={'server error': 'server error'},
                )
            if 'are not valid in the genome' in soup.get_text():
                return dict(
                    seq=self.data['seq'],
                    guide_seqs={
                        'invalid chromosome range': 'invalid chromosome range'
                    },
                )
            raise RuntimeError('Crispor on {}: No output rows in: {}'.format(
                self.data['seq'], soup.find('body')))

        rows = output_table.find_all(class_='guideRow')

        batch_id = soup.find('input', {'name': 'batchId'})['value']
        url = self.endpoint + '?batchId=' + batch_id
        primers_url = self.endpoint + '?batchId={}&pamId={}&pam=NGG'

        # TODO (gdingle): keeping only top three for now... what is best?
        guide_seqs = OrderedDict((t['id'], t.find_next('tt').get_text()) for t in rows[0:3])
        return dict(
            # TODO (gdingle): crispor uses seq to denote chr_loc... rename?
            # TODO (gdingle): why is off by one from input?
            # seq=soup.find(class_='title').find('a').get_text(),
            seq=self.data['seq'],
            url=url,
            batch_id=batch_id,
            guide_seqs=guide_seqs,
            primer_urls=OrderedDict((t['id'], primers_url.format(batch_id, quote(t['id']))) for t in rows),

            # TODO (gdingle): remove unneeded
            # min_freq=float(soup.find('input', {'name': 'minFreq'})['value']),
            # title=soup.title.string,
            # variant_database=soup.find('select', {'name': 'varDb'})
            # .find('option', {'selected': 'selected'})
            # .get_text(),

            fasta_url=self.endpoint + '?batchId={}&download=fasta'.format(batch_id),
            benchling_url=self.endpoint + '?batchId={}&download=benchling'.format(batch_id),
            # TODO (gdingle): other formats?
            guides_url=self.endpoint + '?batchId={}&download=guides&format=tsv'.format(batch_id),
            offtargets_url=self.endpoint + '?batchId={}&download=offtargets&format=tsv'.format(batch_id),
        )


class CrisporGuideRequestById(CrisporGuideRequest):
    """
    Given an existing Crispor batch, gets candidate guides.

    # TODO (gdingle): fix doctests
    > batch_id = 'aev3eGeG2aIZT1hdHIJ1'
    > data = CrisporGuideRequestById(batch_id).run()
    > data['batch_id'] == batch_id
    True

    > batch_id = '5JS3eHUiAeaV6eTSZ9av'
    > data = CrisporGuideRequestById(batch_id).run()
    Traceback (most recent call last):
    ...
    ValueError: Crispor on 5JS3eHUiAeaV6eTSZ9av: Query sequence, not present in the selected genome, Homo sapiens (hg19)
    """

    def __init__(self, batch_id) -> None:
        self.endpoint = 'http://crispor.tefor.net/crispor.py?batchId=' + batch_id
        self.data = {'seq': batch_id, 'pam_id': batch_id}  # hack for error messages

    def run(self, retries: int=0) -> Dict[str, Any]:
        response = requests.get(self.endpoint)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        return self._extract_data(soup)


class CrisporPrimerRequest(AbstractCrisporRequest):
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
            amp_len: int=400,
            tm: int=60,
            pam: str='NGG',
            seq: str='') -> None:

        self.pam_id = pam_id
        quoted_pam_id = quote(pam_id)  # percent encode the '+' symbol
        self.endpoint = 'http://crispor.tefor.net/crispor.py' + \
            '?ampLen={amp_len}&tm={tm}&batchId={batch_id}&pamId={quoted_pam_id}&pam={pam}'.format(**locals())
        self.seq = seq  # just for metadata
        self.request = requests.Request('GET', self.endpoint).prepare()

    def __repr__(self):
        return 'CrisporPrimerRequest({})'.format(self.endpoint)

    def run(self,
            retries: int=1) -> Dict[str, Any]:
        try:
            response = requests.Session().send(self.request)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            return self._extract_data(soup)
        except RuntimeError:
            if retries:
                return self.run(retries - 1)
            else:
                raise

    def _extract_data(self, soup: BeautifulSoup) -> Dict[str, Any]:
        # TODO (gdingle): all primers still wanted?
        # primer_tables, primer_seqs = self._extract_tables(soup)
        if 'exceptions.ValueError' in soup.get_text():
            raise RuntimeError('Crispor exceptions.ValueError')

        return dict(
            pam_id=self.pam_id,
            # TODO (gdingle): crispor uses seq to denote chr_loc
            seq=self.seq,
            url=self.endpoint,
            amplicon_length=soup
            .find('select', {'name': 'ampLen'})
            .find('option', {'selected': 'selected'})['value'],
            primer_temp=soup
            .find('select', {'name': 'tm'})
            .find('option', {'selected': 'selected'})['value'],
            # TODO (gdingle): may not need primer_seqs
            # primer_tables=primer_tables,
            # primer_seqs=primer_seqs,
            ontarget_primers=self._extract_ontarget_primers(soup)
        )

    def _extract_ontarget_primers(self, soup: BeautifulSoup) -> Dict[str, str]:
        table = soup.find(id='ontargetPcr').find_next(class_='primerTable')
        rows = (row.find_all('td') for row in table.find_all('tr'))
        return dict((row[0].get_text().split('_')[-1], row[1].get_text()) for row in rows)

    def _extract_tables(self, soup: BeautifulSoup) -> Tuple[OrderedDict, OrderedDict]:
        table_dict = OrderedDict([])  # type: OrderedDict[str, dict]
        seq_dict = OrderedDict([])  # type: OrderedDict[str, str]
        tables = soup.find_all(class_='primerTable')
        for table in tables:
            table_name = table.find_previous('h3').get_text()
            rows = table.find_all('tr')
            table_data = {}
            for row in rows:
                cells = row.find_all('td')
                if len(cells):  # not all primerTable have data
                    primer_id = cells[0].get_text()
                    primer_id = primer_id.replace('(constant primer used for all guide RNAs)', '').strip()
                    primer_seq = cells[1].get_text()
                    table_data[primer_id] = primer_seq
                    seq_dict[primer_id] = primer_seq
            table_dict[table_name] = table_data
        return table_dict, seq_dict


# TODO (gdingle): rename module if we end up using TagIn
class TagInRequest(AbstractCrisporRequest):
    """
    Given an an Ensembl Transcript or a custom sequence, gets candidate sgRNAs
    and ssDNAs from tagin.stembio.org service for HDR experiments.

    To bypass CSRF protection, we make a initial request to extract the CSRF
    token out of the form, so we can submit it with the main request.

    The sessionid is also needed, for unknown reasons.

    >>> data = TagInRequest('ENST00000330949').run()
    >>> len(data['guide_seqs']) >= 1
    True
    >>> len(data['donor_seqs']) >= 1
    True
    >>> all(s in data['guide_seqs'].values() for s in data['donor_seqs'].keys())
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
        response = requests.get(url, cookies=cookies)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        return self._extract_data(soup, url)

    def _init_session(self) -> None:
        initial = requests.get(self.endpoint)
        initial.raise_for_status()
        self.csrftoken = initial.cookies['csrftoken']

    def _get_design_id(self) -> str:
        self.data['csrfmiddlewaretoken'] = self.csrftoken
        headers = {
            'X-Requested-With': 'XMLHttpRequest',
            'Cookie': 'csrftoken={}'.format(self.csrftoken),
        }

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
        # TODO (gdingle): keeping only top three for now... what is best?
        top_guides = sorted(data_user[1], key=lambda d: d['sgRNA_score'], reverse=True)[:3]
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
