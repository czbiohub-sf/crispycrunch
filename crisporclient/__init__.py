"""
Crispor Client.

This module calls out to crispor.tefor.net endpoints. Crispor does not have
an official public web API. You've been warned!
"""
import requests

import json
import time

from abc import abstractmethod
from bs4 import BeautifulSoup
from collections import OrderedDict
from typing import Any, Dict, Tuple
from urllib.parse import quote


class AbstractCrisporRequest:

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__, self.endpoint, self.__dict__.get('data'))

    def __str__(self):
        return self.__repr__()

    @abstractmethod
    def run(self) -> Dict[str, Any]:
        """Requests self.endpoint and extracts relevant data from the HTML response"""


class CrisporGuideRequest(AbstractCrisporRequest):
    """
    Given a sequence or chromosome location, gets candidate Crispr guides from
    Crispor service.

    >>> seq = 'chr1:11,130,540-11,130,751'
    >>> data = CrisporGuideRequest(name='test-crispr-guides', seq=seq).run()
    >>> len(data['primer_urls']) > 3
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

    def run(self, retries: int=1) -> Dict[str, Any]:
        try:
            response = requests.post(self.endpoint, data=self.data)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, 'html.parser')
            return self._extract_data(soup)
        except TimeoutError:
            if retries:
                time.sleep(8)  # determined by experience
                return self.run(retries - 1)
            else:
                raise

    # TODO (gdingle): need to also handle case of queued request
    def _extract_data(self, soup: BeautifulSoup) -> Dict[str, Any]:
        title = soup.find(class_='title')
        if title and 'not present in the selected genome' in title.get_text():
            raise ValueError('Crispor: ' + title.get_text())

        if 'retry with a sequence range shorter than 2000 bp' in soup.find(class_='contentcentral').get_text():
            raise ValueError('Crispor: retry with a sequence range shorter than 2000 bp')

        if 'This page will refresh every 10 seconds' in soup.find(class_='contentcentral').get_text():
            raise TimeoutError('Stuck in Crispor job queue. Please retry.')

        output_table = soup.find('table', {'id': 'otTable'})
        if not output_table:
            raise RuntimeError('No Crispor output rows in: {}'.format(soup.find('body')))
        rows = output_table.find_all(class_='guideRow')

        batch_id = soup.find('input', {'name': 'batchId'})['value']
        url = self.endpoint + '?batchId=' + batch_id
        primers_url = self.endpoint + '?batchId={}&pamId={}&pam=NGG'

        # TODO (gdingle): keeping only top three for now... what is best?
        guide_seqs = OrderedDict((t['id'], t.find_next('tt').get_text()) for t in rows[0:3])
        return dict(
            seq=soup.find(class_='title').find('a').get_text(), # TODO (gdingle): why is off by one from input?
            url=url,
            batch_id=batch_id,
            title=soup.title.string,
            variant_database=soup.find('select', {'name': 'varDb'})
            .find('option', {'selected': 'selected'})
            .get_text(),
            min_freq=float(soup.find('input', {'name': 'minFreq'})['value']),
            guide_seqs=guide_seqs,
            primer_urls=OrderedDict((t['id'], primers_url.format(batch_id, quote(t['id']))) for t in rows),
            fasta_url=self.endpoint + '?batchId={}&download=fasta'.format(batch_id),
            benchling_url=self.endpoint + '?batchId={}&download=benchling'.format(batch_id),
            # TODO (gdingle): other formats?
            guides_url=self.endpoint + '?batchId={}&download=guides&format=tsv'.format(batch_id),
            offtargets_url=self.endpoint + '?batchId={}&download=offtargets&format=tsv'.format(batch_id),
        )


class CrisporGuideRequestById(CrisporGuideRequest):
    """
    Given an existing Crispor batch, gets candidate guides.

    >>> batch_id = '9cJNEsbfWiSKa8wlaJMZ'
    >>> data = CrisporGuideRequestById(batch_id).run()
    >>> data['batch_id'] == batch_id
    True

    >>> batch_id = '5JS3eHUiAeaV6eTSZ9av'
    >>> data = CrisporGuideRequestById(batch_id).run()
    Traceback (most recent call last):
    ...
    ValueError: Crispor: Query sequence, not present in the selected genome, Homo sapiens (hg19)
    """

    def __init__(self, batch_id) -> None:
        self.endpoint = 'http://crispor.tefor.net/crispor.py?batchId=' + batch_id

    def run(self, retries: int = 0) -> Dict[str, Any]:
        response = requests.get(self.endpoint)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        return self._extract_data(soup)


class CrisporPrimerRequest(AbstractCrisporRequest):
    """
    Gets primers for a guide generated by Crispor.

    NOTE: Crispor uses Primer3 under the covers.

    >>> data = CrisporPrimerRequest('9cJNEsbfWiSKa8wlaJMZ', 's185+').run()
    >>> len(data['ontarget_primers']) == 2
    True
    """

    def __init__(
            self,
            batch_id: str,
            pam_id: str,
            amp_len: int = 400,
            tm: int = 60,
            pam: str = 'NGG',
            seq: str ='') -> None:

        pam_id = quote(pam_id)  # percent encode the '+' symbol
        self.endpoint = 'http://crispor.tefor.net/crispor.py' + \
            '?ampLen={amp_len}&tm={tm}&batchId={batch_id}&pamId={pam_id}&pam={pam}'.format(**locals())
        self.seq = seq  # just for metadata

    def __repr__(self):
        return 'CrisporPrimerRequest({})'.format(self.endpoint)

    def run(self,
            retries: int=1) -> Dict[str, Any]:
        try:
            response = requests.get(self.endpoint)
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
        donor_seqs = OrderedDict(
            (row[0].get_text(), row[1].get_text())
            for row in rows
            if len(row) and row[0] and row[1] and row[1].get_text().strip() != 'sgRNA too far from stop codon'
            and row[0].get_text() in guide_seqs.values()
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
