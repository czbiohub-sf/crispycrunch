import doctest

from django.test import TestCase

from utils import *

from main.models import *
from main.samplesheet import *
from main.samplesheet import _from_analysis, _insert_fastqs, _new_samplesheet


def load_tests(loader, tests, ignore):
    """
    Need to add doctests manually.
    See https://stackoverflow.com/questions/2380527/django-doctests-in-views-py
    """
    modules = [
        conversions,
        validators,
        chrloc,
        hdr,
    ]
    for module in modules:
        tests.addTests(doctest.DocTestSuite(module))
    return tests


class SampleSheetTestCase(TestCase):
    """
    # TODO (gdingle): write more specific tests
    """

    @property
    def _experiment(self):
        return Experiment(researcher=Researcher())

    @property
    def _guide_design(self):
        return GuideDesign(
            experiment=self._experiment,
            targets_raw=["ENST00000621663"],
            targets=["chr2:38377154-38377424"],
            target_seqs=['ATGACGTGGTTAACCGCGGCGCTTGGG'],
            target_genes=["ATL2"],
            guide_data=[{
                "seq": "chr2:38377154-38377424",
                "target": "chr2:38377154-38377424",
                "batch_id": "R1k4GVEcYvRcHOPSDpJk",
                "scores": {
                        "s28+": ['50', '50', '50'],
                        "s29+": ['50', '50', '50'],
                        "s47+": ['50', '50', '50'],
                        },
                "guide_seqs": {
                    "s28+": "ACGTGGTTAACCGCGGCGCT TGG",
                    "s29+": "CGTGGTTAACCGCGGCGCTT GGG",
                    "s47+": "TTGGGTCGCTGGTCCGTCGC CGG",
                },
                'url': 'http://crispor.tefor.net/crispor.py?batchId=gNw3bkdHkk9DEjitEjGd',
                'primer_urls': {
                    's28+': 'http://crispor.tefor.net/crispor.py?batchId=gNw3bkdHkk9DEjitEjGd&pamId=s16+&pam=NGG',
                    's29+': 'http://crispor.tefor.net/crispor.py?batchId=gNw3bkdHkk9DEjitEjGd&pamId=s17+&pam=NGG',
                    's47+': 'http://crispor.tefor.net/crispor.py?batchId=gNw3bkdHkk9DEjitEjGd&pamId=s18-&pam=NGG',
                },
            }],
        )

    @property
    def _guide_selection(self):
        return GuideSelection(
            guide_design=self._guide_design,
            selected_guides={
                "chr2:38377154-38377424": {
                    "s28+": "ACGTGGTTAACCGCGGCGCT TGG",
                    # Commented out to demonstrate selection of subset
                    # "s29+": "CGTGGTTAACCGCGGCGCTT GGG",
                    "s47+": "TTGGGTCGCTGGTCCGTCGC CGG",
                }
            })

    @property
    def _primer_design(self):
        return PrimerDesign(guide_selection=self._guide_selection)

    @property
    def _primer_selection(self):
        return PrimerSelection(
            primer_design=self._primer_design,
            selected_primers={
                "chr2:38377154-38377424 s28+": [
                    # TODO (gdingle): these primers are not realistic
                    ("GGTTCTCCCAGCAGCTACTG",
                     "GGTTCTCCCAGCAGCTACTGACGTGGTTAACCGCGGCGCTCGTGGTTAACCGCGGCGCTTTTGGGTCGCTGGTCCGTCGC"),
                    ("GTTTGACGTCAGTGGGGAGT",
                     "GTTTGACGTCAGTGGGGAGTACGTGGTTAACCGCGGCGCTCGTGGTTAACCGCGGCGCTTTTGGGTCGCTGGTCCGTCGC"),
                ]})

    @property
    def _analysis(self):
        return Analysis(
            experiment=self._experiment,
            researcher=Researcher(),
            fastq_data=[
                ('A1-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz',
                 'A1-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz',),
            ],
            results_data=[{
                "success": True,
                "log_params": "",
                "report_url": "http://crispresso.pinellolab.partners.org/view_report/8VYlbE",
                "report_zip": "http://crispresso.pinellolab.partners.org/reports_data/CRISPRessoRun8VYlbE/CRISPResso_Report_8VYlbE.zip",
                "report_files": [],
                "report_stats": {"Total": 14470}}])

    def test_guide_design_to_df(self):
        df = self._guide_design.to_df()
        num_guides = sum(len(g['guide_seqs']) for g in self._guide_design.guide_data)
        self.assertEqual(len(df), num_guides)

    def test_guide_selection_to_df(self):
        df = self._guide_selection.to_df()
        num_guides = sum(len(list(g)) for g in self._guide_selection.selected_guides.values())
        self.assertEqual(len(df), num_guides)

    def test_primer_selection_to_df(self):
        df = self._primer_selection.to_df()
        self.assertEqual(len(df), len(self._primer_selection.selected_primers))

    def test_from_experiment(self):
        sheet = from_experiment(self._experiment)
        self.assertTrue(len(sheet))

    def test_from_guide_selection(self):
        sheet = from_guide_selection(self._guide_selection)
        self.assertTrue(len(sheet))

    def test_from_primer_selection(self):
        sheet = from_primer_selection(self._primer_selection)
        self.assertTrue(len(sheet))

    def test_from_analysis(self):
        sheet = from_primer_selection(self._primer_selection)
        sheet = _from_analysis(self._analysis, sheet)
        self.assertTrue(len(sheet))

    # TODO (gdingle): _insert_fastqs is deprecated
    def test_insert_fastqs(self):
        sheet = _new_samplesheet()
        self.assertGreaterEqual(len(sheet), 96)
        fastqs = [
            'A3-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz',
            'A1-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz',
            'A3-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz',
            'A1-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz',
        ]
        sheet = _insert_fastqs(sheet, fastqs)
        self.assertEqual(sheet.loc['A1', 'fastq_fwd'],
                         'A1-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz')
        self.assertEqual(sheet.loc['A1', 'fastq_rev'],
                         'A1-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz')
        self.assertEqual(sheet.loc['A3', 'fastq_fwd'],
                         'A3-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz')
        self.assertEqual(len(sheet), 2)
