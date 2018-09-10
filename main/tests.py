import doctest

from django.test import TestCase

from main import conversions, validators

from main.models import *
from main.samplesheet import *
from main.samplesheet import _insert_fastqs, _new_samplesheet


def load_tests(loader, tests, ignore):
    """
    Need to add doctests manually.
    See https://stackoverflow.com/questions/2380527/django-doctests-in-views-py
    """
    tests.addTests(doctest.DocTestSuite(conversions))
    tests.addTests(doctest.DocTestSuite(validators))
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
            targets=["chr2:38377154-38377424"],
            target_seqs=["chr2:38377154-38377424"],
            guide_data=[{
                "seq": "chr2:38377154-38377424",
                "target": "chr2:38377154-38377424",
                "batch_id": "R1k4GVEcYvRcHOPSDpJk",
                "guide_seqs": {
                        "s28+": "ACGTGGTTAACCGCGGCGCT TGG",
                        "s29+": "CGTGGTTAACCGCGGCGCTT GGG",
                        "s47+": "TTGGGTCGCTGGTCCGTCGC CGG",
                        }}])

    @property
    def _guide_selection(self):
        return GuideSelection(
            guide_design=self._guide_design,
            selected_guides={
                "chr2:38377154-38377424": {
                    "s28+": "ACGTGGTTAACCGCGGCGCT TGG",
                    "s29+": "CGTGGTTAACCGCGGCGCTT GGG",
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
                    "GGTTCTCCCAGCAGCTACTG",
                    "GTTTGACGTCAGTGGGGAGT"
                ]})

    @property
    def _analysis(self):
        return Analysis(
            experiment=self._experiment,
            researcher=Researcher(),
            fastqs=[
                'A3-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz',
                'A1-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz',
                'A3-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz',
                'A1-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz',
            ],
            results_data=[{
                "success": True,
                "log_params": "",
                "report_url": "http://crispresso.pinellolab.partners.org/view_report/8VYlbE",
                "report_zip": "http://crispresso.pinellolab.partners.org/reports_data/CRISPRessoRun8VYlbE/CRISPResso_Report_8VYlbE.zip",
                "report_files": [],
                "report_stats": {"Total": 14470}}])

    def test_from_experiment(self):
        sheet = from_experiment(self._experiment)
        self.assertTrue(len(sheet))

    def test_from_guide_selection(self):
        sheet = from_guide_selection(self._guide_selection)
        self.assertTrue(len(sheet))

    def test_from_primer_selection(self):
        sheet = from_primer_selection(self._primer_selection)
        self.assertTrue(len(sheet))

    def test_from_analysis_and_primer_selection(self):
        sheet = from_analysis_and_primer_selection(self._analysis, self._primer_selection)
        self.assertTrue(len(sheet))

    def test_insert_fastqs(self):
        sheet = _new_samplesheet()
        self.assertEqual(len(sheet), 96)
        fastqs = [
            'A3-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz',
            'A1-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz',
            'A3-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz',
            'A1-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz',
        ]
        sheet = _insert_fastqs(sheet, fastqs)
        self.assertEqual(sheet.loc['A1', 'fastq_fwd'], 'A1-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz')
        self.assertEqual(sheet.loc['A1', 'fastq_rev'], 'A1-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz')
        self.assertEqual(sheet.loc['A3', 'fastq_fwd'], 'A3-BCAP31-C-sorted-180212_S3_L001_R1_001.fastq.gz')
        self.assertEqual(len(sheet), 2)
