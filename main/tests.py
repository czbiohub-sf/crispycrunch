import doctest

from django.test import TestCase

from main import conversions, platelayout, validators

from main.models import *
from main.samplesheet import *
from main.samplesheet import _insert_fastqs, _new_samplesheet


def load_tests(loader, tests, ignore):
    """
    Need to add doctests manually.
    See https://stackoverflow.com/questions/2380527/django-doctests-in-views-py
    """
    tests.addTests(doctest.DocTestSuite(platelayout))
    tests.addTests(doctest.DocTestSuite(conversions))
    tests.addTests(doctest.DocTestSuite(validators))
    return tests


class SampleSheetTestCase(TestCase):
    """
    NOTE: These tests were run manually during development of samplesheet.py.
    They are likely to fail now because the data and data model has changed.
    # TODO (gdingle): make tests robust enough to run automatically
    """

    # TODO (gdingle): slim down the fixtures data
    fixtures = ['main/tests.json']

    # def test_from_experiment(self):
    #     experiment = Experiment.objects.get(id=1)
    #     sheet = from_experiment(experiment)
    #     self.assertTrue(len(sheet))

    # def test_from_guide_selection(self):
    #     guide_selection = GuideSelection.objects.get(id=1)
    #     sheet = from_guide_selection(guide_selection)
    #     self.assertTrue(len(sheet))

    # def test_from_primer_selection(self):
    #     primer_selection = PrimerSelection.objects.get(id=1)
    #     sheet = from_primer_selection(primer_selection)
    #     self.assertTrue(len(sheet))

    # def test_from_analysis(self):
    #     analysis = Analysis.objects.get(id=1)
    #     sheet = from_analysis(analysis)
    #     self.assertTrue(len(sheet))

    def test_insert_fastqs(self):
        # TODO (gdingle):
        """
        """
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
