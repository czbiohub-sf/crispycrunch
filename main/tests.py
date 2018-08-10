import doctest

from django.test import TestCase

from main import conversions, platelayout, validators

from main.models import *
from main.samplesheet import *


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

    # TODO (gdingle): slim down the fixtures data
    fixtures = ['main/tests.json']

    def test_from_experiment(self):
        # TODO (gdingle): pick a more interesting experiment
        experiment = Experiment.objects.get(name='testsum3')
        sheet = from_experiment(experiment)
        # TODO (gdingle): some good assertions
        self.assertTrue(len(sheet))

    def test_from_guide_selection(self):
        # TODO (gdingle): pick more intersting
        guide_selection = GuideSelection.objects.get(id=45)
        sheet = from_guide_selection(guide_selection)
        self.assertTrue(len(sheet))

    def test_from_primer_selection(self):
        primer_selection = PrimerSelection.objects.get(id=18)
        sheet = from_primer_selection(primer_selection)
        # print(sheet.head().loc[:, 'target_loc':'primer_melt_temp'])
        self.assertTrue(len(sheet))

    def test_from_analysis(self):
        # TODO (gdingle): something more intereting
        analysis = Analysis.objects.get(id=24)
        sheet = from_analysis(analysis)
        self.assertTrue(len(sheet))
        print(sheet.head().loc[:,'target_seq':])
