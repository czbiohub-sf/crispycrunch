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
        samplesheet = from_experiment(experiment)
        self.assertTrue(len(samplesheet))
