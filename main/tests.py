# from django.test import TestCase
# See https://stackoverflow.com/questions/2380527/django-doctests-in-views-py

import doctest
from main import samplesheet, platelayout, conversions, validators


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(samplesheet))
    tests.addTests(doctest.DocTestSuite(platelayout))
    tests.addTests(doctest.DocTestSuite(conversions))
    tests.addTests(doctest.DocTestSuite(validators))
    return tests
