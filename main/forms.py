import json

from django.contrib.postgres.forms import SimpleArrayField
from django.forms import FileField, Form, ModelForm, widgets
from main.models import *


# TODO (gdingle): set placeholders or helptext


# TODO (gdingle): consider YAML instead for easier editing
# See https://pyyaml.org/wiki/PyYAMLDocumentation
class PrettyJsonWidget(widgets.Textarea):

    def format_value(self, value):
        # Incredibly, there is no easier way to pretty format JSON in Django.
        return json.dumps(json.loads(value), indent=2)
    # TODO (gdingle): can we strip trailing commas here as well?
    # see https://stackoverflow.com/questions/23705304/can-json-loads-ignore-trailing-commas


# TODO (gdingle): ArrayField is not showing full error messages, only first part
# "Item 1 in the array did not validate:". See:
# https://docs.djangoproject.com/en/2.1/_modules/django/contrib/postgres/forms/array/
class NewlineArrayField(SimpleArrayField):

    def __init__(self, *args, **kwargs):
        kwargs['widget'] = widgets.Textarea(attrs={'rows': 10})
        kwargs['delimiter'] = '\n'
        kwargs['max_length'] = 96  # same a standard plate size
        kwargs['min_length'] = 1
        super().__init__(*args, **kwargs)

    # TODO (gdingle): how to strip empty lines?
    # see https://docs.djangoproject.com/en/2.0/_modules/django/contrib/postgres/forms/array/#SimpleArrayField
    # def clean(self, value):
    #     super().__init__(value.strip())


class ResearcherForm(ModelForm):
    class Meta:
        model = Researcher
        fields = '__all__'


class ExperimentForm(ModelForm):
    class Meta:
        model = Experiment
        fields = '__all__'


class GuideDesignForm(ModelForm):

    class Meta:
        model = GuideDesign
        fields = '__all__'
        exclude = ['experiment', 'guide_data', 'target_seqs', 'hdr_seq', 'targets']
        field_classes = {
            'targets_raw': NewlineArrayField,
            'target_fastas': NewlineArrayField,
        }


class GuideSelectionForm(ModelForm):
    class Meta:
        model = GuideSelection
        fields = '__all__'
        exclude = ['guide_design']
        widgets = {
            'selected_guides': PrettyJsonWidget(attrs={'rows': 40}),
            'selected_guides_tagin': PrettyJsonWidget(attrs={'rows': 40}),
        }
        labels = {
            "selected_guides_tagin": 'Selected guides',
        }


class PrimerDesignForm(ModelForm):
    class Meta:
        model = PrimerDesign
        fields = '__all__'
        exclude = ['guide_selection', 'primer_data']


class PrimerSelectionForm(ModelForm):

    class Meta:
        model = PrimerSelection
        fields = '__all__'
        exclude = ['primer_design']
        widgets = {'selected_primers': PrettyJsonWidget(attrs={'rows': 20})}


class AnalysisForm(ModelForm):
    class Meta:
        model = Analysis
        fields = '__all__'
        exclude = ['fastq_data', 'results_data']

    # TODO (gdingle): fetch only for display experiments that have status of "ready"
    def clean_experiment(self):
        experiment = self.cleaned_data['experiment']
        if experiment.is_custom_analysis:
            return experiment
        primer_selection = PrimerSelection.objects.filter(
            primer_design__guide_selection__guide_design__experiment=experiment)
        if not primer_selection:
            raise ValidationError('Experiment "{}" is not ready for analysis. Did you complete all setup steps?'.format(
                experiment.name))
        return experiment


class CustomAnalysisForm(Form):
    file = FileField(help_text='Excel or CSV file', label='')

    def clean_file(self):
        valid = (
            'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
            'application/vnd.ms-excel',
            'text/csv',
        )
        file = self.cleaned_data['file']
        if file.content_type not in valid:
            raise ValidationError(
                'Cannot handle file format: "{}". Must be one of: {}'.format(
                    file.content_type, valid))
        return file
