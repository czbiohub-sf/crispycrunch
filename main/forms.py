import json

from django.contrib.postgres.forms import SimpleArrayField
from django.forms import ModelForm, widgets
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


class NewlineArrayField(SimpleArrayField):

    def __init__(self, *args, **kwargs):
        kwargs['widget'] = widgets.Textarea(attrs={'rows': 10})
        kwargs['delimiter'] = '\n'
        kwargs['max_length'] = 96  # same a standard plate size
        kwargs['min_length'] = 1
        super().__init__(*args, **kwargs)


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
        exclude = ['experiment', 'guide_data', 'donor_data']
        field_classes = {
            'targets': NewlineArrayField,
            'target_fastas': NewlineArrayField,
        }


class GuideSelectionForm(ModelForm):
    class Meta:
        model = GuideSelection
        fields = '__all__'
        exclude = ['guide_design']
        widgets = {
            'selected_guides': PrettyJsonWidget(attrs={'rows': 20}),
            'selected_guides_tagin': PrettyJsonWidget(attrs={'rows': 20}),
            'selected_donors': PrettyJsonWidget(attrs={'rows': 20}),
        }


class GuidePlateLayoutForm(ModelForm):
    class Meta:
        exclude = ['guide_selection']
        model = GuidePlateLayout
        fields = '__all__'


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


class PrimerPlateLayoutForm(ModelForm):
    class Meta:
        model = PrimerPlateLayout
        fields = '__all__'
        exclude = ['primer_selection']
