import json

from django.forms import ModelForm, widgets
from main.models import *


# TODO (gdingle): set placeholders or helptext


class PrettyJsonWidget(widgets.Textarea):

    def format_value(self, value):
        # Incredibly, there is no easier way to pretty format JSON in Django.
        return json.dumps(json.loads(value), indent=2)


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
        exclude = ['experiment', 'guide_data']
        widgets = {'guide_data': PrettyJsonWidget(attrs={'rows': 20})}


class GuideSelectionForm(ModelForm):
    class Meta:
        model = GuideSelection
        fields = '__all__'
        exclude = ['guide_design']
        widgets = {'selected_guides': PrettyJsonWidget(attrs={'rows': 20})}


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
        widgets = {'selected_primers': PrettyJsonWidget(attrs={'rows': 60})}


class PrimerPlateLayoutForm(ModelForm):
    class Meta:
        model = PrimerPlateLayout
        fields = '__all__'
        exclude = ['primer_selection']
