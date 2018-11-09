import json

from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django.contrib.postgres.forms import SimpleArrayField
from django.forms import EmailField, FileField, Form, ModelForm, widgets

from main.models import *

# TODO (gdingle): consider YAML instead for easier editing
# See https://pyyaml.org/wiki/PyYAMLDocumentation


class PrettyJsonWidget(widgets.Textarea):

    def format_value(self, value):
        # Incredibly, there is no easier way to pretty format JSON in Django.
        return json.dumps(json.loads(value), indent=2)


# TODO (gdingle): get this working or use YAML
# class FriendlyJSONField(JSONField):

#     def clean(self, value):
#         """Strip trailing commas. WARNING: This may mangle valid JSON."""
#         value = re.sub(",\s+}", "}", value)
#         value = re.sub(",\s+\]", "]", value)
#         assert False, value
#         return self.clean(value)


class NewlineArrayField(SimpleArrayField):

    def __init__(self, *args, **kwargs):
        kwargs['widget'] = widgets.Textarea(attrs={'rows': 10})
        kwargs['delimiter'] = '\n'
        kwargs['max_length'] = 96  # same a standard plate size
        kwargs['min_length'] = 1
        super().__init__(*args, **kwargs)

    def clean(self, value):
        """Strip empty lines"""
        return super().clean(value.strip())


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
        exclude = ['owner',
                   'experiment',
                   'guide_data',
                   'target_locs',
                   'target_seqs',
                   'target_genes',
                   'target_tags',
                   ]
        field_classes = {
            'targets_raw': NewlineArrayField,
            'target_fastas': NewlineArrayField,
        }


class GuideSelectionForm(ModelForm):
    class Meta:
        model = GuideSelection
        fields = '__all__'
        exclude = ['owner', 'guide_design']
        widgets = {
            'selected_guides': PrettyJsonWidget(
                attrs={'rows': 40, 'spellcheck': "false"}),
        }


class PrimerDesignForm(ModelForm):
    class Meta:
        model = PrimerDesign
        fields = '__all__'
        exclude = ['owner', 'guide_selection', 'primer_data']


class PrimerSelectionForm(ModelForm):

    class Meta:
        model = PrimerSelection
        fields = '__all__'
        exclude = ['owner', 'primer_design']
        widgets = {'selected_primers': PrettyJsonWidget(
            attrs={'rows': 20, 'spellcheck': "false"})}


class AnalysisForm(ModelForm):
    class Meta:
        model = Analysis
        fields = '__all__'
        exclude = ['owner', 'fastq_data', 'results_data']

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


class CustomUserCreationForm(UserCreationForm):
    """
    All to put email in the form :/
    """
    email = EmailField(label='Email address', required=True,
                       help_text='Required.')

    class Meta:
        model = User
        fields = ('username', 'email', 'password1', 'password2')

    def save(self, commit=True):
        user = super().save(commit=False)
        user.email = self.cleaned_data['email']
        if commit:
            user.save()
        return user
