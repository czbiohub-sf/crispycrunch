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


class ExperimentForm(ModelForm):
    class Meta:
        model = Experiment
        fields = '__all__'
        exclude = ['owner', ]


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


class GuideDesignForm2(ModelForm):
    """
    Same as above but excludes hdr fields.
    # TODO (gdingle): refactor both classes to "include" instead of "exclude"
    """

    class Meta:
        model = GuideDesign
        fields = '__all__'
        exclude = GuideDesignForm.Meta.exclude + [
            'hdr_tag',
            'hdr_start_codon_tag_seq',
            'hdr_stop_codon_tag_seq',
            'hdr_homology_arm_length',
        ]
        field_classes = GuideDesignForm.Meta.field_classes

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        targets_raw = self.fields['targets_raw']
        targets_raw.initial = [
            'chr2:136114360-136114419',
            'chr2:136115613-136115672',
            'chr2:136116738-136116797',
            'chr2:136117544-136117603',
        ]


class GuideSelectionForm(ModelForm):
    class Meta:
        model = GuideSelection
        fields = '__all__'
        exclude = ['owner', 'guide_design']
        widgets = {
            'selected_guides': PrettyJsonWidget(
                attrs={'rows': 20, 'spellcheck': "false"}),
        }


class PrimerDesignForm(ModelForm):

    class Meta:
        model = PrimerDesign
        fields = '__all__'
        exclude = ['owner', 'guide_selection', 'primer_data']

    def __init__(self, *args, **kwargs):
        self._is_hdr = GuideSelection.objects.get(id=kwargs['id']).guide_design.is_hdr
        del kwargs['id']
        super().__init__(*args, **kwargs)
        max_amplicon_length = self.fields['max_amplicon_length']
        max_amplicon_length.initial = self._initial
        max_amplicon_length.help_text = self._help_text
        max_amplicon_length.disabled = self._disabled

        # TODO (gdingle): make crispor.py more flexible
        primer_temp = self.fields['primer_temp']
        primer_temp.disabled = self._disabled

        # Only enforced for non-hdr
        max_amplicon_length.widget.attrs['step'] = 100

    @property
    def _help_text(self):
        if self._is_hdr:
            # TODO (gdingle): make amplicon size for hdr variable in crispor.py
            # and then variable here.
            return """Primers will flank each side of insertion point by at
            least 105bp. The resulting amplicon length will be between 250-310bp,
            not including the HDR tag sequence. Shorter amplicons are preferred
            in returned primers."""
        else:
            return """Primers will be centered on the guide for optimal
                sequencing of NHEJ. The resulting amplicon length will be
                between: [max_length - 150, max_length] if max_length > 250, or
                [max_length - 50, max_length] if max_length <= 250."""

    @property
    def _initial(self):
        if self._is_hdr:
            return 310
        else:
            return self.fields['max_amplicon_length'].initial

    @property
    def _disabled(self):
        if self._is_hdr:
            return True
        else:
            return False


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
        fields = (
            'username',
            'email',
            'first_name',
            'last_name',
            'password1',
            'password2',
        )

    def save(self, commit=True):
        user = super().save(commit=False)
        user.email = self.cleaned_data['email']
        if commit:
            user.save()
        return user
