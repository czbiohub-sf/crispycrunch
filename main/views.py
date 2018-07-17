"""
Views here are mostly subclasses of the Django generic CreateView. Each view has
its own template, named in snake-case for the model to-be-created. Each view
presents a Django form based on a model.

The views are linked together in a linear sequence that reflects the
dependencies of the data. Successfully submitting one form redirects the user to
the next form. Each model has a foreign key into the preceding model. One
exception is PrimerDesign which depends on GuideSelection, two steps back in the
sequence.
"""

from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django.views import View
from django.views.generic import DetailView
from django.views.generic.edit import CreateView
# TODO (gdingle): why is lint not working here?
from openpyxl import Workbook, writer  # noqa

import crisporclient

from main.forms import *
from main.models import *
from main.platelayout import *

# TODO (gdingle): create useful index


def index(request):
    return HttpResponse("Hello, world. You're at the main index.")

#
# BEGIN EXPERIMENT CREATION VIEWS
#


class CreatePlusView(CreateView):
    """
    Simplifies adding foreign keys and other pre-determined data to a ModelForm
    before saving it.

    See https://github.com/django/django/blob/master/django/views/generic/edit.py.
    """

    def form_valid(self, form: ModelForm) -> HttpResponse:
        obj = form.save(commit=False)
        obj = self.plus(obj)
        obj.save()
        self.object = obj
        return HttpResponseRedirect(self.get_success_url())

    def plus(self, obj: models.Model) -> models.Model:
        """Add attributes to a model before saving it."""


class ExperimentView(CreateView):
    template_name = 'experiment.html'
    form_class = ExperimentForm
    success_url = '/main/experiment/{id}/guide-design/'


class GuideDesignView(CreatePlusView):
    template_name = 'guide-design.html'
    form_class = GuideDesignForm
    success_url = '/main/guide-design/{id}/guide-selection/'

    def plus(self, obj):
        obj.experiment = Experiment.objects.get(id=self.kwargs['id'])
        obj.guide_data = crisporclient.CrisporGuideRequest(
            obj.targets,
            # TODO (gdingle): does experiment name get us anything useful?
            name=obj.experiment.name,
            org=obj.genome,
            pam=obj.pam).run()
        return obj


class GuideSelectionView(CreatePlusView):
    template_name = 'guide-selection.html'
    form_class = GuideSelectionForm
    success_url = '/main/guide-selection/{id}/guide-plate-layout/'

    def get_initial(self):
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        return {
            'selected_guides': guide_design.guide_data['guide_seqs']
        }

    def plus(self, obj):
        obj.guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        return obj


class GuidePlateLayoutView(CreatePlusView):
    template_name = 'guide-plate-layout.html'
    form_class = GuidePlateLayoutForm
    success_url = '/main/guide-selection/{guide_selection_id}/primer-design/'

    def plus(self, obj):
        obj.guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        return obj

    def get_context_data(self, **kwargs):
        guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        kwargs['plate_layout'] = Plate96Layout(guide_selection.selected_guides)
        return super().get_context_data(**kwargs)


class PrimerDesignView(CreatePlusView):
    template_name = 'primer-design.html'
    form_class = PrimerDesignForm
    success_url = '/main/primer-design/{id}/primer-selection/'

    def get_initial(self):
        guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        return {
            'pam': guide_selection.guide_design.pam,
            # TODO (gdingle): this should be pam_ids... all selected guides must
            # have primers
            'pam_id': guide_selection.selected_guides.popitem()[0],
        }

    def plus(self, obj):
        obj.guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        batch_id = obj.guide_selection.guide_design.guide_data['batch_id']
        obj.primer_data = crisporclient.CrisporPrimerRequest(
            batch_id=batch_id,
            amp_len=obj.maximum_amplicon_length,
            tm=obj.primer_temp,
            pam=obj.pam,
            pam_id=obj.pam_id).run()
        return obj


class PrimerSelectionView(CreatePlusView):
    template_name = 'primer-selection.html'
    form_class = PrimerSelectionForm
    success_url = '/main/primer-selection/{id}/primer-plate-layout/'

    def get_initial(self):
        return {
            'selected_primers': PrimerDesign.objects.get(id=self.kwargs['id']).primer_data['primer_seqs']
        }

    def plus(self, obj):
        obj.primer_design = PrimerDesign.objects.get(id=self.kwargs['id'])
        return obj

    def get_context_data(self, **kwargs):
        kwargs['crispor_url'] = PrimerDesign.objects.get(id=self.kwargs['id']).primer_data['url']
        return super().get_context_data(**kwargs)


class PrimerPlateLayoutView(CreatePlusView):
    template_name = 'primer-plate-layout.html'
    form_class = PrimerPlateLayoutForm
    success_url = '/main/primer-plate-layout/{id}/experiment-summary/'

    def plus(self, obj):
        obj.primer_selection = PrimerSelection.objects.get(id=self.kwargs['id'])
        return obj

    def get_context_data(self, **kwargs):
        primer_selection = PrimerSelection.objects.get(id=self.kwargs['id'])
        kwargs['plate_layout'] = Plate384Layout(primer_selection.selected_primers)
        return super().get_context_data(**kwargs)

#
# END EXPERIMENT CREATION VIEWS
#


class ExperimentSummaryView(View):
    template_name = 'experiment-summary.html'

    def get(self, request, *args, **kwargs):
        # Fetch objects associated with new experiment by traversing the graph.
        primer_plate_layout = PrimerPlateLayout.objects.get(id=kwargs['id'])
        primer_selection = primer_plate_layout.primer_selection
        primer_design = primer_selection.primer_design
        guide_selection = primer_design.guide_selection
        guide_plate_layout = GuidePlateLayout.objects.filter(guide_selection=guide_selection)[0]
        guide_design = guide_selection.guide_design
        experiment = guide_design.experiment

        return render(request, self.template_name, locals())


class OrderFormView(DetailView):
    """
    Produces a downloadable Excel order form for IDT. The model must have a
    plate layout.
    """

    model: models.Model = None

    def _create_excel_file(self, plate_layout: AbstractPlateLayout, title: str):
        wb = Workbook()
        ws = wb.active

        ws.title = title[0:31]  # Excel limits to 30 chars
        ws['A1'] = 'Well Position'
        ws['B1'] = 'Name'
        ws['C1'] = 'Sequence'

        # TODO (gdingle): should print all well positions or just until end of contents?
        for i, pos in enumerate(plate_layout.well_positions):
            index = str(i + 2)
            ws['A' + index] = pos
            ws['B' + index] = plate_layout.well_names[pos]
            ws['C' + index] = plate_layout.well_seqs[pos]

        return writer.excel.save_virtual_workbook(wb)

    def get(self, request, *args, **kwargs):
        instance = self.model.objects.get(id=kwargs['id'])
        # TODO (gdingle): friendlier title?
        title = request.path.replace('/', ' ').replace('main ', '')
        excel_file = self._create_excel_file(instance.layout, title)

        response = HttpResponse(
            excel_file, content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename="{}.xlsx"'.format(title)

        return response


class GuideOrderFormView(OrderFormView):

    model = GuidePlateLayout


class PrimerOrderFormView(OrderFormView):

    model = PrimerPlateLayout
