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
import logging
import os
import time

from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from typing import no_type_check

from django.http import Http404
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django.views import View
from django.views.generic import DetailView, ListView
from django.views.generic.edit import CreateView
from openpyxl import Workbook, writer  # type: ignore

import webscraperequest

from crispresso.s3 import download_fastqs
from main import conversions
from main import samplesheet
from main.forms import *
from main.models import *
from main.validators import is_ensemble_transcript

# TODO (gdingle): move somewhere better
CRISPRESSO_ROOT_URL = 'http://crispresso:5000/'
CRISPRESSO_PUBLIC_ROOT_URL = 'http://0.0.0.0:5000/'

logger = logging.getLogger(__name__)


class IndexView(ListView):
    model = Experiment
    template_name = 'index.html'

    # TODO (gdingle): fill in
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context['analyses'] = Analysis.objects.all()
        return context


def index(request):
    # TODO (gdingle): create useful index
    return HttpResponse("Hello, world. You're at the main index.")


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
    success_url = '/main/guide-design/{id}/progress/'

    def _normalize_targets(self, targets):
        # TODO (gdingle): handle mix of target types
        if all(is_seq(t) for t in targets):
            return targets

        if not all(is_gene(t) for t in targets):
            return targets

        with ThreadPoolExecutor() as pool:
            normalized = list(pool.map(
                # TODO (gdingle): this still needs some work to get best region of gene
                conversions.gene_to_chr_loc,
                targets,
            ))
        # TODO: normalize seqs also
        assert all(is_chr(t) for t in normalized)
        return normalized

    def _get_target_seqs(self, targets, genome):
        if all(is_seq(t) for t in targets):
            # TODO (gdingle): return chr locations for sequences?
            return targets

        with ThreadPoolExecutor() as pool:
            seqs = list(pool.map(
                functools.partial(conversions.chr_loc_to_seq, genome=genome),
                targets,
            ))
        return seqs

    def plus(self, obj):
        """
        If an HDR tag-in experiment, get donor DNA then get guides. If not, just
        get guides.
        """
        obj.experiment = Experiment.objects.get(id=self.kwargs['id'])

        obj.targets = self._normalize_targets(obj.targets)

        obj.target_seqs = self._get_target_seqs(obj.targets, obj.genome)

        # TODO (gdingle): ignore HDR for now
        # if obj.hdr_seq:
        #     # TODO (gdingle): put in form validation somehow
        #     assert all(is_ensemble_transcript(t) and len(t) <= 600 for t in obj.targets), 'Bad input for TagIn'

        batch = webscraperequest.CrisporGuideBatchWebRequest(obj)
        largs = [[target_seq, obj.experiment.name, obj.genome, obj.pam, target]
                 for target_seq, target in zip(obj.target_seqs, obj.targets)]
        batch.start(largs, [-1])

        return obj


class GuideDesignProgressView(View):

    template_name = 'guide-design-progress.html'
    success_url = '/main/guide-design/{id}/guide-selection/'

    def get(self, request, **kwargs):
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        batch_status = webscraperequest.CrisporGuideBatchWebRequest(guide_design).get_batch_status()

        if not batch_status.is_successful:
            return render(request, self.template_name, locals())
        else:
            return HttpResponseRedirect(
                self.success_url.format(id=self.kwargs['id']))


class GuideSelectionView(CreatePlusView):
    template_name = 'guide-selection.html'
    form_class = GuideSelectionForm
    success_url = '/main/guide-selection/{id}/primer-design/'

    @staticmethod
    def _slice(odict, top=3):
        """
        OrderedDict does not support slicing :(
        See https://stackoverflow.com/a/30975520.
        """
        from itertools import islice
        from collections import OrderedDict
        return OrderedDict(islice(odict.items(), top))

    def get_initial(self):
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        # TODO (gdingle): maybe make this a guide design option to "fit to plate"
        wells_per_target = max(1, 96 // len(guide_design.targets))
        return {
            'selected_guides': dict(
                (g['target'], self._slice(g['guide_seqs'], wells_per_target))
                for g in guide_design.guide_data),
            'selected_donors': dict((g['metadata']['chr_loc'], g['donor_seqs'])
                                    for g in guide_design.donor_data),
            'selected_guides_tagin': dict(
                (g['metadata']['chr_loc'], self._slice(g['guide_seqs'], wells_per_target))
                for g in guide_design.donor_data),
        }

    def plus(self, obj):
        obj.guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        return obj

    def get_context_data(self, **kwargs):
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        kwargs['crispor_url'] = [
            gd['url']
            for gd in guide_design.guide_data
            if gd.get('url')][0]
        if guide_design.donor_data:
            kwargs['tagin_url'] = guide_design.donor_data[0]['url']
        return super().get_context_data(**kwargs)


class PrimerDesignView(CreatePlusView):
    template_name = 'primer-design.html'
    form_class = PrimerDesignForm
    success_url = '/main/primer-design/{id}/progress/'

    def plus(self, obj):
        guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        obj.guide_selection = guide_selection

        sheet = samplesheet.from_guide_selection(guide_selection)
        batch = webscraperequest.CrisporPrimerBatchWebRequest(obj)
        largs = [[row['_crispor_batch_id'],
                  row['_crispor_pam_id'],
                  obj.max_amplicon_length,
                  obj.primer_temp,
                  guide_selection.guide_design.pam,
                  row['target_loc']]
                 for row in sheet.to_records()]
        batch.start(largs, [-1])

        # TODO (gdingle): run crispr-primer if HDR experiment
        # https://github.com/chanzuckerberg/crispr-primer
        return obj


class PrimerDesignProgressView(View):
    template_name = 'primer-design-progress.html'
    success_url = '/main/primer-design/{id}/primer-selection/'

    def get(self, request, **kwargs):
        primer_design = PrimerDesign.objects.get(id=kwargs['id'])
        batch_status = webscraperequest.CrisporPrimerBatchWebRequest(primer_design).get_batch_status()

        if not batch_status.is_successful:
            return render(request, self.template_name, locals())
        else:
            return HttpResponseRedirect(
                self.success_url.format(id=self.kwargs['id']))


class PrimerSelectionView(CreatePlusView):
    template_name = 'primer-selection.html'
    form_class = PrimerSelectionForm
    success_url = '/main/primer-selection/{id}/experiment-summary/'

    def get_initial(self):
        primer_data = PrimerDesign.objects.get(id=self.kwargs['id']).primer_data

        def get_fwd_and_rev_primers(ontarget_primers: dict):
            values = list(ontarget_primers.values())
            return values[0], values[1]

        return {
            'selected_primers': dict(
                # TODO (gdingle): is this the best way of identifying guides?
                (p['target'] + ' ' + p['pam_id'],
                    get_fwd_and_rev_primers(p['ontarget_primers'])
                 if p['success'] else p['error'])
                for p in primer_data)
        }

    def plus(self, obj):
        obj.primer_design = PrimerDesign.objects.get(id=self.kwargs['id'])
        return obj

    def get_context_data(self, **kwargs):
        primer_data = PrimerDesign.objects.get(id=self.kwargs['id']).primer_data
        kwargs['example_crispor_url'] = primer_data[0]['url']
        kwargs['example_pam_id'] = primer_data[0]['pam_id']
        return super().get_context_data(**kwargs)


class ExperimentSummaryView(View):
    template_name = 'experiment-summary.html'

    def get(self, request, *args, **kwargs):
        try:
            primer_selection = PrimerSelection.objects.get(id=kwargs['id'])
        except PrimerSelection.DoesNotExist:
            # Alternate path to summary
            try:
                primer_selection = PrimerSelection.objects.filter(
                    primer_design__guide_selection__guide_design__experiment=kwargs['id'])[0]
            except IndexError:
                raise Http404('Experiment summary does not exist')

        sheet = samplesheet.from_primer_selection(primer_selection)
        sheet = self._prepare_sheet(sheet)

        primer_design = primer_selection.primer_design
        guide_selection = primer_design.guide_selection
        guide_design = guide_selection.guide_design
        experiment = guide_design.experiment

        return render(request, self.template_name, locals())

    def _prepare_sheet(self, sheet):
        """Modify sheet for optimal rendering"""
        sheet = sheet.loc[:, 'target_loc':]  # type: ignore
        sheet = sheet.loc[:, [not c.startswith('_') for c in sheet.columns]]
        sheet = sheet.dropna(axis=1, how='all')
        sheet.insert(0, 'well_pos', sheet.index)
        sheet.insert(1, 'well_num', range(1, len(sheet) + 1))
        sheet.columns = [c.replace('_', ' ').title() for c in sheet.columns]
        return sheet


class AnalysisView(CreatePlusView):
    template_name = 'analysis.html'
    form_class = AnalysisForm
    success_url = '/main/analysis/{id}/progress/'

    def plus(self, obj):
        # TODO (gdingle): use predetermined s3 location of fastq
        obj.fastqs = download_fastqs(obj.s3_bucket, obj.s3_prefix, overwrite=False)
        sheet = samplesheet.from_analysis(obj)
        batch = webscraperequest.CrispressoBatchWebRequest(obj)
        largs = [[
            row['target_seq'],
            row['guide_seq'],
            row['fastq_fwd'],
            row['fastq_rev'],
            row['donor_seq'],
            row.index
        ] for row in sheet.to_records()]
        batch.start(largs, [-1])
        return obj


class AnalysisProgressView(View):
    template_name = 'analysis-progress.html'
    success_url = '/main/analysis/{id}/results/'

    def get(self, request, **kwargs):
        analysis = Analysis.objects.get(id=kwargs['id'])
        batch_status = webscraperequest.CrispressoBatchWebRequest(analysis).get_batch_status()

        if not batch_status.is_successful:
            return render(request, self.template_name, locals())
        else:
            return HttpResponseRedirect(
                self.success_url.format(id=self.kwargs['id']))


class ResultsView(View):
    template_name = 'crispresso-results.html'

    def get(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=self.kwargs['id'])
        try:
            sheet = samplesheet.from_analysis(analysis)
        except IndexError:
            raise Http404('Analysis results do not exist')
        return render(request, self.template_name, locals())


# TODO (gdingle): is DetailView needed here?
class OrderFormView(DetailView):
    """
    Produces a downloadable Excel order form for IDT. The model must have a
    plate layout.
    """

    model: models.Model = None
    seq_key: str

    def _create_excel_file(self, sheet: samplesheet.pandas.DataFrame, title: str):
        wb = Workbook()
        ws = wb.active

        ws.title = title[0:31]  # Excel limits to 30 chars
        ws['A1'] = 'Well Position'
        ws['B1'] = 'Name'
        ws['C1'] = 'Sequence'

        for i, well_pos in enumerate(sheet.index):
            index = str(i + 2)
            ws['A' + index] = well_pos
            # TODO (gdingle): what is the best name of each well for order form?
            row = sheet.loc[well_pos]
            ws['B' + index] = '{}{}'.format(row.guide_offset, row.guide_direction)
            ws['C' + index] = row[self.seq_key]

        return writer.excel.save_virtual_workbook(wb)

    def get(self, request, *args, **kwargs):
        instance = self.model.objects.get(id=kwargs['id'])
        # TODO (gdingle): friendlier title?
        title = request.path.replace('/', ' ').replace('main ', '')
        excel_file = self._create_excel_file(instance.samplesheet, title)

        response = HttpResponse(
            excel_file, content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename="{}.xlsx"'.format(title)

        return response


class GuideOrderFormView(OrderFormView):

    model = GuideSelection
    # TODO: include guide_pam or not?
    seq_key = 'guide_seq'


class PrimerOrderFormView(OrderFormView):

    model = PrimerSelection
    # TODO (gdingle): IMPORTANT: how to order fwd and reverse primer at once?
    seq_key = 'primer_seq_fwd'
