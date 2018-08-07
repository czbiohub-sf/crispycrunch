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
import requests

from typing import no_type_check

from concurrent.futures import ThreadPoolExecutor
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render
from django.views import View
from django.views.generic import DetailView
from django.views.generic.edit import CreateView
from itertools import islice
from openpyxl import Workbook, writer  # noqa

import crisporclient

from main.conversions import convert_chr_to_fasta
from main.forms import *
from main.models import *
from main.platelayout import *
from main.validators import is_ensemble_transcript


def index(request):
    # TODO (gdingle): create useful index
    return HttpResponse("Hello, world. You're at the main index.")


def crispresso(request):
    # TODO (gdingle): switch to s3... see https://github.com/etianen/django-s3-storage

    # TODO (gdingle): rename and move
    data = {
        's3_bucket': 'jasonli-bucket',
        's3_prefix': 'JasonHDR/96wp1sorted-fastq/',
        'amplicon_seq': 'cgaggagatacaggcggagggcgaggagatacaggcggagggcgaggagatacaggcggagagcgGCGCTAGGACCCGCCGGCCACCCCGCCGGCTCCCGGGAGGTTGATAAAGCGGCGGCGGCGTTTGACGTCAGTGGGGAGTTAATTTTAAATCGGTACAAGATGGCGGAGGGGGACGAGGCAGCGCGAGGGCAGCAACCGCACCAGGGGCTGTGGCGCCGGCGACGGACCAGCGACCCAAGCGCCGCGGTTAACCACGTCTCGTCCAC',  # noqa
        'guide_seq': 'AATCGGTACAAGATGGCGGA',
        'expected_hdr_amplicon_seq': 'cgaggagatacaggcggagggcgaggagatacaggcggagggcgaggagatacaggcggagagcgGCGCTAGGACCCGCCGGCCACCCCGCCGGCTCCCGGGAGGTTGATAAAGCGGCGGCGGCGTTTGACGTCAGTGGGGAGTTAATTTTAAATCGGTACAAGATGCGTGACCACATGGTCCTTCATGAGTATGTAAATGCTGCTGGGATTACAGGTGGCGGAttggaagttttgtttcaaggtccaggaagtggtGCGGAGGGGGACGAGGCAGCGCGAGGGCAGCAACCGCACCAGGGGCTGTGGCGCCGGCGACGGACCAGCGACCCAAGCGCCGCGGTTAACCACGTCTCGTCCAC',  # noqa
        'dryrun': True,
    }

    url = 'http://crispresso:5000/crispresso'  # host is name of docker service
    # TODO (gdingle): switch to post
    response = requests.get(url, params=data)
    response.raise_for_status()
    return HttpResponse(response.text)

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
    # TODO (gdingle): refactor success_urls to models for reuse
    success_url = '/main/experiment/{id}/guide-design/'


class GuideDesignView(CreatePlusView):
    template_name = 'guide-design.html'
    form_class = GuideDesignForm
    success_url = '/main/guide-design/{id}/guide-selection/'

    def plus(self, obj):
        """
        If an HDR tag-in experiment, get donor DNA then get guides. If not, just
        get guides.
        """

        obj.experiment = Experiment.objects.get(id=self.kwargs['id'])

        def tagin_request(target):
            return crisporclient.TagInRequest(
                target,
                tag=obj.tag_in,
                # species # TODO (gdingle): translate from crispor
            ).run()

        def guide_request(target):
            return crisporclient.CrisporGuideRequest(
                target,
                # TODO (gdingle): does experiment name get us anything useful?
                name=obj.experiment.name,
                org=obj.genome,
                pam=obj.pam).run()

        with ThreadPoolExecutor() as ex:
            if obj.tag_in:
                # TODO (gdingle): put in form validation somehow
                assert all(is_ensemble_transcript(t) and len(t) <= 600 for t in obj.targets), 'Bad input for TagIn'
                obj.donor_data = list(ex.map(tagin_request, obj.targets))
                # Crispor does not accept Ensembl transcript IDs
                # and use guide_chr_range to avoid 2000 bp limit
                # TODO (gdingle): is this wise?
                crispor_targets = [d['metadata']['guide_chr_range'] for d in obj.donor_data]
            else:
                crispor_targets = obj.targets

            obj.guide_data = list(ex.map(guide_request, crispor_targets))
            # TODO (gdingle): is this even useful?
            # obj.target_fastas = list(ex.map(convert_chr_to_fasta, filter(is_chr, obj.targets)))

        return obj


class GuideSelectionView(CreatePlusView):
    template_name = 'guide-selection.html'
    form_class = GuideSelectionForm
    success_url = '/main/guide-selection/{id}/guide-plate-layout/'

    def get_initial(self):
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        return {
            'selected_guides': dict((g['seq'], g['guide_seqs'])
                                    for g in guide_design.guide_data),
            'selected_donors': dict((g['metadata']['chr_loc'], g['donor_seqs'])
                                    for g in guide_design.donor_data),
            # TODO (gdingle): temp for debuggin
            'selected_guides_tagin': dict((g['metadata']['chr_loc'], g['guide_seqs'])
                                          for g in guide_design.donor_data),
        }

    def plus(self, obj):
        obj.guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        # TODO (gdingle): this is awkward
        obj.selected_guides = dict((target, dict(islice(guides.items(), 1)))
                                   for target, guides in obj.selected_guides.items())
        return obj

    def get_context_data(self, **kwargs):
        # TODO (gdingle): make multi guide
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        kwargs['crispor_url'] = guide_design.guide_data[0]['url']
        if guide_design.donor_data:
            kwargs['tagin_url'] = guide_design.donor_data[0]['url']
        return super().get_context_data(**kwargs)


class GuidePlateLayoutView(CreatePlusView):
    template_name = 'guide-plate-layout.html'
    form_class = GuidePlateLayoutForm
    success_url = '/main/guide-selection/{guide_selection_id}/primer-design/'

    def plus(self, obj):
        obj.guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        return obj

    def get_context_data(self, **kwargs):
        guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        # TODO (gdingle): show targets in plate or sgRNA?
        layout_map = dict((list(guides.values())[0], target)
                          for target, guides
                          in guide_selection.selected_guides.items())
        kwargs['plate_layout'] = Plate96Layout(layout_map)
        return super().get_context_data(**kwargs)


class PrimerDesignView(CreatePlusView):
    template_name = 'primer-design.html'
    form_class = PrimerDesignForm
    success_url = '/main/primer-design/{id}/primer-selection/'

    def plus(self, obj):
        guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])

        guide_data = guide_selection.guide_design.guide_data

        def request_primers(args):
            seq, pam_id, batch_id = args
            return crisporclient.CrisporPrimerRequest(
                batch_id=batch_id,
                amp_len=obj.max_amplicon_length,
                tm=obj.primer_temp,
                pam=guide_selection.guide_design.pam,
                pam_id=pam_id,
                seq=seq).run()

        with ThreadPoolExecutor() as ex:
            # Take first selected guide of each target only
            # TODO (gdingle): this is awkward... keep batch_id and seq in a better way
            args = ((target, list(pam_ids.keys())[0], g['batch_id'])
                    for target, pam_ids
                    in guide_selection.selected_guides.items()
                    for g in guide_data
                    if g['seq'] == target)
            obj.primer_data = list(ex.map(request_primers, args))

        obj.guide_selection = guide_selection
        # TODO (gdingle): run crispr-primer if HDR experiment
        # https://github.com/chanzuckerberg/crispr-primer
        return obj


class PrimerSelectionView(CreatePlusView):
    template_name = 'primer-selection.html'
    form_class = PrimerSelectionForm
    success_url = '/main/primer-selection/{id}/primer-plate-layout/'

    def get_initial(self):
        primer_data = PrimerDesign.objects.get(id=self.kwargs['id']).primer_data
        return {
            # Andy's group is only interested in ontarget primer pairs, two per well
            # TODO (gdingle): is this true of other groups?
            # TODO (gdingle): put both primers in one well
            'selected_primers': dict((p['seq'], p['ontarget_primers']) for p in primer_data)
        }

    def plus(self, obj):
        obj.primer_design = PrimerDesign.objects.get(id=self.kwargs['id'])
        return obj

    def get_context_data(self, **kwargs):
        # TODO (gdingle): make multi guide
        primer_data = PrimerDesign.objects.get(id=self.kwargs['id']).primer_data
        kwargs['crispor_url'] = primer_data[0]['url']
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
        primer_pairs = dict((tuple(p.keys()), tuple(p.values()))
                            for target, p in primer_selection.selected_primers.items())
        kwargs['plate_layout'] = Plate96Layout(primer_pairs)
        return super().get_context_data(**kwargs)


class AnalysisView(CreateView):
    template_name = 'analysis.html'
    form_class = AnalysisForm
    success_url = '/main/analysis/{id}/results/'


#
# END EXPERIMENT CREATION VIEWS
#

class ResultsView(View):
    template_name = 'results.html'

    def get(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=self.kwargs['id'])
        return render(request, self.template_name, locals())


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

        rows = self._get_rows(guide_plate_layout, primer_plate_layout)

        return render(request, self.template_name, locals())

    @no_type_check  # TODO (gdingle): why cant work with variable tuple?
    def _get_rows(self, guide_plate_layout, primer_plate_layout) -> List[tuple]:
        rows = []
        for pos, guide in guide_plate_layout.layout.well_names.items():
            if guide is None:
                break
            rows.append((
                pos,
                guide,
                guide_plate_layout.layout.well_seqs[pos],
                primer_plate_layout.layout.well_seqs[pos],
                # TODO (gdingle): donor_seqs
            ))
        return rows


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
            # TODO (gdingle): this is awkward
            pair = plate_layout.well_seqs[pos]
            ws['C' + index] = list(pair.values())[0] if pair is not None else None

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
