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

from main import samplesheet
# TODO (gdingle): remove me?
from main import conversions
from main.forms import *
from main.models import *
from main.validators import is_ensemble_transcript


def index(request):
    # TODO (gdingle): create useful index
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
    # TODO (gdingle): refactor success_urls to models for reuse
    success_url = '/main/experiment/{id}/guide-design/'


class GuideDesignView(CreatePlusView):
    template_name = 'guide-design.html'
    form_class = GuideDesignForm
    # TODO (gdingle): redirect to waiting page
    success_url = '/main/guide-design/{id}/progress/'

    def _normalize_targets(self, targets):
        # TODO (gdingle): handle mix of target types
        if not all(is_gene(t) for t in targets):
            return targets

        with ThreadPoolExecutor() as pool:
            normalized = list(pool.map(
                conversions.gene_to_chr_loc,
                targets,
            ))
        # TODO: normalize seqs also
        assert all(is_chr(t)for t in normalized)
        return normalized

    def _get_target_seqs(self, targets, genome):
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
        # def tagin_request(target):
        #     return crisporclient.TagInRequest(
        #         target,
        #         tag=obj.tag_in,
        #         # species # TODO (gdingle): translate from crispor
        #     ).run()
        # if obj.hdr_seq:
        #     # TODO (gdingle): put in form validation somehow
        #     assert all(is_ensemble_transcript(t) and len(t) <= 600 for t in obj.targets), 'Bad input for TagIn'
        #     obj.donor_data = list(ex.map(tagin_request, obj.targets))
        #     # Crispor does not accept Ensembl transcript IDs
        #     # and use guide_chr_range to avoid 2000 bp limit
        #     # TODO (gdingle): is this wise?
        #     crispor_targets = [d['metadata']['guide_chr_range'] for d in obj.donor_data]

        def guide_request(target):
            return crisporclient.CrisporGuideRequest(
                target,
                # TODO (gdingle): does experiment name get us anything useful?
                name=obj.experiment.name,
                org=obj.genome,
                pam=obj.pam).run()

        # More than 8 threads appears to cause a 'no output' Crispor error
        pool = ThreadPoolExecutor(8)

        def insert_guide_data(future, index=None):
            obj.guide_data[index] = future.result()
            obj.save()

        obj.guide_data = [{}] * len(obj.targets)
        for i, target in enumerate(obj.targets):
            future = pool.submit(guide_request, target)
            future.add_done_callback(
                functools.partial(insert_guide_data, index=i))

        return obj


class GuideDesignProgressView(View):

    template_name = 'guide-design-progress.html'
    success_url = '/main/guide-design/{id}/guide-selection/'

    def get(self, request, **kwargs):
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])

        # See also guide_request above. These should match. TODO: refactor.
        def guide_request(target):
            return crisporclient.CrisporGuideRequest(
                target,
                name=guide_design.experiment.name,
                org=guide_design.genome,
                pam=guide_design.pam)

        statuses = [
            (target, guide_request(target).in_cache())
            for target in guide_design.targets]
        completed = [target for target, status in statuses if status]
        incomplete = [target for target, status in statuses if not status]
        assert len(completed) + len(incomplete) == len(statuses)

        if len(incomplete):
            return render(request, self.template_name, locals())
        else:
            return HttpResponseRedirect(
                self.success_url.format(id=self.kwargs['id']))


class GuideSelectionView(CreatePlusView):
    template_name = 'guide-selection.html'
    form_class = GuideSelectionForm
    success_url = '/main/guide-selection/{id}/primer-design/'

    def get_initial(self):
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        return {
            'selected_guides': dict((g['seq'], g['guide_seqs'])
                                    for g in guide_design.guide_data if g),
            'selected_donors': dict((g['metadata']['chr_loc'], g['donor_seqs'])
                                    for g in guide_design.donor_data),
            # TODO (gdingle): temp for debuggin
            'selected_guides_tagin': dict((g['metadata']['chr_loc'], g['guide_seqs'])
                                          for g in guide_design.donor_data),
        }

    def plus(self, obj):
        obj.guide_design = GuideDesign.objects.get(id=self.kwargs['id'])

        # TODO (gdingle): custom validate selected guides
        # not found, error....all must have sequences
        return obj

    def get_context_data(self, **kwargs):
        # TODO (gdingle): make multi guide
        guide_design = GuideDesign.objects.get(id=self.kwargs['id'])
        kwargs['crispor_url'] = [
            gd['url']
            for gd in guide_design.guide_data
            # TODO (gdingle): this line is failing sometimes
            if gd.get('url')][0]
        if guide_design.donor_data:
            kwargs['tagin_url'] = guide_design.donor_data[0]['url']
        return super().get_context_data(**kwargs)


# TODO (gdingle): do we still want this?
# class GuidePlateLayoutView(CreatePlusView):
#     template_name = 'guide-plate-layout.html'
#     form_class = GuidePlateLayoutForm
#     success_url = '/main/guide-selection/{guide_selection_id}/primer-design/'

#     def plus(self, obj):
#         obj.guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
#         return obj

#     def get_context_data(self, **kwargs):
#         guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
#         # TODO (gdingle): show targets in plate or sgRNA?
#         layout_map = dict((list(guides.values())[0], target)
#                           for target, guides
#                           in guide_selection.selected_guides.items())
#         kwargs['plate_layout'] = Plate96Layout(layout_map)
#         return super().get_context_data(**kwargs)


class PrimerDesignView(CreatePlusView):
    template_name = 'primer-design.html'
    form_class = PrimerDesignForm
    success_url = '/main/primer-design/{id}/progress/'

    def plus(self, obj):
        guide_selection = GuideSelection.objects.get(id=self.kwargs['id'])
        obj.guide_selection = guide_selection

        def primers_request(args):
            seq, pam_id, batch_id = args
            return crisporclient.CrisporPrimerRequest(
                batch_id=batch_id,
                amp_len=obj.max_amplicon_length,
                tm=obj.primer_temp,
                pam=guide_selection.guide_design.pam,
                pam_id=pam_id,
                seq=seq).run()

        def insert_primer_data(future, index=None):
            obj.primer_data[index] = future.result()
            obj.save()

        sheet = samplesheet.from_guide_selection(guide_selection)
        obj.primer_data = [{}] * len(sheet)
        pool = ThreadPoolExecutor()
        # TODO (gdingle): _crispor_batch_id has a lifetime... handle expired
        largs = sheet[['target_loc', '_crispor_pam_id', '_crispor_batch_id']].values
        for i, args in enumerate(largs):
            future = pool.submit(primers_request, args)
            future.add_done_callback(
                functools.partial(insert_primer_data, index=i))

        # TODO (gdingle): run crispr-primer if HDR experiment
        # https://github.com/chanzuckerberg/crispr-primer
        return obj


class PrimerDesignProgressView(View):
    template_name = 'primer-design-progress.html'
    success_url = '/main/primer-design/{id}/primer-selection/'

    def get(self, request, **kwargs):
        primer_design = PrimerDesign.objects.get(id=kwargs['id'])
        sheet = samplesheet.from_guide_selection(primer_design.guide_selection)

        # See also primers_request in PrimerDesignView. TODO: refactor
        def primers_request(row):
            return crisporclient.CrisporPrimerRequest(
                batch_id=row['_crispor_batch_id'],
                amp_len=primer_design.max_amplicon_length,
                tm=primer_design.primer_temp,
                pam=primer_design.guide_selection.guide_design.pam,
                pam_id=row['_crispor_pam_id'],
                seq=row['target_loc'])

        statuses = [(row._crispor_batch_id, row._crispor_pam_id, primers_request(row).in_cache())
                    for _, row in sheet.iterrows()]
        completed = [g for g in statuses if g[1]]
        incomplete = [g for g in statuses if not g[1]]
        assert len(statuses) == len(completed) + len(incomplete)

        if len(incomplete):
            return render(request, self.template_name, locals())
        else:
            return HttpResponseRedirect(
                self.success_url.format(id=kwargs['id']))


class PrimerSelectionView(CreatePlusView):
    template_name = 'primer-selection.html'
    form_class = PrimerSelectionForm
    success_url = '/main/primer-selection/{id}/experiment-summary/'

    def get_initial(self):
        primer_data = PrimerDesign.objects.get(id=self.kwargs['id']).primer_data

        def get_fwd_and_rev_primers(ontarget_primers):
            values = list(ontarget_primers.values())
            return values[0], values[1]

        return {
            # TODO (gdingle): IMPORTANT... need a unique ID ... batch_id + pam_id, or chr_loc + pam_id, or guide_loc
            'selected_primers': dict(
                (p['seq'] + ' ' + p['pam_id'],
                    get_fwd_and_rev_primers(p['ontarget_primers']))
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


# TODO (gdingle): do we still want this?
# class PrimerPlateLayoutView(CreatePlusView):
#     template_name = 'primer-plate-layout.html'
#     form_class = PrimerPlateLayoutForm
#     success_url = '/main/primer-plate-layout/{id}/experiment-summary/'

#     def plus(self, obj):
#         obj.primer_selection = PrimerSelection.objects.get(id=self.kwargs['id'])
#         return obj

#     def get_context_data(self, **kwargs):
#         primer_selection = PrimerSelection.objects.get(id=self.kwargs['id'])
#         primer_pairs = dict((tuple(p.keys()), tuple(p.values()))
#                             for target, p in primer_selection.selected_primers.items())
#         kwargs['plate_layout'] = Plate96Layout(primer_pairs)
#         return super().get_context_data(**kwargs)


class AnalysisView(CreatePlusView):
    template_name = 'analysis.html'
    form_class = AnalysisForm
    success_url = '/main/analysis/{id}/results/'

    def plus(self, obj):
        # TODO (gdingle): replace this with sheet.to_json()
        # which includes target_seqs
        data = {
            'selected_guides': obj.get_selected_guides(),
            'selected_donors': obj.get_selected_donors(),
            's3_bucket': obj.s3_bucket,
            's3_prefix': obj.s3_prefix,
            # TODO (gdingle): turn back on after we have a some real end-to-end
            # results to test
            'dryrun': False,
        }

        url = 'http://crispresso:5000/analyze'  # host is name of docker service
        response = requests.post(url, json=data)
        response.raise_for_status()

        obj.results_data = response.json()
        return obj


#
# END EXPERIMENT CREATION VIEWS
#

class ResultsView(View):
    template_name = 'results.html'

    def get(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=self.kwargs['id'])
        # TODO (gdingle): this implies that crispresso will have a public web server
        # do we want to upload all files to s3 instead? copy them to local django static dir?
        crispresso_root_url = 'http://0.0.0.0:5000/'
        # TODO (gdingle): temp... remove me
        guide_seq = 'AATCGGTACAAGATGGCGGA'
        return render(request, self.template_name, locals())


class ExperimentSummaryView(View):
    template_name = 'experiment-summary.html'

    def get(self, request, *args, **kwargs):
        primer_selection = PrimerSelection.objects.get(id=kwargs['id'])

        sheet = samplesheet.from_primer_selection(primer_selection)
        sheet = self._prepare_sheet(sheet)

        primer_design = primer_selection.primer_design
        guide_selection = primer_design.guide_selection
        guide_design = guide_selection.guide_design
        experiment = guide_design.experiment

        return render(request, self.template_name, locals())

    def _prepare_sheet(self, sheet):
        """Modify sheet for optimal rendering"""
        sheet = sheet.loc[:, 'target_loc':]
        sheet = sheet.loc[:, [not c.startswith('_') for c in sheet.columns]]
        sheet = sheet.dropna(axis=1, how='all')
        sheet.insert(0, 'well_pos', sheet.index)
        sheet.insert(1, 'well_num', range(1, len(sheet) + 1))
        sheet.columns = [c.replace('_', ' ').title() for c in sheet.columns]
        return sheet


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
    # TODO (gdingle): how to order fwd and reverse primer at once?
    # TODO (gdingle): what to do about "not found" primers?
    seq_key = 'primer_seq_fwd'
