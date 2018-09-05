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
import os
import time

import requests
import requests_cache

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

from main import conversions
from main import samplesheet
from main.forms import *
from main.models import *
from main.validators import is_ensemble_transcript

# TODO (gdingle): move this when crispresso1 not needed
from crispresso.s3 import download_fastqs

# TODO (gdingle): move somewhere better
CRISPRESSO_ROOT_URL = 'http://crispresso:5000/'
CRISPRESSO_PUBLIC_ROOT_URL = 'http://0.0.0.0:5000/'


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
                # TODO (gdingle): does experiment name get us anything useful? aside from cache isolation per experiment?
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

        # Give some time for threads to finish to avoid GuideSelectionView too soon
        time.sleep(1)

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
            # Give some time for threads to finish updating database
            time.sleep(1)
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
        largs = sheet[['target_loc', '_crispor_pam_id', '_crispor_batch_id']].values
        for i, args in enumerate(largs):
            future = pool.submit(primers_request, args)
            future.add_done_callback(
                functools.partial(insert_primer_data, index=i))

        # Give some time for threads to finish to avoid PrimerSelectionView too soon
        time.sleep(1)

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
                    for row in sheet.to_records()]
        completed = [g for g in statuses if g[1]]
        incomplete = [g for g in statuses if not g[1]]
        assert len(statuses) == len(completed) + len(incomplete)

        if len(incomplete):
            return render(request, self.template_name, locals())
        else:
            # Give some time for threads to finish updating database
            time.sleep(1)
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


class AnalysisView(CreatePlusView):
    """
    # TODO (gdingle): replace AnalysisViewOLD
    """
    template_name = 'analysis.html'
    form_class = AnalysisForm
    success_url = '/main/analysis/{id}/progress/'

    def plus(self, obj):
        # TODO (gdingle): make use of precise s3 location of fastq, hdr amplicon
        obj.fastqs = download_fastqs(obj.s3_bucket, obj.s3_prefix, overwrite=False)

        # TODO (gdingle): pass in fastqs to sheet for rows
        sheet = samplesheet.from_analysis(obj)
        sheet = sheet[['target_seq', 'guide_seq', 'donor_seq', 'well_name']]

        self._start_all_analyses(sheet, obj)

        return obj

    @staticmethod
    def _start_all_analyses(sheet, obj):

        def analyze_fastq_pair(fastq_r1, fastq_r2, row):
            try:
                return crisporclient.CrispressoScrapeRequest(
                    row['target_seq'],
                    row['guide_seq'],
                    fastq_r1,
                    fastq_r2,
                    row['donor_seq'],
                    row['well_name']
                ).run()
            except Exception as e:
                return {
                    'success': False,
                    'error': e.args[0],
                }

        obj.results_data = [{'success': None}] * len(sheet)

        def insert_results_data(future, index=None):
            obj.results_data[index] = future.result()
            obj.save()

        # TODO (gdingle): optimal number of workers for crispresso2?
        pool = ThreadPoolExecutor(max_workers=4)
        for i, row in enumerate(sheet.to_records()):
            fastq_r1 = obj.fastqs[i * 2]
            fastq_r2 = obj.fastqs[i * 2 + 1]
            assert '_R1_' in fastq_r1 and '_R2_' in fastq_r2, 'Fastq files must be paired and sorted'
            # TODO (gdingle): match fastq files and selected_guides on sample names
            pool.submit(
                analyze_fastq_pair,
                fastq_r1,
                fastq_r2,
                row,
            ).add_done_callback(
                functools.partial(insert_results_data, index=i))


class AnalysisProgressView(View):
    template_name = 'analysis-progress.html'
    success_url = '/main/analysis/{id}/results/'

    # TODO (gdingle): decide whether to use requests_cache for status or database... then update other views
    # # TODO (gdingle): test different statuses
    # TODO (gdingle): why does first error wipe away running?
    def get(self, request, **kwargs):
        analysis = Analysis.objects.get(id=kwargs['id'])
        sheet = samplesheet.from_analysis(analysis)
        # TODO (gdingle): donor_seq
        sheet = sheet[['target_seq', 'guide_seq', 'well_name']]
        sheet.insert(0, 'well_pos', sheet.index)
        statuses, completed, running, errorred = [], [], [], []
        for i, row in enumerate(sheet.to_records()):
            statuses.append(row['well_name'])
            result = analysis.results_data[i]
            if result['success'] is True:
                completed.append(result['result_url'])
            elif result['success'] is False:
                errorred.append(result['error'])
            elif result['success'] is None:
                # TODO (gdingle): why isn't well_name working here?
                running.append((row['well_pos'], row['well_name']))

        return render(request, self.template_name, locals())

# TODO (gdingle): replace me
# class AnalysisProgressViewOLD(View):
#     template_name = 'analysis-progress.html'
#     success_url = '/main/analysis/{id}/results/'

#     def get(self, request, **kwargs):
#         analysis = Analysis.objects.get(id=kwargs['id'])
#         results_paths = analysis.results_data['results']

#         success_file = '/Quantification_of_editing_frequency.txt'
#         statuses = [(self._exists(path + success_file), path) for path in results_paths]
#         completed = [s[1] for s in statuses if s[0]]
#         incomplete = [s[1] for s in statuses if not s[0]]

#         # TODO (gdingle): temp hack for testing... need to replace with well_names goodness
#         def log_path(path):
#             return path.replace('output/', 'output/CRISPResso_on_') \
#                 + '/CRISPResso_RUNNING_LOG.txt'

#         # TODO (gdingle): allow set of running AND errorred
#         running = [path for path in results_paths if self._exists(log_path(path))]
#         if running:
#             errorred = self._get_errorred(results_paths, log_path)
#         return render(request, self.template_name, locals())

#     # TODO (gdingle): clean me up with final path structure
#     def _get_errorred(self, results_paths, log_path):
#         errorred = []
#         for path in results_paths:
#             error_lines = {
#                 line for line in
#                 self._get_log(log_path(path)).split('\n')
#                 if line.startswith('ERROR')
#             }
#             if error_lines:
#                 errorred.append((path, error_lines))
#         return errorred

#     @staticmethod
#     def _exists(path):
#         with requests_cache.disabled():
#             response = requests.head(CRISPRESSO_ROOT_URL + path)
#         assert response.status_code in (200, 404)
#         return response.status_code == 200

#     @staticmethod
#     def _get_log(path):
#         with requests_cache.disabled():
#             return requests.get(CRISPRESSO_ROOT_URL + path).text


# TODO (gdingle): replace me
# class AnalysisViewOLD(CreatePlusView):
#     template_name = 'analysis.html'
#     form_class = AnalysisForm
#     success_url = '/main/analysis/{id}/progress/'

#     def plus(self, obj):
#         sheet = samplesheet.from_analysis(obj)
#         # TODO (gdingle): make use of precise s3 location of fastq
#         sheet = sheet[['target_seq', 'guide_seq', 'well_name']]
#         data = {
#             'sheet': sheet.to_dict(),
#             's3_bucket': obj.s3_bucket,
#             's3_prefix': obj.s3_prefix,
#             'dryrun': False,
#             'async': True,
#         }

#         url = CRISPRESSO_ROOT_URL + 'analyze'  # host is name of docker service
#         # TODO: should we allow caching here as elsewhere?
#         with requests_cache.disabled():
#             response = requests.post(url, json=data)
#             response.raise_for_status()

#         obj.results_data = response.json()
#         return obj


class ResultsView(View):
    template_name = 'results.html'

    def get(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=self.kwargs['id'])
        results_urls = [CRISPRESSO_PUBLIC_ROOT_URL + path
                        for path in analysis.results_data['results']]
        return render(request, self.template_name, locals())


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
    seq_key = 'primer_seq_fwd'
