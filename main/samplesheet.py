"""
This represents the contents of an experiment on a per sample, per well basis,
whereas the Django models represent contents per plate.
"""
import pandas

from main.models import *
from main.validators import get_guide_loc


def from_experiment(experiment: Experiment) -> pandas.DataFrame:
    sheet = _new_samplesheet()

    # Special pandas attribute for preserving metadata across transforms
    sheet._metadata = [
        'experiment_id',
        'experiment_name',
        # TODO (gdingle): use when avail
        # 'experiment_create_time',
        'experimenter_name',
    ]
    sheet.experiment_id = experiment.id
    sheet.experiment_name = experiment.name
    sheet.experimenter_name = experiment.researcher.full_name

    return sheet


def from_guide_selection(guide_selection: GuideSelection) -> pandas.DataFrame:
    guide_design = guide_selection.guide_design
    sheet = from_experiment(guide_design.experiment)

    # Ungroup guide data into rows
    target_loc_to_batch_id = dict((g['seq'], g['batch_id'])
                                  for g in guide_design.guide_data
                                  if g.get('batch_id'))  # filter out errors
    target_loc_to_target_seq = dict(zip(guide_design.targets, guide_design.target_seqs))
    guides = [(target_loc, offset, seq,
               target_loc_to_batch_id[target_loc],
               target_loc_to_target_seq[target_loc])
              for target_loc, selected in guide_selection.selected_guides.items()
              for offset, seq in selected.items()]

    # For assigning to subset of rows of A1 to H12
    lg = slice(0, len(guides))

    sheet['target_genome'] = guide_design.genome
    # TODO (gdingle): is this even worth having?
    sheet['target_pam'] = guide_design.pam

    sheet['target_loc'][lg] = [g[0] for g in guides]
    sheet['target_seq'][lg] = [g[4] for g in guides]

    # TODO: good to assume always 20 bp with trailing PAM as in Crispor?
    # TTCCGGCGCGCCGAGTCCTT AGG
    assert all(' ' in g[2] for g in guides), 'Expecting trailing PAM'
    assert all(len(g[2]) == 24 for g in guides), 'Expecting 20bp guide'
    sheet['guide_seq'][lg] = [g[2].split(' ')[0] for g in guides]
    sheet['guide_pam'][lg] = [g[2].split(' ')[1] for g in guides]

    # Taken direct from Crispor. Example: "s207-" and "s76+".
    # See http://crispor.tefor.net/manual/.
    # TODO (gdingle): recompute from guide_seq
    sheet['guide_offset'][lg] = [int(g[1][1:-1]) for g in guides]
    sheet['guide_direction'][lg] = [
        # TODO (gdingle): best naming?
        'fwd' if g[1][-1] == '+' else 'rev'
        for g in guides]

    # TODO (gdingle): is this correct? off by one?
    sheet['guide_loc'][lg] = sheet[lg].apply(
        lambda row: get_guide_loc(row['target_loc'], row['guide_offset'], len(row['guide_seq'])),
        axis=1,
    )

    sheet['_crispor_batch_id'][lg] = [g[3] for g in guides]
    sheet['_crispor_pam_id'][lg] = [g[1] for g in guides]
    # 'TODO_crispor_stats',

    # TODO (gdingle): donors
    # sheet['donor_seq'][lg] =
    # TODO (gdingle): how to compute this? talk to Jason Li
    # sheet['donor_target_seq'][lg] =

    # TODO (gdingle): is this wise to return subset?
    return sheet[lg]


def from_primer_selection(primer_selection: PrimerSelection) -> pandas.DataFrame:
    sheet = from_guide_selection(primer_selection.primer_design.guide_selection)
    selected_primers = primer_selection.selected_primers
    for primer_id, primer_pair in selected_primers.items():
        # TODO (gdingle): this is awkward
        target_loc, _crispor_pam_id = primer_id.split(' ')
        mask1 = sheet['target_loc'] == target_loc
        mask2 = sheet['_crispor_pam_id'] == _crispor_pam_id
        sheet['primer_seq_fwd'][mask1 & mask2] = primer_pair[0]
        sheet['primer_seq_rev'][mask1 & mask2] = primer_pair[1]
    # TODO (gdingle): is this really the best way to reindex?
    sheet = sheet.dropna(subset=['primer_seq_fwd'])
    sheet.index = _new_index(size=len(sheet))
    return sheet


def from_analysis(analysis: Analysis) -> pandas.DataFrame:
    primer_selection = PrimerSelection.objects.filter(
        primer_design__guide_selection__guide_design__experiment=analysis.experiment)[0]
    sheet = from_primer_selection(primer_selection)

    sheet._metadata = [
        'analysis_id',
        # TODO (gdingle): use when avail
        # 'analysis_create_time',
        'analyst_name',
    ]
    sheet.analysis_id = analysis.id
    sheet.analyst_name = analysis.researcher.full_name

    sheet['s3_bucket'] = analysis.s3_bucket
    sheet['s3_prefix'] = analysis.s3_prefix

    if not len(analysis.results_data):
        return sheet

    # TODO (gdingle): s3_key based on returned 'files'
    fastqs = analysis.results_data['fastqs']
    for i in range(0, len(fastqs), 2):
        # example: A3-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz
        sheet['fastq_fwd'] = fastqs[i].split('/')[-1]
        sheet['fastq_rev'] = fastqs[i + 1].split('/')[-1]

    results = analysis.results_data['results']
    # TODO (gdingle): match results by well location
    sheet['results_success'][0:len(results)] = [r[0] for r in results]
    sheet['results_path'][0:len(results)] = [r[1] for r in results]

    # Update target sequences only now because they are found by crispresso service
    amplicon_seqs = analysis.results_data['amplicon_seqs']
    sheet['target_seq'][0:len(amplicon_seqs)] = amplicon_seqs

    return sheet


def _new_index(size=96,
               end_char='H',
               end_int=12):
    chars = [chr(i) for i in range(ord('A'), ord(end_char) + 1)]
    ints = list(range(1, end_int + 1))
    assert len(chars) * len(ints) >= size
    return [c + str(i) for c in chars for i in ints][:size]


def _new_samplesheet():
    return pandas.DataFrame(
        index=_new_index(),
        # TODO (gdingle): use a MultiIndex?
        columns=[
            'target_genome',
            'target_pam',
            'target_loc',
            'target_seq',
            'guide_offset',
            'guide_loc',
            'guide_direction',
            'guide_seq',
            'guide_pam',
            '_crispor_batch_id',
            '_crispor_pam_id',
            # 'TODO_crispor_stats', off targets, etc
            'donor_seq',
            'donor_target_seq',
            # TODO (gdingle): does primer_loc have any meaning?
            # 'primer_loc',
            'primer_seq_fwd',
            'primer_seq_rev',
            # 'primer_melt_temp', # TODO (gdingle): is this useful?
            # TODO (gdingle): what's the source of truth for these?
            # 'well_pos',
            # 'well_num',
            's3_bucket',
            's3_prefix',
            'fastq_fwd',
            'fastq_rev',
            'results_success',
            'results_path',
            # results_stats_TODO
        ])


def to_illumina_sheet(samplesheet):
    # TODO (gdingle):
    return
