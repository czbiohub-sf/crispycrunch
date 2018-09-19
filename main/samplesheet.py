"""
This represents the contents of an experiment on a per sample, per well basis,
whereas the Django models represent contents per plate.

For tests, see main/tests.py.
"""
import os

from typing import Dict, Optional

import pandas

from main import conversions
from main.models import *
from main.validators import get_guide_loc, get_primer_loc


def from_experiment(experiment: Experiment) -> pandas.DataFrame:
    sheet = _new_samplesheet()

    # Special pandas attribute for preserving metadata across transforms
    # TODO (gdingle): this is not actually working :'(
    # sheet._metadata = [
    #     'experiment_id',
    #     'experiment_name',
    #     'experiment_description',
    #     'experiment_create_time',
    #     'experimenter_name',
    # ]
    # sheet.experiment_id = experiment.id
    # sheet.experiment_name = experiment.name
    # sheet.experiment_description = experiment.description
    # sheet.experiment_create_time = experiment.create_time
    # sheet.experimenter_name = experiment.researcher.full_name

    return sheet


def from_guide_selection(guide_selection: GuideSelection) -> pandas.DataFrame:
    guide_design = guide_selection.guide_design
    sheet = from_experiment(guide_design.experiment)

    # Ungroup guide data into rows
    # TODO (gdingle): this is a complicated beast that should at least have its own tests
    target_loc_to_batch_id = dict((g['target'], g['batch_id'])
                                  for g in guide_design.guide_data
                                  if g.get('batch_id'))  # filter out errors
    target_loc_to_target_seq = dict(zip(guide_design.targets, guide_design.target_seqs))
    selected_guides_ordered = [(g['target'], guide_selection.selected_guides[g['target']])
                               for g in guide_design.guide_data]
    guides = [(target_loc, offset, seq,
               target_loc_to_batch_id[target_loc],
               target_loc_to_target_seq[target_loc])
              for target_loc, selected in selected_guides_ordered
              for offset, seq in selected.items()]

    # For assigning to subset of rows of A1 to H12
    lg = slice(0, len(guides))

    sheet['target_genome'] = guide_design.genome
    sheet['target_pam'] = guide_design.pam

    sheet['target_loc'][lg] = [g[0] for g in guides]
    sheet['target_seq'][lg] = [g[4] for g in guides]

    # TODO: assume always 20 bp with trailing PAM as in Crispor?
    # TTCCGGCGCGCCGAGTCCTT AGG
    assert all(' ' in g[2] for g in guides), 'Expecting trailing PAM'
    assert all(len(g[2]) == 24 for g in guides), 'Expecting 20bp guide'
    sheet['guide_seq'][lg] = [g[2].split(' ')[0] for g in guides]
    sheet['guide_pam'][lg] = [g[2].split(' ')[1] for g in guides]

    # Taken direct from Crispor. Example: "s207-" and "s76+".
    # See http://crispor.tefor.net/manual/.
    # TODO (gdingle): recompute from guide_seq
    sheet['guide_offset'][lg] = [int(g[1][1:-1]) for g in guides]
    sheet['guide_direction'][lg] = [g[1][-1] for g in guides]

    # TODO (gdingle): is this correct? off by one?
    # TODO (gdingle): indicate direction by high to low for reverse?
    sheet['guide_loc'][lg] = sheet[lg].apply(
        lambda row: get_guide_loc(row['target_loc'], row['guide_offset'], len(row['guide_seq'])),
        axis=1,
    )

    # TODO (gdingle): is this really the best unique name?
    # Example: "hg38:chr2:136116735-136116754:-"
    sheet['well_name'][lg] = sheet[lg].apply(
        lambda row: '{}:{}:{}'.format(
            row['target_genome'],
            row['guide_loc'],
            row['guide_direction']),
        axis=1,
    )

    sheet['_crispor_batch_id'][lg] = [g[3] for g in guides]
    sheet['_crispor_pam_id'][lg] = [g[1] for g in guides]
    # 'TODO_crispor_stats',

    # TODO (gdingle): donor dna
    # sheet['donor_seq'][lg] =
    # TODO (gdingle): how to compute this? wait for protospacex
    # sheet['donor_target_seq'][lg] =

    return sheet[lg]


def from_primer_selection(primer_selection: PrimerSelection) -> pandas.DataFrame:
    sheet = from_guide_selection(primer_selection.primer_design.guide_selection)
    selected_primers = primer_selection.selected_primers
    for primer_id, primer_pair in selected_primers.items():
        # TODO (gdingle): this is awkward... should we use well_name?
        target_loc, _crispor_pam_id = primer_id.split(' ')
        mask1 = sheet['target_loc'] == target_loc
        mask2 = sheet['_crispor_pam_id'] == _crispor_pam_id
        sheet['primer_seq_fwd'][mask1 & mask2] = primer_pair[0][0]
        sheet['primer_seq_rev'][mask1 & mask2] = primer_pair[1][0]
        assert primer_pair[0][1].startswith(primer_pair[0][0]), 'Primer product should start with primer'
        # TODO (gdingle): do we need this after all?
        sheet['primer_product'][mask1 & mask2] = primer_pair[0][1]

    sheet = sheet.dropna(subset=['primer_seq_fwd'])
    sheet.index = _new_index(size=len(sheet))

    # TODO (gdingle): this is nearly the only IO in this file... do we really need it?
    sheet['primer_product'] = sheet.apply(
        lambda row: conversions.chr_loc_to_seq(
            get_primer_loc(row['primer_product'], row['guide_seq'], row['guide_loc']),
            row['target_genome'])
        # Only look up product from chr loc if crispor returns mysterious Ns
        if 'N' in row['primer_product'] else row['primer_product'],
        axis=1,
    )
    return sheet


def from_analysis(analysis: Analysis) -> pandas.DataFrame:
    primer_selection = PrimerSelection.objects.filter(
        primer_design__guide_selection__guide_design__experiment=analysis.experiment)[0]
    return from_analysis_and_primer_selection(analysis, primer_selection)


# TODO (gdingle): opportunity for multiple dispatch
def from_analysis_and_primer_selection(analysis: Analysis, primer_selection: PrimerSelection) -> pandas.DataFrame:
    sheet = from_primer_selection(primer_selection)

    sheet._metadata = [
        'analysis_id',
        'analysis_create_time',
        'analyst_name',
    ]
    sheet.analysis_id = analysis.id
    sheet.analysis_create_time = analysis.create_time
    sheet.analyst_name = analysis.researcher.full_name

    sheet['s3_bucket'] = analysis.s3_bucket
    sheet['s3_prefix'] = analysis.s3_prefix

    sheet = _insert_fastqs(sheet, analysis.fastqs)

    reports = [r for r in analysis.results_data]
    if not any(reports):
        return sheet

    sheet['report_url'] = [r.get('report_url') for r in reports]
    sheet['report_zip'] = [r.get('report_zip') for r in reports]

    sheet['report_stats'] = _drop_empty_report_stats(reports)

    # TODO (gdingle): strangely, we lose the order of headers as when compared to...
    # sheet['report_stats'] = [r.get('report_stats') for r in reports]
    # is it another instance of postgres json type ordered alphabetical?
    # TODO (gdingle): reapply original ordering:
    # ['Total', 'Unmodified', 'Modified', 'Discarded', 'Insertions', 'Deletions', 'Substitutions', 'Only Insertions', 'Only Deletions', 'Only Substitutions', 'Insertions and Deletions', 'Insertions and Substitutions', 'Deletions and Substitutions', 'Insertions Deletions and Substitutions']

    return sheet


def _drop_empty_report_stats(reports: list) -> Optional[Dict[str, int]]:
    """
    See https://stackoverflow.com/questions/21164910/delete-column-in-pandas-if-it-is-all-zeros
    """
    report_stats = [r.get('report_stats', {}) for r in reports]
    if not any(report_stats):
        return None
    temp_sheet = pandas.DataFrame.from_records(report_stats)
    nonzero_cols = (temp_sheet != 0).any(axis='rows')
    temp_sheet = temp_sheet.loc[:, nonzero_cols]
    return temp_sheet.to_dict(orient='records')


def _new_index(size=96,
               end_char='H',
               end_int=12) -> list:
    chars = [chr(i) for i in range(ord('A'), ord(end_char) + 1)]
    ints = list(range(1, end_int + 1))
    assert len(chars) * len(ints) >= size
    return [c + str(i) for c in chars for i in ints][:size]


def _new_samplesheet() -> pandas.DataFrame:
    return pandas.DataFrame(
        index=_new_index(),
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
            'donor_seq',
            'donor_target_seq',
            'primer_seq_fwd',
            'primer_seq_rev',
            # TODO (gdingle): rename to reference amplicon?
            'primer_product',
            'well_name',
            's3_bucket',
            's3_prefix',
            'fastq_fwd',
            'fastq_rev',
            'report_url',
            'report_zip',
            'report_stats',
        ])


def _insert_fastqs(sheet: pandas.DataFrame, fastqs: list) -> pandas.DataFrame:
    """
    Insert pairs of fastq sequence files into their corresponding rows.
    NOTE: This assumes a naming convention of A1, A2, ...H12 in the filename,
    separated by dashes (-), and Illumina "R1", "R2" for foward and reverse reads.

    For example: "A3-BCAP31-C-sorted-180212_S3_L001_R2_001.fastq.gz".
    """
    for well_pos in sheet.index:
        matches = sorted([
            filename for filename in fastqs
            if well_pos in os.path.basename(filename).split('-')
        ])
        if len(matches) == 0:
            # TODO (gdingle): allow missing fastqs?
            continue
        assert len(matches) == 2, 'Exactly two fastq files expected per well'
        sheet.loc[well_pos, 'fastq_fwd'] = matches[0]
        sheet.loc[well_pos, 'fastq_rev'] = matches[1]

    sheet = sheet.dropna(subset=['fastq_fwd', 'fastq_rev'])
    assert len(sheet), 'Some fastqs must match rows'
    return sheet


def to_illumina_sheet(samplesheet):
    # TODO (gdingle):
    return
