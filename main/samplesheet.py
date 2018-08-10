"""
This represents the contents of an experiment on a per sample, per well basis,
whereas the Django models represent contents per plate.
"""

import pandas

# from crispresso.seqs import get_reference_amplicon
from main.models import GuideDesign

"""
TODO: use it in crispresso and crispycrunch services, or maybe just JSON crispresso

metadata schema:

analysis_id
analysis_name
analysis_create_time
analyst_name
"""


def from_experiment(experiment):
    # TODO (gdingle): assign metadata
    """
    experiment_id
    experiment_name
    experiment_create_time
    experimenter_name
    """
    sheet = _new_samplesheet()

    guide_design = GuideDesign.objects.get(experiment=experiment)
    # TODO (gdingle): assert finished experiment?
    sheet['target_genome'] = guide_design.genome
    targets = guide_design.targets
    sheet['target_loc'][0:len(targets)] = targets
    # TODO (gdingle): this does not work because pysam fails to compile... need conda :(
    # sheet['target_seq'][0:len(targets)] = (get_reference_amplicon(chr_loc) for chr_loc in targets)
    sheet['target_seq'] = 'TODO'
    # TODO (gdingle): is this worth having?
    sheet['target_pam'] = guide_design.pam
    return sheet


def from_guide_selection(guide_selection):
    sheet = from_experiment(guide_selection.guide_design.experiment)
    guides = [(chr_loc, offset, seq)
              for chr_loc, selected in guide_selection.selected_guides
              for offset, seq in selected]
    # TODO (gdingle): need to explode here!!!
    'guide_loc',
    'guide_offset',
    'guide_seq',
    'guide_pam',
    'guide_direction',
    # 'TODO_crispor_stats',
    'donor_seq',
    'donor_target_seq',

    # TODO (gdingle): shorter way to assign?
    sheet['guide_loc'][0:len(guides)] = [g[0] for g in guides]
    sheet['guide_offset'][0:len(guides)] = [g[1] for g in guides]
    sheet['guide_seq'][0:len(guides)] = [g[2] for g in guides]

    sheet['donor_loc'][0:len(guides)] = guides.keys()
    sheet['guide_seq'][0:len(guides)] = guides.values()
    return sheet


# def from_primer_selection(experiment):
#     return sheet


# def from_analysis(experiment):
#     return sheet


def _new_samplesheet(
        size=96,
        end_char='H',
        end_int=12):

    chars = [chr(i) for i in range(ord('A'), ord(end_char) + 1)]
    ints = list(range(1, end_int + 1))
    assert len(chars) * len(ints) == size

    return pandas.DataFrame(
        index=(c + str(i) for c in chars for i in ints),
        # TODO (gdingle): use a MultiIndex?
        columns=[
            'target_genome',
            'target_loc',
            'target_seq',
            'guide_loc',
            'guide_offset',
            'guide_seq',
            'guide_pam',
            'guide_direction',
            # 'TODO_crispor_stats',
            'donor_seq',
            'donor_target_seq',
            'primer_loc',
            'primer_seq_fwd',
            'primer_seq_rev',
            'primer_melt_temp',
            'well_pos',
            'well_num',
            'aws_bucket',
            'aws_prefix',
            'fastq_basename',
            'fastq_fwd',
            'fastq_rev',
            'results_status',
            'results_url',
            # results_stats_TODO
        ])


def to_plate(samplesheet):
    # TODO (gdingle):
    return 'Plate96Layout'


def to_order_form(samplesheet):
    # TODO (gdingle):
    return 'OrderFormView'


def to_illumina_sheet(samplesheet):
    # TODO (gdingle):
    return


def to_experiment_summary(samplesheet):
    # TODO (gdingle):
    return
