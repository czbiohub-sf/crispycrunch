"""
This represents the contents of an experiment on a per sample, per well basis,
whereas the Django models represent contents per plate.
"""

import doctest
import pandas

"""
TODO: use it in crispresso and crispycrunch services, or maybe just JSON crispresso

metadata schema:

experiment_id
experiment_name
experiment_create_time
experimenter_name
analysis_id
analysis_name
analysis_create_time
analyst_name
"""

def from_experiment(experiment):
    """
    >>> experiment = Experiment.objects.get(name='testsum3')
    >>> samplesheet = from_experiment(experiment)
    """
    sheet = _new_samplesheet()
    sheet['target_genome'] = 'TODO'
    return sheet


# def from_guide_selection(experiment):
#     sheet = _new_samplesheet()
#     sheet['target_genome'] = 'TODO'
#     return sheet


# def from_primer_selection(experiment):
#     sheet = _new_samplesheet()
#     sheet['target_genome'] = 'TODO'
#     return sheet


# def from_analysis(experiment):
#     sheet = _new_samplesheet()
#     sheet['target_genome'] = 'TODO'
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
            'target_region',
            'target_seq',
            'guide_region',
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


if __name__ == '__main__':
    doctest.testmod()
