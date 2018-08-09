"""
This represents the contents of an experiment on a per sample, per well basis,
whereas the Django models represent contents per plate.

TODO: make a pandas dataframe
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

MultiIndex for metadata?

datafame schema

target_genome,
target_region,
target_seq,
guide_region,
guide_offset,
guide_seq,
guide_pam,
guide_direction,
TODO_crispor_stats,
donor_seq,
donor_target_seq,
primer_seq_fwd,
primer_seq_rev,
primer_melt_temp,
well_pos,
well_num,
aws_bucket,
aws_prefix,
fastq_basename,
fastq_r1,
fastq_r2,
results_status,
results_url,
results_stats_TODO

factory methods

from_experiment
from_guide_selection
from_primer_selection
from_analysis

to_plate
to_order_form
to_sample_sheet
to_experiment_summary

"""