"""
This represents the contents of an experiment on a per sample, per well basis,
whereas the Django models represent contents per plate.

For tests, see main/tests.py.
"""
import logging
import os

from io import BytesIO
from typing import Dict, List, Optional, Tuple

import pandas

from pandas import DataFrame

from django.core.files.uploadedfile import UploadedFile

from main.models import Analysis, Experiment, GuideDesign, GuideSelection, PrimerSelection
from utils import conversions
from utils import hdr
from utils.chrloc import ChrLoc, get_guide_cut_to_insert, get_guide_loc, get_primer_loc

# TODO (gdingle): move to conversions
from crispresso.fastqs import reverse_complement

logger = logging.getLogger(__name__)


def from_experiment(experiment: Experiment) -> DataFrame:
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


def from_guide_selection(guide_selection: GuideSelection) -> DataFrame:
    guide_design = guide_selection.guide_design
    sheet = from_experiment(guide_design.experiment)

    guides = _flatten_guide_data(guide_selection)

    # Trim sheet to available guides
    sheet = sheet[0:len(guides)]

    sheet['target_genome'] = guide_design.genome
    sheet['target_pam'] = guide_design.pam

    sheet['target_loc'] = [ChrLoc(g[0]) for g in guides]
    sheet['target_seq'] = [g[4] for g in guides]

    # TTCCGGCGCGCCGAGTCCTT AGG
    assert all(' ' in g[2] for g in guides), 'Expecting trailing PAM'
    assert all(len(g[2]) == 24 for g in guides), 'Expecting 20bp guide'
    sheet['guide_seq'] = [g[2].split(' ')[0] for g in guides]
    sheet['guide_pam'] = [g[2].split(' ')[1] for g in guides]

    # Taken direct from Crispor. Example: "s207-" and "s76+".
    # See http://crispor.tefor.net/manual/.
    # DEFINITION: number of chars before the first char of NGG PAM
    sheet['guide_offset'] = [int(g[1][1:-1]) for g in guides]
    sheet['_guide_direction'] = [g[1][-1] for g in guides]

    # TODO (gdingle): is this correct? off by one? reverse strand?
    sheet['guide_loc'] = sheet.apply(
        lambda row: get_guide_loc(
            row['target_loc'],
            row['guide_offset'],
            len(row['guide_seq']),
            row['_guide_direction']),
        axis=1,
    )

    sheet = _set_scores(sheet, guide_design, guides)

    # TODO (gdingle): put into unit test
    if guide_design.hdr_seq:
        sheet = _set_hdr_cols(sheet, guide_design.hdr_seq, guide_design.hdr_tag)

    # TODO (gdingle): is this really the best unique name?
    # Example: "hg38:chr2:136116735-136116754:-"
    # TODO (gdingle): ever useful?
    # sheet['well_name'] = sheet.apply(
    #     lambda row: '{}:{}:{}'.format(
    #         row['target_genome'],
    #         row['guide_loc'],
    #         row['_guide_direction']),
    #     axis=1,
    # )

    sheet['_crispor_batch_id'] = [g[3] for g in guides]
    sheet['_crispor_pam_id'] = [g[1] for g in guides]
    # TODO (gdingle): is this the best way of identifying guides?
    sheet['_crispor_guide_id'] = [g[0] + ' ' + g[1] for g in guides]

    return sheet


def from_primer_selection(primer_selection: PrimerSelection) -> DataFrame:
    guide_selection = primer_selection.primer_design.guide_selection
    sheet = from_guide_selection(guide_selection)
    selected_primers = primer_selection.selected_primers
    for guide_id, primer_pair in selected_primers.items():
        # TODO (gdingle): this is awkward... should we use well_name?
        target_loc, _crispor_pam_id = guide_id.split(' ')
        mask1 = sheet['target_loc'] == target_loc
        mask2 = sheet['_crispor_pam_id'] == _crispor_pam_id
        assert any(mask1) and any(mask2)

        # TODO (gdingle): A value is trying to be set on a copy of a slice from a DataFrame
        # See the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
        sheet['primer_seq_fwd'][mask1 & mask2] = primer_pair[0][0]
        sheet['primer_seq_rev'][mask1 & mask2] = primer_pair[1][0]
        assert primer_pair[0][1].startswith(
            primer_pair[0][0]), 'Primer product should start with forward primer'
        sheet['primer_product'][mask1 & mask2] = primer_pair[0][1]

    sheet = sheet.dropna(subset=['primer_seq_fwd'])
    assert len(sheet)
    sheet.index = _new_index(size=len(sheet))

    sheet['primer_product'] = sheet.apply(_transform_primer_product, axis=1)

    if guide_selection.guide_design.hdr_seq:
        sheet = _set_hdr_primer(
            sheet,
            guide_selection.guide_design,
            primer_selection.primer_design.max_amplicon_length)

    return sheet


def _set_hdr_primer(sheet: DataFrame, guide_design: GuideDesign, max_amplicon_length: int):
    # TODO (gdingle): ryan leenay says primer needs to be codon-aware
    # TODO (gdingle): figure out primers and HDR... align to guide then to hdr_dist?
    # sheet['primer_product'] = 'need to align to codon'
    # return sheet

    def get_primer_product(row):
        primer_product = row['primer_product']

        if ' ' in primer_product:
            # previous warning
            return primer_product

        guide_seq_aligned = hdr.HDR(
            row['target_seq'],
            guide_design.hdr_seq,
            guide_design.hdr_tag,
            row['hdr_dist'],
            row['_guide_direction']).guide_seq_aligned

        # Crispor returns primer products by strand. Normalize to positive strand.
        if row['target_loc'].strand == '-':
            primer_product = reverse_complement(primer_product)

        guide_offset = primer_product.find(guide_seq_aligned)

        if guide_offset == -1:
            logger.warning('Could not find guide {} in primer {} for target {}'.format(
                row['guide_seq'], primer_product, row['target_loc']))
            return 'guide not found: ' + primer_product

        start = guide_offset % 3
        before, primer_product_aligned = \
            primer_product[:start], primer_product[start:]
        assert before + primer_product_aligned == primer_product

        hdr_primer_product = hdr.HDR(
            primer_product_aligned,
            guide_design.hdr_seq,
            guide_design.hdr_tag,
            row['hdr_dist'],
            row['_guide_direction']).template
        assert len(before) + len(hdr_primer_product) == len(primer_product) + \
            len(guide_design.hdr_seq)

        # TODO: Normalize back to direction of Crispor return?
        return before + hdr_primer_product

    # TODO (gdingle): fix root cause of issues
    def warn_hdr_primer(row) -> str:
        primer_product = row['primer_product']

        if ' ' in primer_product:
            # previous warning
            return primer_product

        plen = len(primer_product)
        if plen > max_amplicon_length:
            return f'too long, {plen}bp: {primer_product}'

        arms = primer_product.upper().split(guide_design.hdr_seq)
        larm, rarm = len(arms[0]), len(arms[1])
        if min(larm, rarm) < 55:
            return 'homology arm too short, {}bp: {}'.format(min(larm, rarm), primer_product)

        return primer_product

    sheet['primer_product'] = sheet.apply(get_primer_product, axis=1)
    sheet['primer_product'] = sheet.apply(warn_hdr_primer, axis=1)
    return sheet


def _flatten_guide_data(
    guide_selection: GuideSelection
) -> List[Tuple[str, str, str, str, str]]:

    # TODO (gdingle): try using bolton itertools remap?

    guide_design = guide_selection.guide_design
    # Ungroup guide data into rows
    # TODO (gdingle): this is a complicated beast that should at least have its own tests
    target_loc_to_batch_id = dict((g['target'], g['batch_id'])
                                  for g in guide_design.guide_data
                                  if g.get('batch_id'))  # filter out errors
    target_loc_to_target_seq = dict(zip(guide_design.targets, guide_design.target_seqs))
    selected_guides_ordered = [(g['target'], guide_selection.selected_guides[g['target']])
                               for g in guide_design.guide_data
                               if g['target'] in guide_selection.selected_guides]
    guides = [(target_loc, offset, seq,
               target_loc_to_batch_id[target_loc],
               target_loc_to_target_seq[target_loc])
              for target_loc, selected in selected_guides_ordered
              for offset, seq in selected.items()]
    return guides


def _transform_primer_product(row) -> str:
    """
    This is only necessary because Crispor returns NNNs in primer product.

    Max says: "I am masking repeats to N when sending the sequence to
    primer3. This was a major request by some labs, as they found that the
    primer3 primers were sometimes not specific at all. What I SHOULD do
    one day would be to run the primers through bwa to check their
    uniqueness, but since i'm not doing that right now, I simply mask the
    repeats, which is at least something to reduce the amount of
    nonspecific binding. "
    """

    # Only look up product from chr loc if crispor returns mysterious Ns
    if 'N' not in row['primer_product']:
        return row['primer_product']

    if row['_guide_direction'] == '+':
        guide_seq = row['guide_seq']
    else:
        guide_seq = reverse_complement(row['guide_seq'])

    if guide_seq not in row['primer_product']:
        # TODO (gdingle): use yellow highlighting of crispor to determine
        # guide_seq location. See http://crispor.tefor.net/crispor.py?ampLen=400&tm=60&batchId=fL1KMBReetZZeDh1XBkm&pamId=s29-&pam=NGG
        return 'too many repeats: ' + row['primer_product']

    logger.warning('Replacing Ns in primer product for {} guide {}'.format(
        row['target_loc'], row['guide_seq']))

    # TODO (gdingle): this is nearly the only IO in this file... do we really need it?
    primer_loc = get_primer_loc(
        row['primer_product'],
        guide_seq,
        row['guide_loc'])
    converted = conversions.chr_loc_to_seq(
        str(primer_loc),
        row['target_genome'])
    assert converted.startswith(row['primer_seq_fwd'])
    assert converted.endswith(reverse_complement(row['primer_seq_rev']))

    return converted


def from_analysis(analysis: Analysis) -> DataFrame:
    primer_selection = PrimerSelection.objects.filter(
        primer_design__guide_selection__guide_design__experiment=analysis.experiment)[0]
    sheet = from_primer_selection(primer_selection)
    return _from_analysis(analysis, sheet)


def from_custom_analysis(analysis: Analysis) -> DataFrame:
    sheet = from_experiment(analysis.experiment)[:len(analysis.fastq_data)]
    return _from_analysis(analysis, sheet)


def from_excel(file: UploadedFile) -> DataFrame:
    if file.content_type in (
        'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
            'application/vnd.ms-excel'):
        sheet = pandas.read_excel(file, sheet_name=0)
    else:
        sheet = pandas.read_csv(file)
    # TODO (gdingle): trim and validate all cells?
    return sheet


def _from_analysis(analysis: Analysis, sheet: DataFrame) -> DataFrame:

    # TODO (gdingle): does this metadata stick?
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

    # TODO (gdingle): remove _insert_fastqs when new method proved
    # sheet = _insert_fastqs(sheet, analysis.fastqs)
    if not analysis.fastq_data:
        return sheet

    sheet['fastq_fwd'] = [pair[0] for pair in analysis.fastq_data]
    sheet['fastq_rev'] = [pair[1] for pair in analysis.fastq_data]

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


def to_excel(sheet: DataFrame) -> BytesIO:
    # Workaround for saving in-memory. See:
    # https://stackoverflow.com/questions/28058563/write-to-stringio-object-using-pandas-excelwriter
    excel_file = BytesIO()
    xlw = pandas.ExcelWriter('temp.xlsx', engine='openpyxl')
    sheet.to_excel(xlw)
    xlw.book.save(excel_file)
    excel_file.seek(0)
    return excel_file


# We actually assume a larger well-plate here so we can assign more than
# 96 wells before there have been any drop-outs.
def _new_index(size=192,
               end_char='P',
               end_int=12) -> list:
    chars = [chr(i) for i in range(ord('A'), ord(end_char) + 1)]
    ints = list(range(1, end_int + 1))
    assert len(chars) * len(ints) >= size
    return [c + str(i) for c in chars for i in ints][:size]


def _new_samplesheet() -> DataFrame:
    return DataFrame(
        index=_new_index(),
        columns=[
            'target_genome',
            'target_pam',
            'target_loc',
            'target_seq',
            'guide_offset',
            'guide_loc',
            '_guide_direction',
            'guide_seq',
            'guide_pam',
            'guide_score',
            '_crispor_batch_id',
            '_crispor_pam_id',
            '_crispor_guide_id',
            'hdr_dist',
            'hdr_template',
            'hdr_rebind',
            'hdr_mutated',
            # TODO (gdingle): temp
            'hdr_guide_match',
            'primer_seq_fwd',
            'primer_seq_rev',
            # TODO (gdingle): rename to reference amplicon?
            'primer_product',
            # TODO (gdingle): ever useful?
            # 'well_name',
            's3_bucket',
            's3_prefix',
            'fastq_fwd',
            'fastq_rev',
            'report_url',
            'report_zip',
            'report_stats',
        ])

# TODO (gdingle): _insert_fastqs is deprecated pending whether
# illumina sequencer sample sheet will correctly communicate names
# of fastq files.
# TODO (gdingle): remove me when matching proven


def _insert_fastqs(sheet: DataFrame, fastqs: list) -> DataFrame:
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


def _drop_empty_report_stats(reports: list) -> Optional[Dict[str, int]]:
    """
    See https://stackoverflow.com/questions/21164910/delete-column-in-pandas-if-it-is-all-zeros
    """
    report_stats = [r.get('report_stats', {}) for r in reports]
    if not any(report_stats):
        return None
    temp_sheet = DataFrame.from_records(report_stats)
    nonzero_cols = (temp_sheet != 0).any(axis='rows')
    temp_sheet = temp_sheet.loc[:, nonzero_cols]
    return temp_sheet.to_dict(orient='records')


def _set_hdr_cols(sheet: DataFrame, hdr_seq: str, hdr_tag: str) -> DataFrame:
    sheet['hdr_dist'] = sheet.apply(
        lambda row: get_guide_cut_to_insert(
            row['target_loc'],
            row['guide_loc'],
            hdr_tag,
        ),
        axis=1,
    )
    sheet['hdr_template'] = sheet.apply(
        # TODO (gdingle): do reverse complement of hdr_seq if guide is negative strand after construction
        # TODO (gdingle): what's the max length of the template? is it the full cds?
        lambda row: hdr.HDR(
            row['target_seq'],
            hdr_seq,
            hdr_tag,
            row['hdr_dist'],
            row['_guide_direction']).template,
        axis=1,
    )
    sheet['hdr_rebind'] = sheet.apply(
        # less than 14 nucleotides* of the original protospacer remaining to be safe
        lambda row: abs(row['hdr_dist']) >= 14,
        axis=1,
    )

    # TODO (gdingle): move to hdr.py and write test cases

    # TODO (gdingle): override hdr_template when good enough
    sheet['hdr_mutated'] = sheet.apply(
        lambda row: hdr.HDR(
            row['target_seq'],
            hdr_seq,
            hdr_tag,
            row['hdr_dist'],
            row['_guide_direction']).template_mutated if row['hdr_rebind'] else '',
        axis=1)

    # TODO (gdingle): temp for comparing guide seq
    # TODO (gdingle): add in PAM?
    # sheet['hdr_guide_match'] = sheet.apply(
    #     lambda row: [
    #         hdr.HDR(
    #             row['target_seq'],
    #             hdr_seq,
    #             hdr_tag,
    #             row['hdr_dist'],
    #             row['_guide_direction']).guide_seq,
    #         row['guide_seq'] if row['_guide_direction'] == '+' else reverse_complement(
    #             row['guide_seq']),
    #     ],
    #     axis=1)

    return sheet


# TODO (gdingle): try using bolton itertools remap?
def _set_scores(
        sheet: DataFrame,
        guide_design: GuideDesign,
        # TODO (gdingle): figure out shorter type sig for guides
        guides: List[Tuple[str, str, str, str, str]]) -> DataFrame:
    # Get first score by target and offset, MIT score
    scores = dict((row['target'], row['scores'])
                  for row in guide_design.guide_data
                  if row.get('scores'))
    # Scores are stored as chr_loc -> offset -> list of numbers
    # TODO (gdingle): this is pretty ugly
    sheet['guide_score'] = [int(scores[g[0]][g[1]][0]) for g in guides]
    return sheet
