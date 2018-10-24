"""
This represents the contents of an experiment on a per sample, per well basis,
whereas the Django models represent contents per plate.

For tests, see main/tests.py.
"""
import logging
import os

from io import BytesIO
from typing import Dict, Optional

import pandas

from pandas import DataFrame

from django.core.files.uploadedfile import UploadedFile

from main.models import Analysis, Experiment, GuideDesign, GuideSelection, PrimerSelection
from protospacex import get_cds_seq
from utils import conversions
from utils import hdr
from utils.chrloc import ChrLoc, get_guide_cut_to_insert, get_guide_loc, get_primer_loc

# TODO (gdingle): move to conversions
from crispresso.fastqs import reverse_complement

logger = logging.getLogger(__name__)


def from_experiment(experiment: Experiment) -> DataFrame:
    # TODO (gdingle): how to assign metadata? nothing seems to work... see stackoverflow
    return _new_samplesheet()


def from_guide_selection(guide_selection: GuideSelection) -> DataFrame:
    guide_design = guide_selection.guide_design
    sheet = from_experiment(guide_design.experiment)

    guides = _join_guide_data(guide_selection)

    # Trim sheet to available guides
    sheet = sheet[0:len(guides)]

    sheet['target_genome'] = guide_design.genome
    sheet['target_pam'] = guide_design.pam

    sheet['target_loc'] = [
        g['target_loc'] if isinstance(g['target_loc'], ChrLoc) else ChrLoc(g['target_loc'])
        for g in guides.to_records()]
    sheet['target_seq'] = list(guides['target_seq'])
    sheet['target_gene'] = list(guides['target_gene'])

    # Add original target string only if different
    sheet['target_input'] = [(g['target_input'] if g['target_input'] != str(g['target_loc']) else None)
                             for g in guides.to_records()]

    # TTCCGGCGCGCCGAGTCCTT AGG
    assert all(' ' in g['guide_seq'] for g in guides.to_records()), 'Expecting trailing PAM'
    assert all(len(g['guide_seq']) == 24 for g in guides.to_records()), 'Expecting 20bp guide'
    sheet['guide_seq'] = [g['guide_seq'].split(' ')[0] for g in guides.to_records()]
    sheet['guide_pam'] = [g['guide_seq'].split(' ')[1] for g in guides.to_records()]

    # Taken direct from Crispor. Example: "s207-" and "s76+".
    # See http://crispor.tefor.net/manual/.
    # DEFINITION: number of chars before the first char of NGG PAM
    sheet['guide_offset'] = [int(g['_crispor_pam_id'][1:-1]) for g in guides.to_records()]

    # Currently, this is relative to target strand, because that's how
    # Crispor defines it. So _guide_strand == '+' means "same strand as target".
    sheet['_guide_strand'] = [g['_crispor_pam_id'][-1] for g in guides.to_records()]

    sheet['guide_loc'] = sheet.apply(
        lambda row: get_guide_loc(
            row['target_loc'],
            row['guide_offset'],
            len(row['guide_seq']),
            row['_guide_strand']),
        axis=1,
    )

    # Take the MIT score
    sheet['guide_score'] = [g['scores'][0] for g in guides.to_records()]

    if guide_design.hdr_tag:
        sheet = _set_hdr_cols(sheet, guide_design, guides)

    sheet['_crispor_batch_id'] = list(guides['_crispor_batch_id'])
    sheet['_crispor_pam_id'] = list(guides['_crispor_pam_id'])
    # TODO (gdingle): why cannot use guide_id? see _join_guide_data
    sheet['_crispor_guide_id'] = [
        str(g['target_loc']) + ' ' + g['_crispor_pam_id']
        for g in guides.to_records()]

    return sheet


def from_primer_selection(primer_selection: PrimerSelection) -> DataFrame:
    guide_selection = primer_selection.primer_design.guide_selection
    sheet = from_guide_selection(guide_selection)

    ps_df = primer_selection.to_df()
    sheet = sheet.set_index('_crispor_guide_id').join(ps_df, how='inner')
    assert len(sheet) <= len(ps_df)

    sheet = sheet.dropna(subset=['primer_seq_fwd'])
    assert len(sheet)
    sheet.index = _new_index(size=len(sheet))

    sheet['primer_product'] = sheet.apply(_transform_primer_product, axis=1)

    if guide_selection.guide_design.hdr_tag:
        sheet = _set_hdr_primer(
            sheet,
            guide_selection.guide_design,
            primer_selection.primer_design.max_amplicon_length)

    return sheet


def _set_hdr_primer(sheet: DataFrame, guide_design: GuideDesign, max_amplicon_length: int):

    def get_primer_product(row):
        primer_product = row['primer_product']

        if ' ' in primer_product:
            # previous warning
            return primer_product

        guide_seq_aligned = hdr.HDR(
            row['target_seq'],
            row['_hdr_seq'],
            row['_hdr_tag'],
            row['hdr_dist'],
            row['_guide_strand']).guide_seq_aligned

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
            row['_hdr_seq'],
            row['_hdr_tag'],
            row['hdr_dist'],
            row['_guide_strand']).template
        assert len(before) + len(hdr_primer_product) == len(primer_product) + \
            len(row['_hdr_seq'])

        # TODO: Normalize back to direction of Crispor return?
        return before + hdr_primer_product

    # TODO (gdingle): fix root cause of issues... need more control over primer3
    def warn_hdr_primer(row) -> str:
        primer_product = row['primer_product']

        if ' ' in primer_product:
            # previous warning
            return primer_product

        plen = len(primer_product)
        if plen > max_amplicon_length:
            return f'too long, {plen}bp: {primer_product}'

        arms = primer_product.upper().split(row['_hdr_seq'])
        larm, rarm = len(arms[0]), len(arms[1])
        if min(larm, rarm) < 105:
            return 'homology arm too short, {}bp: {}'.format(min(larm, rarm), primer_product)

        return primer_product

    sheet['primer_product'] = sheet.apply(get_primer_product, axis=1)
    sheet['primer_product'] = sheet.apply(warn_hdr_primer, axis=1)
    return sheet


def _join_guide_data(guide_selection: GuideSelection) -> DataFrame:
    gd_df = guide_selection.guide_design.to_df()
    gs_df = guide_selection.to_df()

    # TODO (gdingle): why can't keep guide_id? ValueError: name already used as a name or title
    guides_df = gd_df.set_index('guide_id').join(
        gs_df.set_index('guide_id'), how='inner')
    assert len(guides_df) >= len(gs_df)

    return guides_df


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

    if row['_guide_strand'] == '+':
        guide_seq = row['guide_seq']
    else:
        guide_seq = reverse_complement(row['guide_seq'])

    if guide_seq not in row['primer_product']:
        # TODO (gdingle): use yellow highlighting of crispor to determine
        # guide_seq location. See http://crispor.tefor.net/crispor.py?ampLen=400&tm=60&batchId=fL1KMBReetZZeDh1XBkm&pamId=s29-&pam=NGG
        return 'too many repeat Ns: ' + row['primer_product']

    logger.warning('Replacing Ns in primer product for {} guide {}'.format(
        row['target_loc'], row['guide_seq']))

    # TODO (gdingle): this is nearly the only IO in this file... do we really need it?
    primer_loc = get_primer_loc(
        row['primer_product'],
        guide_seq,
        row['guide_loc'])
    # TODO (gdingle): this does not take into account strand! VERY BAD!
    converted = conversions.chr_loc_to_seq(
        str(primer_loc),
        row['target_genome'])

    if row['target_loc'].strand == '+':
        assert row['primer_product'].startswith(row['primer_seq_fwd'])
        assert row['primer_product'].endswith(reverse_complement(row['primer_seq_rev']))
    else:
        assert row['primer_product'].startswith(row['primer_seq_rev'])
        assert row['primer_product'].endswith(reverse_complement(row['primer_seq_fwd']))

    return 'Ns converted: ' + converted


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
    # TODO (gdingle): use json pandas see https://stackoverflow.com/questions/43572359/keep-column-and-row-order-when-storing-pandas-dataframe-in-json

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
            'target_input',
            'target_terminus',
            'target_gene',
            'target_loc',
            'target_seq',
            'guide_offset',
            'guide_loc',
            '_guide_strand',
            'guide_seq',
            'guide_pam',
            'guide_score',
            '_crispor_batch_id',
            '_crispor_pam_id',
            '_crispor_guide_id',
            '_hdr_tag',
            '_hdr_seq',
            'hdr_dist',
            'hdr_template',
            'hdr_mutated',
            # TODO (gdingle): we can't have these and join with ps_df
            # We only need this here for ordering... is it okay because those are last?
            # 'primer_seq_fwd',
            # 'primer_seq_rev',
            # 'primer_product',
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


def _set_hdr_cols(sheet: DataFrame, guide_design: GuideDesign, guides: DataFrame) -> DataFrame:
    hdr_tag = guide_design.hdr_tag
    if hdr_tag == 'per_target':
        sheet['_hdr_tag'] = list(guides['target_tag'])
        sheet['_hdr_seq'] = list(guides['hdr_seq'])
        sheet['target_terminus'] = list(guides['target_terminus'])
    else:
        sheet['_hdr_tag'] = hdr_tag
        sheet['_hdr_seq'] = guide_design.hdr_seq

    sheet['hdr_dist'] = sheet.apply(
        lambda row: get_guide_cut_to_insert(
            row['target_loc'],
            row['guide_loc'],
            row['_hdr_tag'],
        ),
        axis=1,
    )
    sheet['hdr_template'] = sheet.apply(
        lambda row: hdr.HDR(
            row['target_seq'],
            row['_hdr_seq'],
            row['_hdr_tag'],
            row['hdr_dist'],
            row['_guide_strand']).template,
        axis=1,
    )

    # TODO (gdingle): override hdr_template when good enough
    def mutate(row):
        """
        Most of the logic here checks for mutation or cutting on the exon/intron
        boundary, as required by https://czi.quip.com/YbAhAbOV4aXi/.

        The best response is not clear, so we return a warning. In future,
        we may filter guides further upstream, and-or mutate around junction.
        """

        # HACK ALERT! Get the CDS seq to check for mutation on exon boundary.
        # It's another instance of IO, but should be cached always.
        # TODO (gdingle): return more info from protospacex, and store throughout
        if guide_design.hdr_tag == 'per_target':
            cds_index = GuideDesign.HDR_TAG_TO_CDS_INDEX[row['_hdr_tag']]
        else:
            cds_index = guide_design.cds_index
        assert isinstance(cds_index, int)
        cds_seq = get_cds_seq(row['target_input'].split(',')[0], cds_index, -1)

        row_hdr = hdr.HDR(
            row['target_seq'],
            row['_hdr_seq'],
            row['_hdr_tag'],
            row['hdr_dist'],
            row['_guide_strand'],
            cds_seq)

        if not row_hdr.should_mutate:
            return 'not needed'

        if not row_hdr.junction:
            return row_hdr.template_mutated

        if row_hdr.cut_in_junction:
            return 'cut in intron/exon junction: ' + row_hdr.template_mutated

        # Lowercase means mutated
        if row_hdr.mutation_in_junction:
            return 'mutation in intron/exon junction: ' + row_hdr.template_mutated

        return row_hdr.template_mutated

    sheet['hdr_mutated'] = sheet.apply(mutate, axis=1)

    def check_hdr_guide_match(row):
        ghdr = hdr.HDR(
            row['target_seq'],
            row['_hdr_seq'],
            row['_hdr_tag'],
            row['hdr_dist'],
            row['_guide_strand'])
        if row['_guide_strand'] == '+':
            guide_seq = ghdr.guide_seq
        else:
            guide_seq = reverse_complement(ghdr.guide_seq)
        assert guide_seq == row['guide_seq'] + row['guide_pam']

    sheet.apply(check_hdr_guide_match, axis=1)

    return sheet
