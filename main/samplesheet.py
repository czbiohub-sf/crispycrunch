"""
This represents all the contents of an experiment on a per sample, per well
basis, whereas the Django models represent contents per plate.

Most of the code simply transforms values to useful representations for display.

For tests, see main/tests.py.
"""
import logging
import os

from io import BytesIO
from typing import Callable, Dict, Optional

import pandas

from pandas import DataFrame

from django.conf import settings
from django.core.files.uploadedfile import UploadedFile

from main.models import Analysis, GuideDesign, GuideSelection, PrimerSelection
from protospacex import get_cds_codon_at, get_cds_seq, get_ultramer_seq
from utils import conversions
from utils import hdr
from utils import primerchecks
from utils import manuscore
from utils.chrloc import *

from crispresso.fastqs import reverse_complement

logger = logging.getLogger(__name__)

NOT_FOUND = 'not found'

# TODO (gdingle): change defaults in module
hdr.HDR.guide_seq_aligned_length = 27
hdr.HDR.use_cfd_score = True
hdr.HDR.stop_mutating_at_first_success = False
# TODO (gdingle): increase perf by caching
hdr.HDR.mutate_all_permutations = True


def from_guide_selection(guide_selection: GuideSelection) -> DataFrame:
    guide_design = guide_selection.guide_design

    sheet = _join_guide_data(guide_selection)

    sheet['_target_genome'] = guide_design.genome
    sheet['_target_pam'] = guide_design.pam

    sheet['target_loc'] = [g['target_loc']
                           if isinstance(g['target_loc'], ChrLoc)
                           else ChrLoc(g['target_loc'])
                           for g in sheet.to_records()
                           ]

    sheet = _set_guide_cols(sheet)

    if guide_design.is_hdr:
        sheet = _set_hdr_cols(sheet, guide_design, sheet)

    return sheet


def _set_guide_cols(sheet: DataFrame) -> DataFrame:

    # TTCCGGCGCGCCGAGTCCTT AGG
    assert all(' ' in g['guide_seq'] for g in sheet.to_records()
               if pandas.notna(g['guide_seq'])), 'Expecting trailing PAM'
    assert all(len(g['guide_seq']) == 24 for g in sheet.to_records()
               if pandas.notna(g['guide_seq'])), 'Expecting 20bp guide'

    def _apply(l: Callable):
        """Set empty string for all missing guide cols for friendly display"""
        return sheet.apply(
            lambda g: l(g) if pandas.notna(g['guide_seq']) else '',
            axis=1)

    sheet['guide_pam'] = _apply(lambda g: g['guide_seq'].split(' ')[1])
    sheet['guide_seq'] = _apply(lambda g: g['guide_seq'].split(' ')[0])

    # Taken direct from Crispor. Example: "s207-" and "s76+".
    # See http://crispor.tefor.net/manual/.
    # DEFINITION: number of chars before the first char of NGG PAM
    sheet['guide_offset'] = _apply(lambda g: int(g['_crispor_pam_id'][1:-1]))

    # Crispor returns guide strand relative to target, NOT genome.
    sheet['_guide_strand_same'] = _apply(
        lambda g: g['_crispor_pam_id'][-1] == '+')

    sheet['guide_loc'] = sheet.apply(
        lambda row: get_guide_loc(
            row['target_loc'],
            row['guide_offset'],
            len(row['guide_seq']),
            row['_guide_strand_same'])
        if row['guide_seq'] else '',
        axis=1,
    )

    # Take the Doench 2016 score
    sheet['guide_score'] = _apply(lambda g: int(g['_scores'][1]))

    return sheet


def _join_guide_data(guide_selection: GuideSelection) -> DataFrame:
    gd_df = guide_selection.guide_design.to_df()
    if not len(gd_df):
        return gd_df
    dupes = gd_df['_guide_id'][gd_df['_guide_id'].duplicated()]
    assert not len(dupes), dupes

    gs_df = guide_selection.to_df()
    dupes = gs_df['_guide_id'][gs_df['_guide_id'].duplicated()]
    assert not len(dupes), dupes

    guides_df = gd_df.set_index('_guide_id').join(
        gs_df.set_index('_guide_id'), how='inner')
    assert len(guides_df) >= len(gs_df), (len(guides_df), len(gs_df))

    return guides_df


def _well_positions(size=96 * 12,
                    end_char='ÿ',
                    end_int=12) -> list:
    """
    end_char and end_int determine the shape of the plate together.

    Size is the number of wells returned. We actually assume a larger than 96
    well-plate here so we can have drop-outs.
    """
    chars = [chr(i) for i in range(ord('A'), ord(end_char) + 1)]
    ints = list(range(1, end_int + 1))
    assert len(chars) * len(ints) >= size
    return [c + str(i) for c in chars for i in ints][:size]


def _set_hdr_cols(sheet: DataFrame, guide_design: GuideDesign, guides: DataFrame) -> DataFrame:
    sheet['target_terminus'] = list(guides['target_terminus'])

    sheet['hdr_dist'] = sheet.apply(
        lambda row: get_guide_cut_to_insert(
            row['target_loc'],
            row['guide_loc'],
            row['_hdr_tag'],
        ) if row['guide_seq'] else '',
        axis=1,
    )

    # HACK ALERT! Get the CDS seq to check for mutation on exon boundary.
    # It's another instance of IO, but should be cached always, first in sqlite
    # then in local mem.
    # TODO (gdingle): return more info from protospacex, and store throughout
    sheet['_seq_cds'] = sheet.apply(
        lambda row: get_cds_seq(row['target_input'], row['_cds_index'], -1),
        axis=1)
    sheet['_seq_codon_at'] = sheet.apply(
        lambda row: get_cds_codon_at(
            row['target_input'],
            row['_cds_index'],
            len(row['target_seq'])
        ),
        axis=1)

    sheet['_hdr_insert_at'] = sheet.apply(
        lambda row: get_insert(row['target_loc'], row['_hdr_tag']),
        axis=1,
    )
    sheet['hdr_inserted'] = sheet.apply(
        lambda row:
        '' if not row['guide_seq'] else
        _get_hdr_row(row).inserted,
        axis=1,
    )

    sheet['hdr_score'] = sheet.apply(
        lambda row:
        '' if not row['guide_seq'] else
        round(manuscore.manu_score(row['guide_score'], row['hdr_dist']) * 100),
        axis=1,
    )

    # TODO (gdingle): override hdr_inserted when good enough
    def mutate(row):
        """
        Most of the logic here checks for mutation or cutting on the exon/intron
        boundary, as required by https://czi.quip.com/YbAhAbOV4aXi/.

        The best response is not clear, so we return a warning. In future,
        we may filter guides further upstream, and-or mutate around junction.
        """
        if not row['guide_seq']:
            return ''

        row_hdr = _get_hdr_row(row)

        if not row_hdr.should_mutate:
            return 'not needed'

        if not row_hdr.junction:
            return row_hdr.inserted_mutated

        if row_hdr.cut_in_junction:
            return 'cut in intron/exon junction: ' + row_hdr.inserted_mutated

        # Lowercase means mutated
        if row_hdr.mutation_in_junction:
            return 'mutation in intron/exon junction: ' + row_hdr.inserted_mutated

        return row_hdr.inserted_mutated

    sheet['hdr_mutated'] = sheet.apply(mutate, axis=1)

    def check_hdr_guide_match(row):
        if not row['guide_seq']:
            return
        row_hdr = _get_hdr_row(row)
        if row['_guide_strand_same']:
            guide_seq = row_hdr.guide_seq
        else:
            guide_seq = reverse_complement(row_hdr.guide_seq)
        assert guide_seq == row['guide_seq'] + row['guide_pam']

    sheet.apply(check_hdr_guide_match, axis=1)

    def set_ultramer(row):
        if not row['guide_seq']:
            return ''

        # HACK ALERT! Get the ultramer from the ENST.
        # It's another instance of IO, but should be cached always.
        # TODO (gdingle): return more info from protospacex, and store throughout
        try:
            ultramer_seq, codon_at = get_ultramer_seq(
                row['target_input'],
                row['_cds_index'],
                110,
            )
        except ValueError:
            return NOT_FOUND

        uhdr = hdr.HDR(
            ultramer_seq,
            row['_hdr_seq'],
            row['_hdr_tag'],
            row['hdr_dist'],
            row['_guide_strand_same'],
            row['_seq_cds'],
            codon_at)
        try:
            ultramer_mutated = uhdr.inserted_mutated
        except Exception:
            # TODO (gdingle): figure out why!
            return 'error in ultramer: please contact ' + settings.ADMIN_EMAIL

        # TODO (gdingle): non-IDT ultramers?
        # assert len(ultramer_mutated) <= 200, '200bp is max for IDT ultramer'
        assert len(ultramer_mutated) >= 150, '150bp is min for IDT ultramer'

        if not row['_guide_strand_same']:
            return reverse_complement(ultramer_mutated)
        else:
            return ultramer_mutated

    sheet['_hdr_ultramer'] = sheet.apply(set_ultramer, axis=1)

    return sheet


def _get_hdr_row(row) -> hdr.HDR:
    return hdr.HDR(
        row['target_seq'],
        row['_hdr_seq'],
        row['_hdr_tag'],
        row['hdr_dist'],
        row['_guide_strand_same'],
        row['_seq_cds'],
        row['_seq_codon_at']
    )


def from_primer_selection(primer_selection: PrimerSelection) -> DataFrame:
    guide_selection = primer_selection.primer_design.guide_selection
    sheet = from_guide_selection(guide_selection)

    ps_df = primer_selection.to_df()
    sheet = sheet.join(ps_df, how='left')
    assert len(sheet) >= len(ps_df)

    sheet['primer_product'] = sheet.apply(_transform_primer_product, axis=1)

    if guide_selection.guide_design.is_hdr:
        sheet = _set_hdr_primer(
            sheet,
            guide_selection.guide_design,
            primer_selection.primer_design.max_amplicon_length)

    sheet.insert(0, 'well_pos', _well_positions(size=len(sheet)))
    # TODO (gdingle): is this wanted?
    # sheet.insert(1, 'well_num', range(1, len(sheet) + 1))

    return sheet


def _set_hdr_primer(sheet: DataFrame, guide_design: GuideDesign, max_amplicon_length: int):

    def get_primer_product(row):
        primer_product = row['primer_product']

        if ' ' in primer_product:
            # previous warning
            return primer_product

        guide_seq_aligned = _get_hdr_row(row).guide_seq_aligned.upper()

        # Crispor returns primer products by strand. Normalize to positive strand.
        if row['target_loc'].strand == '-':
            primer_product = reverse_complement(primer_product)

        guide_offset = primer_product.find(guide_seq_aligned)

        if guide_offset == -1:
            logger.warning('Could not find guide {} in primer {} for target {}'.format(
                row['guide_seq'], primer_product, row['target_loc']))
            return 'guide not found: ' + primer_product

        start = guide_offset % 3

        # TODO (gdingle): HACK ALERT!!! Because the target codon seq can appear
        # in frame but outside the CCDS, the insert is misidentified. We set a buffer here to avoid the worst.
        # The buffer length is a multiple of 3 less than the min homology len,
        # and larger than the largest observed misidentification.

        # TODO (gdingle): offset by row['_seq_codon_at']
        start += 90

        before, primer_product_aligned = \
            primer_product[:start], primer_product[start:]
        assert before + primer_product_aligned == primer_product

        phdr = hdr.HDR(
            primer_product_aligned,
            row['_hdr_seq'],
            row['_hdr_tag'],
            row['hdr_dist'],
            row['_guide_strand_same'],
            row['_seq_cds'],
            # TODO (gdingle): make this work with more thought
            # row['_seq_codon_at']
        )

        try:
            hdr_primer_product = phdr.inserted_mutated if phdr.should_mutate else phdr.inserted
        except AssertionError:
            return 'error in HDR, no insert: ' + primer_product

        assert len(before) + len(hdr_primer_product) == len(primer_product) + \
            len(row['_hdr_seq'])
        return before + hdr_primer_product

    def warn_hdr_primer(row) -> str:
        """
        These should happen rarely when using Crispor modified for HDR.
        """
        primer_product = row['primer_product']

        if ' ' in primer_product:
            # previous warning
            return primer_product

        if row['_hdr_seq'] not in primer_product.upper():
            return 'no HDR insertion: ' + primer_product

        plen = len(primer_product)
        if plen > max_amplicon_length:
            return f'too long, {plen}bp: {primer_product}'

        arms = primer_product.upper().split(row['_hdr_seq'])
        larm, rarm = len(arms[0]), len(arms[1])
        if min(larm, rarm) < 105:
            return 'homology arm too short, {}bp: {}'.format(min(larm, rarm), primer_product)

        return primer_product

    def warn_primer_self_bind(row) -> DataFrame:
        primer_product = row['primer_product']

        if ' ' in primer_product:
            # previous warning
            return primer_product

        if primerchecks.is_self_binding(row['primer_seq_fwd'], row['primer_seq_rev']):
            return 'self binding: ' + primer_product

        if primerchecks.is_self_binding_with_adapters(row['primer_seq_fwd'], row['primer_seq_rev']):
            return 'binds to adapters: ' + primer_product

        return primer_product

    sheet['primer_product'] = sheet.apply(get_primer_product, axis=1)
    sheet['primer_product'] = sheet.apply(warn_hdr_primer, axis=1)
    sheet['primer_product'] = sheet.apply(warn_primer_self_bind, axis=1)

    return sheet


# TODO (gdingle): should not be needed anymore after crispor mods
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
    if not row['guide_seq']:
        return ' '  # one space as "warning"

    if not isinstance(row['primer_product'], str):
        return NOT_FOUND

    # Only look up product from chr loc if crispor returns mysterious Ns
    if 'N' not in row['primer_product']:
        return row['primer_product']

    if row['_guide_strand_same']:
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
    # TODO (gdingle): this does not take into account strand!
    converted = conversions.chr_loc_to_seq(
        str(primer_loc),
        row['_target_genome'])

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
    # TODO (gdingle): re-test me
    sheet = DataFrame()
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
    # TODO (gdingle): does this do anything?
    # sheet.analysis_id = analysis.id
    # sheet.analysis_create_time = analysis.create_time
    # sheet.analyst_name = analysis.owner.username

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
