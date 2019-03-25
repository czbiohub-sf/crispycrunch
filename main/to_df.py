"""
Provides to-dataframe functionality to app models.
"""

from pandas import Categorical, DataFrame, Series

# See https://stackoverflow.com/questions/39740632/python-type-hinting-without-cyclic-imports/43041058
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from main import models  # noqa

# TODO (gdingle): refactor constant out of scraperequest
NOT_FOUND = 'not found'


def gd_to_df(gd: 'models.GuideDesign') -> DataFrame:
    """
    Returns a flattened representation of the instance data. In particular,
    the JSON-stored data in guide_data is joined row-by-row with the array-
    stored targets data.
    """
    target_inputs, target_tags = gd.parse_targets_raw()
    tag_to_terminus = dict((v, k) for k, v in gd.TERMINUS_TO_TAG.items())

    # Keep original import order
    target_inputs = Categorical(
        target_inputs,
        categories=Series(target_inputs).unique(),
        ordered=True)

    df_targets = DataFrame(data={
        'target_input': target_inputs,
        'target_loc': gd.target_locs,
        'target_seq': gd.target_seqs,
        'target_gene': gd.target_genes,
        # Below depends on per_target
        '_hdr_tag': target_tags or gd.hdr_tag,
        '_cds_index': [gd.HDR_TAG_TO_CDS_INDEX[t] for t in target_tags] or gd.cds_index,
        '_hdr_seq': [gd.hdr_start_codon_tag_seq if t == 'start_codon' else gd.hdr_stop_codon_tag_seq
                     for t in target_tags] or gd.hdr_seq,
        'target_terminus': [tag_to_terminus[t] for t in target_tags] or None,
    })

    df_guides = DataFrame()
    for gd in gd.guide_data:
        if NOT_FOUND in gd['guide_seqs']:
            df_guides = df_guides.append(DataFrame(data={
                'target_loc': gd['target'],
                '_url': gd['url'],
                '_guide_id': gd['target'] + ' ' + NOT_FOUND,
            }, index=[gd['target']]), sort=False)
        else:
            df_guides = df_guides.append(DataFrame(data={
                # scalars
                'target_loc': gd['target'],
                '_url': gd['url'],
                '_crispor_batch_id': gd['batch_id'],
                # collections
                '_crispor_pam_id': list(gd['guide_seqs'].keys()),
                '_guide_id': [gd['target'] + ' ' + _crispor_pam_id for
                              _crispor_pam_id in gd['guide_seqs']],
                'guide_seq': list(gd['guide_seqs'].values()),
                '_scores': list(gd['scores'].values()),  # list of lists
            }))

    if not len(df_guides) or not 'guide_seq' in df_guides.columns:  # Edge case
        # TODO (gdingle): handle zero guides case better
        raise ValueError('No guides found for any targets')

    return df_targets.set_index('target_loc', drop=False).join(
        df_guides.set_index('target_loc'), how='inner')


def sg_to_df(gs: 'models.GuideSelection') -> DataFrame:
    """
    Returns a dataframe representation of selected_guides that can be joined
    easily with to_df of GuideDesign.
    """
    df = DataFrame()
    for target_loc, sgs in gs.selected_guides.items():
        df = df.append(DataFrame({
            # All we need here is an ID for filtering GuideDesign to_df
            '_guide_id': [target_loc + ' ' + _crispor_pam_id for
                          _crispor_pam_id in sgs],
            # TODO (gdingle): do we want to allow manual override of guide seq?
            # 'guide_seq_selected': list(sgs.values()),
        }))
    return df


def ps_to_df(ps: 'models.PrimerSelection') -> DataFrame:
    """
    Returns a dataframe representation of selected_primers that can be joined
    easily with to_df of GuideDesign.
    """
    df = DataFrame()
    for guide_id, primers in ps.selected_primers.items():
        df = df.append(DataFrame({
            'primer_seq_fwd': primers[0] if NOT_FOUND not in primers else NOT_FOUND,
            'primer_seq_rev': primers[1] if NOT_FOUND not in primers else NOT_FOUND,
            'primer_product': primers[2] if NOT_FOUND not in primers else NOT_FOUND,
            '_primer_adapt_seq_fwd': ps.primer_design.adapter_seq_left,
            '_primer_adapt_seq_rev': ps.primer_design.adapter_seq_right,
            '_primer_adapt_name': ps.primer_design.adapter_name,
        }, index=[guide_id]))
    return df
