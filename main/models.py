"""
Database models for CrispyCrunch app.

See SampleSheetTestCase for example data.
"""

import functools

from pandas import DataFrame

from django.conf import settings
from django.contrib.postgres import fields
from django.contrib.postgres.fields import JSONField
from django.core.validators import MaxLengthValidator, MaxValueValidator, MinLengthValidator, MinValueValidator
from django.db import models
from django.utils.functional import cached_property
from django.utils.text import slugify

from main import to_df
from utils.chrloc import ChrLoc
from utils.validators import *

# TODO (gdingle): refactor constant out of scraperequest
NOT_FOUND = 'not found'

# TODO (gdingle): ArrayField is not showing full error messages, only first part
# "Item 1 in the array did not validate:". See:
# https://docs.djangoproject.com/en/2.1/_modules/django/contrib/postgres/forms/array/

VAR_TERMINUS_EXAMPLES = [
    # ENST00000447866,C length of cds only 36',
    'ENST00000617316,N',
    'ENST00000278840,C',
    'ENST00000638572,N',
    'ENST00000361781,C',
    'ENST00000317551,C',
    'ENST00000460006,N',
    'ENST00000323646,N',
    'ENST00000300737,C',
    'ENST00000356978,N',
    'ENST00000325110,C',
    'ENST00000287936,C',
    'ENST00000228510,C',
    'ENST00000368467,N',
    'ENST00000301012,C',
    'ENST00000381344,C',
    'ENST00000282841,C',
    'ENST00000220584,C',
    'ENST00000265896,N',
    'ENST00000430767,C',
    'ENST00000356396,C',
    'ENST00000450723,C',
    'ENST00000279263,N',
    'ENST00000261507,C',
    'ENST00000352397,C',
    'ENST00000370274,N',
    'ENST00000495186,C',
    'ENST00000264027,C',
    'ENST00000355527,N',
    'ENST00000371269,C',
    'ENST00000216484,C',
    'ENST00000406396,C',
    'ENST00000623882,C',
    'ENST00000271688,C',
    'ENST00000251363,C',
    'ENST00000323699,C',
    'ENST00000374279,C',
    'ENST00000306851,N',
    'ENST00000341156,N',
    'ENST00000394684,C',
    'ENST00000351288,N',
    'ENST00000323374,C',
    'ENST00000245222,C',
    'ENST00000247225,N',
    'ENST00000321276,N',
    'ENST00000373202,C',
    'ENST00000216264,N',
    'ENST00000352035,C',
    'ENST00000306749,C',
    'ENST00000372458,C',
    'ENST00000354666,C',
    'ENST00000369816,C',
    'ENST00000304434,C',
    'ENST00000394607,C',
    'ENST00000508821,C',
    'ENST00000370355,C',
    'ENST00000319540,C',
    'ENST00000350997,C',
    'ENST00000278829,C',
    'ENST00000611707,C',
    'ENST00000396987,C',
    'ENST00000371696,C',
    'ENST00000398063,C',
    'ENST00000320285,C',
    'ENST00000285518,C',
    'ENST00000302182,N',
    'ENST00000283415,N',
    'ENST00000262134,N',
    'ENST00000261407,C',
    'ENST00000314891,N',
    'ENST00000305997,C',
    'ENST00000245615,C',
    'ENST00000366997,C',
    'ENST00000264775,C',
    'ENST00000371250,C',
    'ENST00000424479,C',
    'ENST00000381883,C',
    'ENST00000295887,N',
    'ENST00000545121,N',
    'ENST00000219789,C',
    'ENST00000229266,N',
    'ENST00000517309,N',
    'ENST00000308020,N',
    'ENST00000367466,C',
    'ENST00000539276,C',
    'ENST00000374316,N',
    'ENST00000283006,C',
    'ENST00000447750,C',
    'ENST00000380191,c',
    'ENST00000306480,C',
    'ENST00000268607,N',
    'ENST00000377190,C',
    'ENST00000358901,N',
    'ENST00000274140,C',
]


N_TERMINUS_EXAMPLES = [
    # TODO (gdingle): understand why all commented out are "not found"
    # 'ENST00000066544',
    # TODO (gdingle): why error during processing?
    # 'ENST00000222990',
    'ENST00000225728',
    'ENST00000250498',
    'ENST00000252992',
    'ENST00000256474',
    # TODO (gdingle): why does get_ultramer_seq not work for this?
    # 'ENST00000258648',
    'ENST00000264935',
    'ENST00000264982',
    'ENST00000268711',
    'ENST00000269392',
    'ENST00000275603',
    'ENST00000278935',
    'ENST00000283122',
    'ENST00000286788',
    'ENST00000292035',
    'ENST00000293777',
    'ENST00000295688',
    'ENST00000295770',
    'ENST00000309439',
    'ENST00000313525',
    'ENST00000321685',
    # TODO (gdingle): why does get_ultramer_seq not work for this?
    # 'ENST00000323927',
    'ENST00000324817',
    'ENST00000331821',
    'ENST00000339839',
    'ENST00000342382',
    'ENST00000345519',
    # 'ENST00000361564',
    'ENST00000367607',
    'ENST00000370277',
    'ENST00000394903',
    'ENST00000431606',
    # TODO (gdingle): why does get_ultramer_seq not work for this?
    # 'ENST00000455511',
    # 'ENST00000529196',
    'ENST00000579978',
    'ENST00000610888',
    # 'ENST00000615911',
    'ENST00000621663',
]
C_TERMINUS_EXAMPLES = [
    # TODO (gdingle): primer error: too many perfect matches
    'ENST00000221138',
    'ENST00000227618',
    'ENST00000237380',
    'ENST00000237530',
    'ENST00000251871',
    'ENST00000254950',
    'ENST00000255764',
    'ENST00000257287',
    'ENST00000258091',
    'ENST00000261819',
    'ENST00000263205',
    'ENST00000265350',
    'ENST00000267935',
    'ENST00000282892',
    'ENST00000283195',
    'ENST00000287598',
    # TODO (gdingle): bad length... only 68
    'ENST00000295682',
    'ENST00000296255',
    'ENST00000299300',
    'ENST00000306467',
    'ENST00000310955',
    'ENST00000311832',
    # TODO (gdingle): bad length... only 71
    'ENST00000312865',
    'ENST00000315368',
    'ENST00000325542',
    'ENST00000341068',
    'ENST00000345046',
    'ENST00000352035',
    'ENST00000354910',
    'ENST00000355801',
    'ENST00000356221',
    'ENST00000366542',
    'ENST00000366999',
    'ENST00000371485',
    'ENST00000372457',
    'ENST00000373842',
    'ENST00000374080',
    'ENST00000374206',
    'ENST00000376227',
    'ENST00000376300',
    'ENST00000378230',
    'ENST00000394128',
    'ENST00000394440',
    'ENST00000394886',
    'ENST00000397527',
    'ENST00000397786',
    # TODO (gdingle): this is causing error in hdr guide_seq... bad codon?
    'ENST00000424347',
    'ENST00000430055',
    'ENST00000456793',
    'ENST00000473414',
    'ENST00000481195',
    'ENST00000503026',
    'ENST00000506447',
    'ENST00000556440',
    'ENST00000602624',
]

# First 4 only
RYAN_LEENAY_EXAMPLE = [
    'chr2:136114360-136114419',
    'chr2:136115613-136115672',
    'chr2:136116738-136116797',
    'chr2:136117544-136117603',
]


class BaseModel(models.Model):

    class Meta:
        abstract = True
        ordering = ['-id']

    create_time = models.DateTimeField(auto_now_add=True)
    update_time = models.DateTimeField(auto_now=True)
    owner = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE)


class ChrLocField(models.CharField):
    """
    Custom field for structuring 'chr7:5569177-5569415' strings.
    See https://docs.djangoproject.com/en/2.1/howto/custom-model-fields/#converting-values-to-python-objects
    """

    # TODO (gdingle): why isn't higher level value_to_string called when this is?
    # see https://docs.djangoproject.com/en/2.1/_modules/django/contrib/postgres/fields/array/
    def get_db_prep_value(self, value, connection, prepared):
        return str(value)

    def from_db_value(self, item, expression, connection, context):
        return ChrLoc(item)


# TODO (gdingle): join with user accounts, one to one, or remove
# class Researcher(BaseModel):
#     first_name = models.CharField(max_length=40)
#     last_name = models.CharField(max_length=40)

#     class Meta:
#         unique_together = ('first_name', 'last_name')

#     @cached_property
#     def full_name(self):
#         return self.first_name + ' ' + self.last_name

#     def __str__(self):
#         return self.full_name


class Experiment(BaseModel):
    name = models.CharField(
        max_length=40,
        help_text='Identifying name of the experiment'
    )
    description = models.CharField(
        max_length=65536,
        blank=True,
        help_text='Friendly description for other people',
    )
    is_hdr = models.BooleanField(
        default=True,
        verbose_name='HDR experiment?',
        help_text='Will guides be used with HDR (Homology Directed Repair)?',
    )

    def __str__(self):
        return self.name.title()

    @cached_property
    def is_custom_analysis(self):
        # TODO (gdingle): change to special id=1 value when ready
        return self.name == 'No experiment -- Custom analysis'

    @cached_property
    def short_name(self):
        return slugify(self.name)

    def get_current_step(self):
        """
        Return the rel URL to the current step to be completed.
        """
        try:
            # TODO (gdingle): DRY and extract paths somewhere with success_url
            obj = self
            view = '/main/experiment/{id}/guide-design/'

            obj = GuideDesign.objects.filter(experiment=self.id)[0]
            view = '/main/guide-design/{id}/guide-selection/'

            obj = GuideSelection.objects.filter(
                guide_design__experiment=self.id)[0]
            view = '/main/guide-selection/{id}/primer-design/'

            obj = PrimerDesign.objects.filter(
                guide_selection__guide_design__experiment=self.id)[0]
            view = '/main/primer-design/{id}/primer-selection/'

            obj = PrimerSelection.objects.filter(
                primer_design__guide_selection__guide_design__experiment=self.id)[0]
            view = '/main/primer-selection/{id}/experiment-summary/'
        except (IndexError):
            pass
        return view.format(id=obj.id)


class GuideDesign(BaseModel):
    GENOMES = [
        ('hg38', 'Homo sapiens - Human - UCSC Dec. 2013 (GRCh38/hg38)'),
        # ('hg19', 'Homo sapiens - Human - UCSC Feb. 2009 (GRCh37/hg19)'),
        # ('todo', 'TODO: more genomes'),
    ]
    GENOME_TO_ORGANISM = {
        'hg38': 'Human',
        'hg19': 'Human',
    }
    HDR_TAG_TO_CDS_LENGTH = {
        'per_target': None,
        'start_codon': 96,
        'stop_codon': 96,
    }
    HDR_TAG_TERMINUSES = [
        # See cds_length below
        # TODO (gdingle): need both start_codon and start_codon2?
        # ('start_codon', 'Within {} after start codon (N-terminus)'.format(self.HDR_TAG_TO_CDS_LENGTH['start_codon'])),
        ('start_codon', 'Within {} before or after start codon (N-terminus)'.format(
            HDR_TAG_TO_CDS_LENGTH['start_codon'])),
        ('stop_codon', 'Within {} before or after stop codon (C-terminus)'.format(
            HDR_TAG_TO_CDS_LENGTH['stop_codon'])),
        ('per_target', 'Within {} before or after, terminus specified per target ("N" or "C")'.format(
            HDR_TAG_TO_CDS_LENGTH['start_codon'])),
    ]
    # TODO (gdingle): need to verify when each of these are appropriate
    HDR_TAG_TERMINUS_TO_HDR_SEQ = {
        'start_codon': [
            # TODO (gdingle): what is first called?
            ('ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT',
             'Neon Green - ACCGAGCTCAA...'),
            ('GACTACAAAGACGATGACGACAAG', 'FLAG - GACTACAAAGA...'),
            ('GACTACAAGGACCACGACGGTGACTACAAGGACCACGACATCGACTACAAGGACGACGACGACAAG', 'XFLAG - GACTACAAGGA...'),
            ('GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG', 'V5 - GGTAAGCCTAT...'),
            ('TACCCATACGATGTTCCAGATTACGCT', 'HA - TACCCATACGA...'),
            ('GAACAAAAACTCATCTCAGAAGAGGATCTG', 'MYC - GAACAAAAACT...'),
            ('CGTGACCACATGGTCCTTCATGAGTATGTAAATGCTGCTGGGATTACAGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT',
             'GFP11 - CGTGACCACAT...'),
        ],
        'stop_codon': [
            # TODO (gdingle): what is first called?
            ('GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG',
             'Neon Green - GGTGGCGGATT...'),
            # TODO (gdingle): these are just reverse complement of above. Should that always be true?
            # See https://github.com/czbiohub/Tagin_web/blob/9fb6fe4ee86f0b99db60545472d0d884bd174071/CRISPRtag/CRISPRtag/CRISPRtag_py3.py#L960
            ('CTTGTCGTCATCGTCTTTGTAGTC', 'FLAG - CTTGTCGTCAT...'),
            ('CTTGTCGTCGTCGTCCTTGTAGTCGATGTCGTGGTCCTTGTAGTCACCGTCGTGGTCCTTGTAGTC', 'XFLAG - CTTGTCGTCGT...'),
            ('CGTAGAATCGAGACCGAGGAGAGGGTTAGGGATAGGCTTACC', 'V5 - CGTAGAATCGA...'),
            ('AGCGTAATCTGGAACATCGTATGGGTA', 'HA - AGCGTAATCTG...'),
            ('CAGATCCTCTTCTGAGATGAGTTTTTGTTC', 'MYC - CAGATCCTCTT...'),
            ('ACCACTTCCTGGACCTTGAAACAAAACTTCCAATCCGCCACCTGTAATCCCAGCAGCATTTACATACTCATGAAGGACCATGTGGTCACG',
             'GFP11 - ACCACTTCCTG...'),
        ],
    }
    TERMINUS_TO_TAG = {
        'N': 'start_codon',
        'C': 'stop_codon',
    }

    experiment = models.ForeignKey(Experiment, on_delete=models.CASCADE)

    # TODO (gdingle): enforce match of ENST transcript and genome
    genome = models.CharField(max_length=80, choices=GENOMES, default='hg38')

    pam = models.CharField(
        verbose_name='PAM',
        help_text='Protospacer Adjacent Motif',
        max_length=80,
        choices=[('NGG', '20bp-NGG (SpCas9, SpCas9-HF1, eSpCas9, ...)'), ],
        default='NGG')

    # TODO (gdingle): enforce must be ENST for HDR,
    # explain ENST is always converted to forward strand
    targets_raw = fields.ArrayField(
        models.CharField(max_length=65536, validators=[validate_chr_or_seq_or_enst_or_gene]),
        validators=[validate_unique_set],
        # TODO (gdingle): support FASTA with description line
        verbose_name='Target regions',
        help_text="""Chromosome location, fasta sequence, ENST transcript ID, or
            gene name. One per line. Append ",N" or ",C" to a line to tag at N or
            C terminus.""",
        # TODO (gdingle): temp default for testing
        # default=['chr7:5569176-5569415', 'chr1:11,130,540-11,130,751'],
        # default=['ENST00000330949'],
        # default=['ATL2', 'ATL3'],
        # default=JASON_LI_EXAMPLE,
        # default=RYAN_LEENAY_EXAMPLE,
        # default=ENST_EXAMPLE,
        default=VAR_TERMINUS_EXAMPLES,
        # default=N_TERMINUS_EXAMPLES,
        # default=C_TERMINUS_EXAMPLES,
    )

    target_locs = fields.ArrayField(
        ChrLocField(max_length=80, validators=[validate_chr]),
        validators=[validate_unique_set],
        verbose_name='Target chromosome locations',
    )
    target_seqs = fields.ArrayField(
        models.CharField(max_length=65536, validators=[validate_seq]),
        validators=[validate_unique_set],
        verbose_name='Target sequences',
    )
    target_genes = fields.ArrayField(
        models.CharField(max_length=40, validators=[validate_gene], blank=True),
        # Genes are not necessarily unique set in an experiment
        verbose_name='Target gene symbols',
    )
    target_tags = fields.ArrayField(
        models.CharField(
            choices=HDR_TAG_TERMINUSES,
            blank=True,
            max_length=40),
        verbose_name='Target HDR tags',
        blank=True,
    )

    guides_per_target = models.IntegerField(
        help_text='The top N number of guides per target to select.',
        default=1,
        validators=[
            MinValueValidator(1),
            MaxValueValidator(96),
        ])

    @property
    def is_hdr(self):
        # TODO (gdingle): simplify when legacy data gone
        return self.hdr_tag is not None or self.experiment.is_hdr

    # TODO (gdingle): rename to hdr_tag_terminus?
    hdr_tag = models.CharField(
        choices=HDR_TAG_TERMINUSES,
        # Should be blank to avoid confusion in non-HDR experiments
        default=None,
        max_length=40,
        verbose_name='HDR tag terminus',
        help_text='Where to insert the tag in each gene',
        # TODO (gdingle): how to remove this in case of GuideDesignForm2?
        blank=True,
        null=True,
    )
    hdr_start_codon_tag_seq = models.CharField(
        # TODO (gdingle): need to verify when each of these are appropriate
        # choices=HDR_TAG_TERMINUS_TO_HDR_SEQ['start_codon'],
        default=HDR_TAG_TERMINUS_TO_HDR_SEQ['start_codon'][0][0],
        max_length=65536,
        validators=[
            validate_seq,
            MinLengthValidator(3),
            MaxLengthValidator(1000),
        ],
        verbose_name='Tag sequence for start codon',
        help_text='Sequence of tag to insert just after start codon',
        # TODO (gdingle): how to remove this in case of GuideDesignForm2?
        blank=True,
        null=True,
    )
    # TODO (gdingle): is it actually needed to show this? can we determine by hdr_start_codon_tag_seq?
    hdr_stop_codon_tag_seq = models.CharField(
        # TODO (gdingle): need to verify when each of these are appropriate
        # choices=HDR_TAG_TERMINUS_TO_HDR_SEQ['stop_codon'],
        default=HDR_TAG_TERMINUS_TO_HDR_SEQ['stop_codon'][0][0],
        # TODO (gdingle): move to HDR_TAG_TERMINUS_TO_HDR_SEQ
        # > 3xHA tag- P2A-mCherry used by TagIn team: GGATCCTACCCCTACGACGTGCCCGACTACGCCGCCGGAGCATACCCATACGATGTTCCAGATTACGCTGCTGGCGCATACCCATACGACGTACCAGATTACGCTGGATCCGGCGCAACAAACTTCTCTCTGCTGAAACAAGCCGGAGATGTCGAAGAGAATCCTGGACCGATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG
        # > Glycine linker-mCherry used by TagIn team: GGTGGCAGCGGTGGCAGCATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAG
        max_length=65536,
        validators=[
            validate_seq,
            # See https://www.idtdna.com/pages/support/faqs/when-designing-donor-dna-for-use-in-homology-directed-repair-(hdr)-what-are-the-optimal-lengths-of-the-left-and-right-homology-arms-and-what-is-the-maximum-size-of-sequence-that-can-be-efficiently-inserted-in-mammalian-cells
            MinLengthValidator(3),
            MaxLengthValidator(1000),
        ],
        verbose_name='Tag sequence for stop codon',
        help_text='Sequence of tag to insert just before stop codon',
        # TODO (gdingle): how to remove this in case of GuideDesignForm2?
        blank=True,
        null=True,
    )

    # TODO (gdingle): custom encoder/decoder for custom dict wrapper object
    guide_data = JSONField(default=list, blank=True,
                           help_text='Data returned by external service')

    def parse_targets_raw(self) -> tuple:
        """
        Parse out optional terminuses from input such as "ENST00000621663,N"
        """
        parsed = [t.split(',') for t in self.targets_raw]
        targets_raw = [p[0].strip() for p in parsed]
        target_tags = [self.TERMINUS_TO_TAG[p[1].strip().upper()]
                       for p in parsed if len(p) > 1]

        if self.hdr_tag != 'per_target' and target_tags:
            raise ValueError(
                'HDR tags entered per target but also "{}". Did you mean to select "per target" HDR?'.format(
                    self.hdr_tag_verbose or '--------'))

        if self.hdr_tag == 'per_target' and not target_tags:
            raise ValueError('You must specify "N" or "C" for each target')

        assert not target_tags or len(target_tags) == len(targets_raw)
        return targets_raw, target_tags

    def to_df(self) -> DataFrame:
        return to_df.gd_to_df(self)

    def __str__(self):
        return 'GuideDesign({}, {}, {}, ...)'.format(
            self.genome, self.pam, self.target_locs)

    @cached_property
    def hdr_seq(self):
        """
        This funky property will either return a scalar or a tuple.
        """
        if not self.hdr_tag:
            return None
        elif self.hdr_tag == 'per_target':
            return (self.hdr_start_codon_tag_seq, self.hdr_stop_codon_tag_seq)
        else:
            return self.hdr_start_codon_tag_seq if self.hdr_tag == 'start_codon' else self.hdr_stop_codon_tag_seq

    @cached_property
    def hdr_seq_name(self):
        """
        For friendly display.
        """
        start_codon_choices = dict(self.HDR_TAG_TERMINUS_TO_HDR_SEQ['start_codon'])
        stop_codon_choices = dict(self.HDR_TAG_TERMINUS_TO_HDR_SEQ['stop_codon'])
        if not self.hdr_tag:
            return None
        else:
            choices = set([
                start_codon_choices.get(self.hdr_start_codon_tag_seq, 'UNKNOWN')
                .split('-')[0].strip(),
                stop_codon_choices.get(self.hdr_stop_codon_tag_seq, 'UNKNOWN')
                .split('-')[0].strip(),
            ])
            if len(choices) > 1:
                return tuple(choices)
            else:
                return choices.pop()

    HDR_TAG_TO_CDS_INDEX = {
        'per_target': None,
        'start_codon': 0,
        'stop_codon': -1,
    }

    @cached_property
    def pre_filter(self):
        """
        The number of guides to check for primer existance *before* saving in guide_data.
        """
        if self.hdr_tag:
            # Don't prefilter because primer design is too different with hdr_dist in crispor
            return 0
        else:
            # 5 based on safe-harbor experiment
            return self.guides_per_target * 5

    @cached_property
    def cds_index(self):
        if not self.hdr_tag:
            return None
        elif self.hdr_tag == 'per_target':
            return tuple(self.HDR_TAG_TO_CDS_INDEX[tag] for tag in self.target_tags)
        else:
            return self.HDR_TAG_TO_CDS_INDEX[self.hdr_tag]

    @cached_property
    def cds_length(self):
        if not self.hdr_tag:
            return None
        elif self.hdr_tag == 'per_target':
            return tuple(self.HDR_TAG_TO_CDS_LENGTH[tag] for tag in self.target_tags)
        return self.HDR_TAG_TO_CDS_LENGTH[self.hdr_tag]

    @cached_property
    def hdr_tag_verbose(self):
        if not self.hdr_tag:
            return None
        return dict(self.HDR_TAG_TERMINUSES)[self.hdr_tag]

    @cached_property
    def crispor_urls(self):
        return dict((gd['target'], gd['url'])
                    for gd in self.guide_data
                    if gd.get('url'))


class GuideSelection(BaseModel):
    guide_design = models.ForeignKey(GuideDesign, on_delete=models.CASCADE)

    def _validate_selected_guides(val):
        return [validate_seq(seq)  # type: ignore
                for seqs in val.values()
                for seq in seqs.values()
                # Allow special string through
                if seq != NOT_FOUND]

    selected_guides = JSONField(
        default=dict,
        validators=[
            # See guides_per_target
            functools.partial(validate_num_wells, max=96 * 3),
            _validate_selected_guides],
        help_text='Guides returned by Crispor. Filtered and ranked.')

    def __str__(self):
        return 'GuideSelection({}, ...)'.format(self.selected_guides)

    @cached_property
    def samplesheet(self):
        # Import here to avoid circular import
        from main.samplesheet import from_guide_selection
        return from_guide_selection(self)

    @cached_property
    def order_form_url(self):
        return '/main/guide-selection/{}/order-form'.format(self.id)

    def to_df(self) -> DataFrame:
        return to_df.sg_to_df(self)


class PrimerDesign(BaseModel):
    guide_selection = models.ForeignKey(
        GuideSelection, on_delete=models.CASCADE)
    primer_temp = models.IntegerField(
        verbose_name='Primer melting temperature',
        default=60,
        validators=[
            MinValueValidator(58),
            MaxValueValidator(62),
        ])

    # TODO (gdingle): crispor has a preset list of values... mirror?
    max_amplicon_length = models.IntegerField(
        # TODO (gdingle): this value is always set to 400 in HDR experiments because of
        # Crispor customizaitons... how to show that in UI?
        verbose_name='Maximum amplicon length',
        help_text='Amplicon = primer product. Length after HDR insertion.',
        default=400,
        validators=[
            # Constrain range to Biohub plausible experiments
            MinValueValidator(200),
            MaxValueValidator(400),
        ])
    primer_data = JSONField(default=list, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'PrimerDesign({}, {}, ...)'.format(self.primer_temp, self.max_amplicon_length)

    @cached_property
    def amplicon_length(self):
        """
        Knocks down size a notch to make space for hdr_seq in primer
        # TODO (gdingle): this is no longer needed with mods to crispor
        # need to indicate this somehow in UI
        """
        hdr_tag = self.guide_selection.guide_design.hdr_tag
        if hdr_tag:
            # TODO (gdingle): generalize up to len(hdr_seq) 200
            return self.max_amplicon_length - 100
        else:
            return self.max_amplicon_length

    @cached_property
    def crispor_urls(self):
        return dict(
            (p['target'], p['url'] + '#ontargetPcr')
            for p in self.primer_data)


class PrimerSelection(BaseModel):

    primer_design = models.ForeignKey(PrimerDesign, on_delete=models.CASCADE)

    def _validate_selected_primers(val):
        return [validate_seq(seq[0])  # type: ignore
                for seqs in val.values()
                for seq in seqs
                if seqs != NOT_FOUND]

    selected_primers = JSONField(
        default=dict,
        validators=[
            # TODO (gdingle): temporary up for comparison with other method
            functools.partial(validate_num_wells, max=96 * 12),
            _validate_selected_primers,
        ],
        help_text='Primers returned by Crispor, grouped by guide, forward primer then reverse primer')

    def __str__(self):
        return 'PrimerSelection({}, ...)'.format(self.selected_primers)

    @cached_property
    def samplesheet(self):
        # Import here to avoid circular import
        from main.samplesheet import from_primer_selection
        return from_primer_selection(self)

    @cached_property
    def order_form_url(self):
        return '/main/primer-selection/{}/order-form'.format(self.id)

    @cached_property
    def illumina_sheet_url(self):
        return '/main/primer-selection/{}/illumina-sheet'.format(self.id)

    @cached_property
    def hdr_order_form_url(self):
        return '/main/primer-selection/{}/hdr-order-form'.format(self.id)

    def to_df(self) -> DataFrame:
        return to_df.ps_to_df(self)


class Analysis(BaseModel):
    # TODO (gdingle): how to filter only to user's own experiments?
    experiment = models.ForeignKey(
        Experiment, on_delete=models.CASCADE,
        help_text='The CrispyCrunch experiment to be analyzed')

    # TODO (gdingle): switch to czb-seqbot/fastqs/180802_M05295_0148_000000000-D49T2/?region=us-east-1&tab=overview
    # or czbiohub-seqbot/fastqs/?region=us-east-1&tab=overview
    s3_bucket = models.CharField(max_length=80,
                                 # default='jasonli-bucket',
                                 default='ryan.leenay-bucket',
                                 help_text='The Amazon S3 bucket that contains the FastQ files to be analyzed'
                                 )
    s3_prefix = models.CharField(max_length=160,
                                 # default='JasonHDR/96wp1sorted-fastq/'
                                 default='Greg_CXCR4_iPSC',
                                 help_text='The S3 directory that contains the FastQ files to be analyzed'
                                 )

    results_data = JSONField(default=list, blank=True,
                             help_text='Data returned by external service')
    fastq_data = JSONField(default=list, blank=True)

    def __str__(self):
        # return 'Analysis({}, {} ...)'.format(self.s3_bucket, self.s3_prefix)
        return '{} samples of {}'.format(
            len(self.results_data), self.experiment)

    @cached_property
    def is_custom(self):
        return self.experiment.is_custom_analysis

    @cached_property
    def is_complete(self):
        return len(self.results_data) > 0

    @cached_property
    def s3_url(self):
        return 'https://console.aws.amazon.com/s3/buckets/{}/{}/'.format(
            self.s3_bucket, self.s3_prefix)

    @cached_property
    def s3_address(self):
        return f's3://{self.s3_bucket}/{self.s3_prefix}'
