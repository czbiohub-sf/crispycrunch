"""
Database models for Crispycrunch app.

See SampleSheetTestCase for example data.
"""

import functools

from pandas import DataFrame

from django.contrib.postgres import fields
from django.contrib.postgres.fields import JSONField
from django.core.validators import MaxValueValidator, MinValueValidator
from django.db import models
from django.utils.text import slugify

from utils.chrloc import ChrLoc
from utils.validators import *


# TODO (gdingle): ArrayField is not showing full error messages, only first part
# "Item 1 in the array did not validate:". See:
# https://docs.djangoproject.com/en/2.1/_modules/django/contrib/postgres/forms/array/


# TODO (gdingle): temp

N_TERMINUS_EXAMPLES = [
    # TODO (gdingle): understand why all commented out are "not found"
    # 'ENST00000066544',
    'ENST00000222990',
    'ENST00000225728',
    'ENST00000250498',
    'ENST00000252992',
    'ENST00000256474',
    'ENST00000258648',
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
    'ENST00000323927',
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
    'ENST00000455511',
    # 'ENST00000529196',
    'ENST00000579978',
    'ENST00000610888',
    # 'ENST00000615911',
    'ENST00000621663',
]
C_TERMINUS_EXAMPLES = [
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
    # 'ENST00000295682',
    'ENST00000296255',
    'ENST00000299300',
    'ENST00000306467',
    'ENST00000310955',
    'ENST00000311832',
    # TODO (gdingle): bad length... only 71
    # 'ENST00000312865',
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
    'ENST00000424347',
    'ENST00000430055',
    'ENST00000456793',
    # 'ENST00000473414',
    'ENST00000481195',
    'ENST00000503026',
    'ENST00000506447',
    'ENST00000556440',
    'ENST00000602624',
]

ENST_EXAMPLE = [
    'ENST00000221801',
    'ENST00000411809',
    'ENST00000398844',
]

# TODO (gdingle): temp
JASON_LI_EXAMPLE = [
    'chr2:38377424-38377154',
    'chr11:63671469-63671209',
    'chrX:153701090-153700826',
    'chr14:75657605-75657356',
    'chr14:75651723-75651474',
    'chr14:23097577-23097833',
    'chr14:23097913-23098179',
    'chr15:67521137-67521414',
    'chr15:67526436-67526700',
    'chr18:35972616-35972885',
    'chr18:35978838-35979107',
    'chr19:1275474-1275724',
    'chr19:1277214-1277466',
    'chr19:12734747-12734478',
    'chr19:12731088-12730836',
    'chr19:13774309-13774571',
    'chr19:13778035-13778297',
    'chr1:53220731-53220476',
    'chr1:53214782-53214523',
    'chr21:32612372-32612122',
    'chr21:32602012-32601760',
    'chr7:1127479-1127213',
    'chr7:997743-997489',
    'chr7:135662408-135662673',
    'chr7:135674007-135674269',
    'chr8:145052416-145052678',
    'chr8:145054080-145054341',
    'chr8:85217669-85217405',
    'chr8:85214648-85214399',
    'chr9:128160324-128160593',
    'chr9:128163437-128163702',
    'chr9:27567241-27566981',
    'chr12:55728446-55728192',
    'chr11:2377415-2377673',
    'chr12:6200379-6200636',
    'chr9:36190899-36191152',
    'chr5:176416500-176416251',
    'chr17:59693711-59693960',
    'chr5:75379451-75379177',
    'chr12:56128228-56128494',
    'chr6:29945285-29945542',
    'chr1:11258733-11258474',
    'chr2:189771869-189771618',
    'chr12:55820298-55820562',
    'chr17:39922679-39922414',
    'chr11:59615837-59615570',
    'chr3:125594954-125594704',
    'chr12:76487706-76487443',
    'chr2:26034454-26034708',
    'chr15:65869483-65869736',
    'chr19:8390321-8390591',
    'chr18:8609599-8609871',
    'chr1:153986315-153986059',
    'chr9:121193525-121193263',
    'chr10:27504201-27504465',
    'chr2:65130095-65129826',
    'chr11:66268524-66268773',
    'chr12:71754963-71755230',
    'chr6:57210557-57210282',
    'chr5:177303410-177303154',
    'chr4:13484255-13483981',
    'chr1:205775137-205774862',
    'chr8:60517056-60517321',
    'chr14:21477037-21476777',
    'chr6:146543708-146543974',
    'chr4:139454092-139454365',
    'chr12:120116726-120116450',
    'chr19:11337504-11337255',
    'chr1:220272500-220272247',
    'chr17:82698775-82698506',
    'chr16:590214-590472',
    'chr3:19950819-19951070',
    'chr12:55986864-55987121',
    'chr3:133895614-133895363',
    'chr3:128795279-128795541',
    'chr15:63189444-63189718',
    'chrX:13708616-13708872',
    'chr19:41959464-41959204',
    'chr17:5282310-5282588',
    'chr9:125202921-125203181',
    'chr9:122956944-122957210',
    'chr1:174219010-174219265',
    'chr14:24271200-24270935',
    'chr1:75786089-75786358',
    'chr1:202889173-202888904',
    'chr3:120742674-120742415',
    'chr7:151519582-151519312',
    'chr2:55050468-55050218',
    'chr3:45703309-45703572',
    'chr9:99222443-99222712',
    'chr4:89726729-89726465',
    'chr9:92032579-92032330',
    'chr14:34462280-34462021',
    'chr18:9914152-9914417',
    'chr20:58389364-58389613',
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


class Researcher(BaseModel):
    first_name = models.CharField(max_length=40)
    last_name = models.CharField(max_length=40)

    class Meta:
        unique_together = ('first_name', 'last_name')

    @property
    def full_name(self):
        return self.first_name + ' ' + self.last_name

    def __str__(self):
        return self.full_name


class Experiment(BaseModel):
    name = models.CharField(max_length=40, unique=True)
    researcher = models.ForeignKey(Researcher, on_delete=models.CASCADE)
    description = models.CharField(max_length=65536, blank=True)

    def __str__(self):
        return self.name.title()

    @property
    def is_custom_analysis(self):
        # TODO (gdingle): change to special id=1 value when ready
        return self.name == 'No experiment -- Custom analysis'

    @property
    def short_name(self):
        return slugify(self.name)


class GuideDesign(BaseModel):
    GENOMES = [
        ('hg38', 'Homo sapiens - Human - UCSC Dec. 2013 (GRCh38/hg38)'),
        ('hg19', 'Homo sapiens - Human - UCSC Feb. 2009 (GRCh37/hg19)'),
        ('todo', 'TODO: more genomes'),
    ]
    GENOME_TO_ORGANISM = {
        'hg38': 'Human',
        'hg19': 'Human',
    }
    HDR_TAG_TERMINUSES = [
        # See cds_length below
        ('start_codon', 'Within 36bp after start codon (N-terminus)'),
        ('stop_codon', 'Within 36bp before or after stop codon (C-terminus)'),
        ('per_target', 'As indicated per target ("N" or "C")'),
    ]
    HDR_TAG_TERMINUS_TO_HDR_SEQ = {
        'start_codon': 'ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT',
        'stop_codon': 'GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG',
        # TODO (gdingle): also offer GFP?
        # CGTGACCACATGGTCCTTCATGAGTATGTAAATGCTGCTGGGATTACAGGTGGCGGAttggaagttttgtttcaaggtccaggaagtggt
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
        # default=N_TERMINUS_EXAMPLES,
        default=C_TERMINUS_EXAMPLES,
    )

    targets = fields.ArrayField(
        ChrLocField(max_length=80, validators=[validate_chr], blank=True),
        verbose_name='Target chromosome locations',
    )
    target_seqs = fields.ArrayField(
        models.CharField(max_length=65536, validators=[validate_seq]),
        verbose_name='Target sequences',
        blank=True,
        default=[],
    )
    target_tags = fields.ArrayField(
        models.CharField(
            choices=HDR_TAG_TERMINUSES,
            blank=True,
            max_length=40),
        verbose_name='Target HDR tags',
        blank=True,
        default=[],
    )
    hdr_tag = models.CharField(
        choices=HDR_TAG_TERMINUSES,
        blank=True,
        max_length=40,
        verbose_name='Insert tag by HDR',
        # TODO (gdingle): no longer GFP... switch to other name Neon something?
        help_text='Insert GFP (Green Fluorescent Protein) by HDR (Homology Directed Repair). Requires ENST transcript IDs.')

    # TODO (gdingle): custom encoder/decoder for custom dict wrapper object
    guide_data = JSONField(default=list, blank=True,
                           help_text='Data returned by external service')

    def parse_targets_raw(self) -> tuple:
        """
        Parse out optional terminuses from input such as "ENST00000621663,N"
        """
        terminus_to_tag = {
            'N': 'start_codon',
            'C': 'stop_codon',
        }
        parsed = [t.split(',') for t in self.targets_raw]
        targets_raw = [p[0].strip() for p in parsed]
        target_tags = [terminus_to_tag[p[1].strip()]
                       for p in parsed if len(p) > 1]
        if self.hdr_tag != 'per_target' and target_tags:
            raise ValueError(
                'HDR tags entered per target but also "{}". Did you mean to select "per target" HDR?'.format(
                    self.hdr_tag_verbose))
        assert not target_tags or len(target_tags) == len(targets_raw)
        return targets_raw, target_tags

    # TODO (gdingle): use in _flatten_guide_data
    def to_df(self) -> DataFrame:
        target_inputs, target_tags = self.parse_targets_raw()
        df_targets = DataFrame(data={
            'target_input': target_inputs,
            'target_loc': self.targets,
            'target_seq': self.target_seqs,
            'target_tag': target_tags or None,
            'hdr_seq': self.hdr_seq or None,
            'cds_index': self.cds_index or None,
            'cds_length': self.cds_length or None,
        })
        df_guides = DataFrame()
        for gd in self.guide_data:
            # TODO (gdingle): should this be in scraperequest?
            df_guides = df_guides.append(DataFrame(data={
                # scalars
                'target_loc': gd['target'],
                'url': gd['url'],
                'batch_id': gd['batch_id'],
                # collections
                'score': gd['scores'],
                'guide_seq': gd['guide_seqs'],
                'primer_url': gd['primer_urls'],
            }))
        return df_targets.set_index('target_loc').join(
            df_guides.set_index('target_loc'))

    def __str__(self):
        return 'GuideDesign({}, {}, {}, ...)'.format(
            self.genome, self.pam, self.targets)

    @property
    def hdr_seq(self):
        """
        This funky property will either return a scalar or a tuple.
        """
        if not self.hdr_tag:
            return None
        elif self.hdr_tag == 'per_target':
            return tuple(self.HDR_TAG_TERMINUS_TO_HDR_SEQ[hdr_tag]
                         for hdr_tag in self.target_tags)
        else:
            return self.HDR_TAG_TERMINUS_TO_HDR_SEQ[self.hdr_tag]

    @property
    def wells_per_target(self):
        return max(1, 96 // len(self.targets))

    HDR_TAG_TO_CDS_INDEX = {
        'per_target': None,
        'start_codon': 0,
        'stop_codon': -1,
    }

    @property
    def cds_index(self):
        if not self.hdr_tag:
            return None
        elif self.hdr_tag == 'per_target':
            return tuple(self.HDR_TAG_TO_CDS_INDEX[tag] for tag in self.target_tags)
        else:
            return self.HDR_TAG_TO_CDS_INDEX[self.hdr_tag]

    HDR_TAG_TO_CDS_LENGTH = {
        'per_target': None,
        'start_codon': 36,
        'stop_codon': 72,
    }

    @property
    def cds_length(self):
        if not self.hdr_tag:
            return None
        elif self.hdr_tag == 'per_target':
            return tuple(self.HDR_TAG_TO_CDS_LENGTH[tag] for tag in self.target_tags)
        return self.HDR_TAG_TO_CDS_LENGTH[self.hdr_tag]

    @property
    def hdr_tag_verbose(self):
        if not self.hdr_tag:
            return None
        elif self.hdr_tag == 'per_target':
            return set(self.HDR_TAG_TERMINUSES[tag] for tag in self.target_tags)
        return dict(self.HDR_TAG_TERMINUSES)[self.hdr_tag]


class GuideSelection(BaseModel):
    guide_design = models.ForeignKey(GuideDesign, on_delete=models.CASCADE)

    def _validate_selected_guides(val):
        return [validate_seq(seq)  # type: ignore
                for seqs in val.values()
                for seq in seqs.values()]

    selected_guides = JSONField(
        default=dict,
        blank=True,
        validators=[
            functools.partial(validate_num_wells, max=96 * 2),
            _validate_selected_guides],
        help_text='Guides returned by Crispor. Filtered and ranked.')

    def __str__(self):
        return 'GuideSelection({}, ...)'.format(self.selected_guides)

    @property
    def samplesheet(self):
        # Import here to avoid circular import
        from main.samplesheet import from_guide_selection
        return from_guide_selection(self)

    @property
    def order_form_url(self):
        return '/main/guide-selection/{}/order-form'.format(self.id)

    # TODO (gdingle): use in _flatten_guide_data
    def to_df(self) -> DataFrame:
        df = DataFrame()
        for target_loc, sgs in self.selected_guides.items():
            df = df.append(DataFrame({
                'target_loc': target_loc,
                'offset': sgs.keys(),
                'seq': sgs.values(),
            }))
        return df


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
        verbose_name='Maximum amplicon length',
        help_text='amplicon = primer product',
        default=400,
        validators=[
            # Constrain range to Biohub plausible experiments
            MinValueValidator(200),
            MaxValueValidator(400),
        ])
    primer_data = JSONField(default=list, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'PrimerDesign({}, {}, ...)'.format(self.primer_temp, self.max_amplicon_length)

    @property
    def amplicon_length(self):
        """
        Knocks down size a notch to make space for hdr_seq in primer
        """
        hdr_tag = self.guide_selection.guide_design.hdr_tag
        if hdr_tag:
            # TODO (gdingle): generalize up to len(hdr_seq) 200
            return self.max_amplicon_length - 100
        else:
            return self.max_amplicon_length


class PrimerSelection(BaseModel):

    primer_design = models.ForeignKey(PrimerDesign, on_delete=models.CASCADE)

    def _validate_selected_primers(val):
        return [validate_seq(seq[0])  # type: ignore
                for seqs in val.values()
                for seq in seqs]

    selected_primers = JSONField(
        default=dict,
        blank=True,
        validators=[
            # TODO (gdingle): check explicitly for "not found" or else different user instructions
            functools.partial(validate_num_wells, max=96 * 2),
            _validate_selected_primers,
        ],
        help_text='Primers returned by Crispor, grouped by guide, forward primer then reverse primer')

    def __str__(self):
        return 'PrimerSelection({}, ...)'.format(self.selected_primers)

    @property
    def samplesheet(self):
        # Import here to avoid circular import
        from main.samplesheet import from_primer_selection
        return from_primer_selection(self)

    @property
    def order_form_url(self):
        return '/main/primer-selection/{}/order-form'.format(self.id)

    @property
    def illumina_sheet_url(self):
        return '/main/primer-selection/{}/illumina-sheet'.format(self.id)

    @property
    def hdr_order_form_url(self):
        return '/main/primer-selection/{}/hdr-order-form'.format(self.id)

    # TODO (gdingle): use in _flatten_guide_data
    def to_df(self) -> DataFrame:
        df = DataFrame()
        for guide_id, primer_pair in self.selected_primers.items():
            target_loc, _crispor_pam_id = guide_id.split(' ')
            df = df.append(DataFrame({
                'target_loc': target_loc,
                '_crispor_pam_id': _crispor_pam_id,
                'primer_seq_fwd': primer_pair[0][0],
                'primer_seq_rev': primer_pair[1][0],
                'primer_product': primer_pair[0][1],
                # TODO (gdingle): use index of two cols?
            }, index=[guide_id]))
        return df


class Analysis(BaseModel):
    experiment = models.ForeignKey(
        Experiment, on_delete=models.CASCADE,
        help_text='The Crispycrunch experiment to be analyzed')
    # TODO (gdingle): default this to experiment researcher
    researcher = models.ForeignKey(
        Researcher, on_delete=models.CASCADE,
        help_text='The researcher doing the analysis')

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

    @property
    def is_custom(self):
        return self.experiment.is_custom_analysis

    @property
    def is_complete(self):
        return len(self.results_data) > 0

    @property
    def s3_url(self):
        return 'https://console.aws.amazon.com/s3/buckets/{}/{}/'.format(
            self.s3_bucket, self.s3_prefix)

    @property
    def s3_address(self):
        return f's3://{self.s3_bucket}/{self.s3_prefix}'
