import functools

from django.contrib.postgres import fields
from django.contrib.postgres.fields import JSONField
from django.db import models

# TODO (gdingle): still needed?
# from main.platelayout import Plate96Layout
from main.validators import *

# TODO (gdingle): mark some model fields as unique or unique_together, and some as editable=False
# TODO (gdingle): review on_delete behaviors

# TODO (gdingle): temp
JASON_LI_EXAMPLE = [
    'chr2:38377365-38377424',
    'chr11:63671410-63671469',
    'chrX:153701031-153701090',
    'chr14:75657546-75657605',
    'chr14:75651664-75651723',
    'chr14:23097577-23097636',
    'chr14:23097913-23097972',
    'chr15:67521137-67521196',
    'chr15:67526436-67526495',
    'chr18:35972616-35972675',
    'chr18:35978838-35978897',
    'chr19:1275474-1275533',
    'chr19:1277214-1277273',
    'chr19:12734688-12734747',
    'chr19:12731029-12731088',
    'chr19:13774309-13774368',
    'chr19:13778035-13778094',
    'chr1:53220672-53220731',
    'chr1:53214723-53214782',
    'chr21:32612313-32612372',
    'chr21:32601953-32602012',
    'chr7:1127420-1127479',
    'chr7:997684-997743',
    'chr7:135662408-135662467',
    'chr7:135674007-135674066',
    'chr8:145052416-145052475',
    'chr8:145054080-145054139',
    'chr8:85217610-85217669',
    'chr8:85214589-85214648',
    'chr9:128160324-128160383',
    'chr9:128163437-128163496',
    'chr9:27567182-27567241',
    'chr12:55728387-55728446',
    'chr11:2377415-2377474',
    'chr12:6200379-6200438',
    'chr9:36190899-36190958',
    'chr5:176416441-176416500',
    'chr17:59693711-59693770',
    'chr5:75379392-75379451',
    'chr12:56128228-56128287',
    'chr6:29945285-29945344',
    'chr1:11258674-11258733',
    'chr2:189771810-189771869',
    'chr12:55820298-55820357',
    'chr17:39922620-39922679',
    'chr11:59615778-59615837',
    'chr3:125594895-125594954',
    'chr12:76487647-76487706',
    'chr2:26034454-26034513',
    'chr15:65869483-65869542',
    'chr19:8390321-8390380',
    'chr18:8609599-8609658',
    'chr1:153986256-153986315',
    'chr9:121193466-121193525',
    'chr10:27504201-27504260',
    'chr2:65130036-65130095',
    'chr11:66268524-66268583',
    'chr12:71754963-71755022',
    'chr6:57210498-57210557',
    'chr5:177303351-177303410',
    'chr4:13484196-13484255',
    'chr1:205775078-205775137',
    'chr8:60517056-60517115',
    'chr14:21476978-21477037',
    'chr6:146543708-146543767',
    'chr4:139454092-139454151',
    'chr12:120116667-120116726',
    'chr19:11337445-11337504',
    'chr1:220272441-220272500',
    'chr17:82698716-82698775',
    'chr16:590214-590273',
    'chr3:19950819-19950878',
    'chr12:55986864-55986923',
    'chr3:133895555-133895614',
    'chr3:128795279-128795338',
    'chr15:63189444-63189503',
    'chrX:13708616-13708675',
    'chr19:41959405-41959464',
    'chr17:5282310-5282369',
    'chr9:125202921-125202980',
    'chr9:122956944-122957003',
    'chr1:174219010-174219069',
    'chr14:24271141-24271200',
    'chr1:75786089-75786148',
    'chr1:202889114-202889173',
    'chr3:120742615-120742674',
    'chr7:151519523-151519582',
    'chr2:55050409-55050468',
    'chr3:45703309-45703368',
    'chr9:99222443-99222502',
    'chr4:89726670-89726729',
    'chr9:92032520-92032579',
    'chr14:34462221-34462280',
    'chr18:9914152-9914211',
    'chr20:58389364-58389423',
]


class BaseModel(models.Model):

    class Meta:
        abstract = True
        ordering = ['-id']

    create_time = models.DateTimeField(auto_now_add=True)
    update_time = models.DateTimeField(auto_now=True)


class Researcher(BaseModel):
    first_name = models.CharField(max_length=40)
    last_name = models.CharField(max_length=40)

    class Meta:
        unique_together = ('first_name', 'last_name')

    @property
    def full_name(self):
        return self.first_name + ' ' + self.last_name

    def __str__(self):
        return 'Researcher({}, {}, ...)'.format(
            self.first_name, self.last_name)


class Experiment(BaseModel):
    # TODO (gdingle): help_text and naming convention
    name = models.CharField(max_length=40, unique=True)
    researcher = models.ForeignKey(Researcher, on_delete=models.PROTECT)
    description = models.CharField(max_length=65536, blank=True)
    # TODO (gdingle): status field

    def __str__(self):
        return 'Experiment({}, ...)'.format(self.name)


class GuideDesign(BaseModel):
    experiment = models.ForeignKey(Experiment, on_delete=models.PROTECT)

    # TODO (gdingle): make sure this works for TagIn as well
    genome = models.CharField(max_length=80, choices=[
        ('hg19', 'Homo sapiens - Human - UCSC Feb. 2009 (GRCh37/hg19) + SNPs: 1000Genomes, ExaC'),
        ('todo', 'TODO: more genomes'),
    ], default='hg19')

    pam = models.CharField(max_length=80, choices=[
        # TODO (gdingle): what does this description all mean?
        ('NGG', '20bp-NGG - SpCas9, SpCas9-HF1, eSpCas9 1.1'),
        ('todo', 'TODO: more pams'),
    ], default='NGG')

    # TODO (gdingle): crispor has a max length of 2000 bp... validate here?
    targets = fields.ArrayField(
        models.CharField(max_length=65536, validators=[validate_chr_or_seq_or_enst_or_gene]),
        # TODO (gdingle): support FASTA with description line
        help_text='Chr location, seq, ENST, or gene. One per line.',
        # TODO (gdingle): temp default for testing
        # default=['chr7:5569176-5569415', 'chr1:11,130,540-11,130,751'],
        # default=['ENST00000330949'],
        # default=['ATL2', 'ATL3'],
        default=JASON_LI_EXAMPLE,
    )
    target_seqs = fields.ArrayField(
        models.CharField(max_length=65536, validators=[validate_seq]),
        blank=True,
        default=[],
    )

    hdr_seq = models.CharField(
        max_length=65536,
        validators=[validate_chr_or_seq_or_enst_or_gene],
        blank=True,
        help_text='Sequence for Homology Directed Repair',
        # TODO (gdingle): is this default correct? it was taken from Jason Li sample sheet example
        default='CGTGACCACATGGTCCTTCATGAGTATGTAAATGCTGCTGGGATTACAGGTGGCGGAttggaagttttgtttcaaggtccaggaagtggt')

    guide_data = JSONField(default=list, blank=True, help_text='Data returned by external service')
    donor_data = JSONField(default=list, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'GuideDesign({}, {}, {}, ...)'.format(
            self.genome, self.pam, self.targets)


class GuideSelection(BaseModel):
    guide_design = models.ForeignKey(GuideDesign, on_delete=models.PROTECT)
    selected_guides_tagin = JSONField(
        default=dict,
        blank=True,
        help_text='sgRNAs from tagin.stembio.org')

    def _validate_selected_guides(val):
        return [validate_seq(seq)
                for seqs in val.values()
                for seq in seqs.values()]

    selected_guides = JSONField(
        default=dict,
        blank=True,
        validators=[validate_num_wells, _validate_selected_guides],
        help_text='sgRNAs from crispor.tefor.net')
    # TODO (gdingle): best name: donor or HDR?
    selected_donors = JSONField(default=dict, blank=True,
                                help_text='ssDNAs from tagin.stembio.org')

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


class PrimerDesign(BaseModel):
    guide_selection = models.ForeignKey(
        GuideSelection, on_delete=models.PROTECT)
    # TODO (gdingle): any point in specifying temp?
    primer_temp = models.IntegerField(default=60)
    # TODO (gdingle): this needs to change based on HDR
    max_amplicon_length = models.IntegerField(default=400)
    primer_data = JSONField(default=list, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'PrimerDesign({}, {}, ...)'.format(self.primer_temp, self.max_amplicon_length)


class PrimerSelection(BaseModel):

    primer_design = models.ForeignKey(PrimerDesign, on_delete=models.PROTECT)

    def _validate_selected_primers(val):
        return [validate_seq(seq)
                for seqs in val.values()
                for seq in seqs]

    selected_primers = JSONField(
        default=dict,
        blank=True,
        validators=[
            functools.partial(validate_num_wells, max=96 * 2),
            _validate_selected_primers,
        ],
        help_text='Primers from crispor.tefor.net')

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


class Analysis(BaseModel):
    experiment = models.ForeignKey(
        Experiment, on_delete=models.PROTECT)
    # TODO (gdingle): default this to experiment researcher
    researcher = models.ForeignKey(
        Researcher, on_delete=models.PROTECT,
        help_text='The researcher doing the analysis')
    # TODO (gdingle): remove me on next migration
    name = models.CharField(max_length=40, blank=True)

    # TODO (gdingle): switch to czb-seqbot/fastqs/180802_M05295_0148_000000000-D49T2/?region=us-east-1&tab=overview
    # or czbiohub-seqbot/fastqs/?region=us-east-1&tab=overview
    s3_bucket = models.CharField(max_length=80,
                                 default='jasonli-bucket')
    s3_prefix = models.CharField(max_length=160,
                                 default='JasonHDR/96wp1sorted-fastq/')

    results_data = JSONField(default=dict, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'Analysis({}, {} ...)'.format(self.s3_bucket, self.s3_prefix)
