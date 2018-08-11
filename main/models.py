# TODO (gdingle): use non-jsonb JSONField so we preserve key order
# see https://github.com/dmkoch/django-jsonfield
from django.contrib.postgres import fields
from django.contrib.postgres.fields import JSONField
from django.db import models

from main.platelayout import Plate96Layout
from main.validators import *


# TODO (gdingle): mark some model fields as unique or unique_together, and some as editable=False
# TODO (gdingle): review on_delete behaviors

class BaseModel(models.Model):

    class Meta:
        abstract = True

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
    class Meta:
        ordering = ['-id']
    # TODO (gdingle): help_text and naming convention
    name = models.CharField(max_length=40, unique=True)
    researcher = models.ForeignKey(Researcher, on_delete=models.PROTECT)
    description = models.CharField(max_length=65536, blank=True)
    # TODO (gdingle): status field
    #
    # TODO (gdingle): ordering
    # see https://docs.djangoproject.com/en/2.0/ref/models/options/

    def __str__(self):
        return 'Experiment({}, ...)'.format(self.name)


class GuideDesign(BaseModel):
    # TODO (gdingle): review guide to experiment cardinality
    experiment = models.ForeignKey(Experiment, on_delete=models.PROTECT)

    # TODO (gdingle): make sure this works for TagIn as well
    genome = models.CharField(max_length=80, choices=[
        ('hg19', 'Homo sapiens - Human - UCSC Feb. 2009 (GRCh37/hg19) + SNPs: 1000Genomes, ExaC'),
        ('todo', 'TODO: more genomes'),
    ], default='hg19')
    # TODO (gdingle): for pysam, we are actually using hg18!!!

    pam = models.CharField(max_length=80, choices=[
        # TODO (gdingle): what does this description all mean?
        ('NGG', '20bp-NGG - SpCas9, SpCas9-HF1, eSpCas9 1.1'),
        ('todo', 'TODO: more pams'),
    ], default='NGG')
    # TODO (gdingle): normalize to samtools region chr loc string
    targets = fields.ArrayField(
        # TODO (gdingle): crispor has a max length of 2000 bp... is that a problem?
        models.CharField(max_length=65536, validators=[validate_chr_or_seq_or_enst]),
        help_text='Chr location or seq or ENST, one per line. Tag-in experiments require ENST.',
        # TODO (gdingle): temp default for testing
        default=['chr7:5569176-5569415', 'chr1:11,130,540-11,130,751'],
        # default=['ENST00000330949'],
    )
    # TODO (gdingle): need to use pysam here
    # TODO (gdingle): do we even want to convert now?
    # target_fastas = fields.ArrayField(
    #     models.CharField(max_length=65536, validators=[validate_chr_or_seq_or_enst]),
    # )

    hdr_seq = models.CharField(
        max_length=65536,
        validators=[validate_chr_or_seq_or_enst],
        blank=True,
        help_text='Sequence for Homology Directed Repair',
        # TODO (gdingle): is this default correct? it was taken from Jason Li sample sheet example
        default='CGTGACCACATGGTCCTTCATGAGTATGTAAATGCTGCTGGGATTACAGGTGGCGGAttggaagttttgtttcaaggtccaggaagtggt')

    # TODO (gdingle): is this useful?
    # tag_in = models.CharField(max_length=40, blank=True, choices=(
    #     ('FLAG', 'FLAG'),
    #     ('3XFLAG', '3XFLAG'),
    #     ('V5', 'V5'),
    #     ('HA', 'HA'),
    #     ('MYC', 'MYC'),
    #     ('TODO', 'TODO: tag used by Manu group'),
    # ))

    guide_data = JSONField(default=dict, blank=True, help_text='Data returned by external service')
    donor_data = JSONField(default=dict, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'GuideDesign({}, {}, {}, ...)'.format(
            self.genome, self.pam, self.targets)


class GuideSelection(BaseModel):
    guide_design = models.ForeignKey(GuideDesign, on_delete=models.PROTECT)
    selected_guides_tagin = JSONField(default=dict, blank=True,
                                      help_text='sgRNAs from tagin.stembio.org')
    selected_guides = JSONField(default=dict, blank=True,
                                help_text='sgRNAs from crispor.tefor.net')
    # TODO (gdingle): best name: donor or HDR?
    selected_donors = JSONField(default=dict, blank=True,
                                help_text='ssDNAs from tagin.stembio.org')

    def __str__(self):
        return 'GuideSelection({}, ...)'.format(self.selected_guides)


# TODO (gdingle): remove me completely
class GuidePlateLayout(BaseModel):
    guide_selection = models.ForeignKey(
        GuideSelection,
        on_delete=models.PROTECT)
    # TODO (gdingle): add name when multiple plates
    # name = models.CharField(max_length=40)
    group_by = models.CharField(max_length=40, choices=[
        ('cell_type', 'Cell Type'),
        ('random', 'Random'),
        ('todo', 'TODO: more plate groupings'),
    ], default='cell_type')
    order_by = models.CharField(max_length=40, choices=[
        ('alphabetical', 'Alphabetical'),
        ('todo', 'TODO: more plate orderings'),
    ], default='alphabetical')

    def __str__(self):
        return 'GuidePlateLayout({}, {}, ...)'.format(
            self.group_by, self.order_by)

    @property
    def layout(self):
        return Plate96Layout(self.guide_selection.selected_guides)


class PrimerDesign(BaseModel):
    guide_selection = models.ForeignKey(
        GuideSelection, on_delete=models.PROTECT)
    # TODO (gdingle): Addgene plasmid type?
    # TODO (gdingle): any point in specifying temp?
    # primer_temp = models.IntegerField(default=60)
    # TODO (gdingle): this needs to change based on HDR
    max_amplicon_length = models.IntegerField(default=400)
    primer_data = JSONField(default=dict, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'PrimerDesign({}, {}, ...)'.format(self.primer_temp, self.max_amplicon_length)


class PrimerSelection(BaseModel):
    primer_design = models.ForeignKey(PrimerDesign, on_delete=models.PROTECT)
    selected_primers = JSONField(default=dict, blank=True,
                                 help_text='Primers from crispor.tefor.net')

    def __str__(self):
        return 'PrimerSelection({}, ...)'.format(self.selected_primers)

# TODO (gdingle): remove me completely


class PrimerPlateLayout(models.Model):
    primer_selection = models.ForeignKey(
        PrimerSelection, on_delete=models.PROTECT)
    # TODO (gdingle): add name when multiple plates
    # name = models.CharField(max_length=40)
    # TODO (gdingle): plate_size
    # plate_size = models.IntegerField()
    group_by = models.CharField(max_length=40, choices=[
        ('cell_type', 'Cell Type'),
        ('random', 'Random'),
        ('todo', 'TODO: more plate groupings'),
    ], default='cell_type')
    order_by = models.CharField(max_length=40, choices=[
        ('alphabetical', 'Alphabetical'),
        ('todo', 'TODO: more plate orderings'),
    ], default='alphabetical')

    def __str__(self):
        return 'PrimerPlateLayout({}, {}, ...)'.format(
            self.group_by, self.order_by)

    @property
    def layout(self):
        return Plate96Layout(self.primer_selection.selected_primers)


class Analysis(BaseModel):
    experiment = models.ForeignKey(
        Experiment, on_delete=models.PROTECT)
    # TODO (gdingle): default this to experiment
    researcher = models.ForeignKey(
        Researcher, on_delete=models.PROTECT,
        help_text='The researcher doing the analysis')
    # TODO (gdingle): remove me on next migration
    name = models.CharField(max_length=40)

    s3_bucket = models.CharField(max_length=80,
                                 default='jasonli-bucket')
    s3_prefix = models.CharField(max_length=160,
                                 default='JasonHDR/96wp1sorted-fastq/')

    results_data = JSONField(default=dict, blank=True, help_text='Data returned by external service')

    def __str__(self):
        return 'Analysis({}, {} ...)'.format(self.s3_bucket, self.s3_prefix)

    def get_selected_guides(self):
        return GuideSelection.objects.get(guide_design__experiment=self.experiment).selected_guides

    def get_selected_donors(self):
        return GuideSelection.objects.get(guide_design__experiment=self.experiment).selected_donors
