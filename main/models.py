# TODO (gdingle): use non-jsonb JSONField so we preserve key order
# see https://github.com/dmkoch/django-jsonfield
from django.contrib.postgres import fields
from django.contrib.postgres.fields import JSONField
from django.db import models

from main.platelayout import Plate96Layout
from main.validators import *


# TODO (gdingle): mark some model fields as unique or unique_together, and some as editable=False
# TODO (gdingle): truncate database and remove NULLs and null=True
# TODO (gdingle): review on_delete behaviors

# TODO (gdingle): add to all models
# created_time = models.DateTimeField(auto_now_add=True)
# updated_time = models.DateTimeField(auto_now=True)

class Researcher(models.Model):
    first_name = models.CharField(max_length=40)
    last_name = models.CharField(max_length=40)

    def __str__(self):
        return 'Researcher({}, {}, ...)'.format(
            self.first_name, self.last_name)


class Experiment(models.Model):
    class Meta:
        ordering = ['-id']
    # TODO (gdingle): help_text and naming convention
    name = models.CharField(max_length=40)
    researcher = models.ForeignKey(Researcher, on_delete=models.PROTECT)
    description = models.CharField(max_length=65536, blank=True, null=True)
    # TODO (gdingle): status field
    #
    # TODO (gdingle): ordering
    # see https://docs.djangoproject.com/en/2.0/ref/models/options/

    def __str__(self):
        return 'Experiment({}, ...)'.format(self.name)


class GuideDesign(models.Model):
    # TODO (gdingle): review guide to experiment cardinality
    experiment = models.ForeignKey(Experiment, on_delete=models.PROTECT)

    # TODO (gdingle): make sure this works for TagIn as well
    genome = models.CharField(max_length=80, choices=[
        ('hg19', 'Homo sapiens - Human - UCSC Feb. 2009 (GRCh37/hg19) + SNPs: 1000Genomes, ExaC'),
        ('todo', 'TODO: more genomes'),
    ], default='hg19')
    pam = models.CharField(max_length=80, choices=[
        ('NGG', '20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
        ('todo', 'TODO: more pams'),
    ], default='NGG')
    targets = fields.ArrayField(
        # TODO (gdingle): crispor has a max length of 2000 bp... is that a problem?
        models.CharField(max_length=65536, validators=[validate_chr_or_seq_or_enst]),
        help_text='Chr location or seq or ENST, one per line. Tag-in experiments require ENST.',
        # TODO (gdingle): temp default for testing
        default=['chr7:5569176-5569415', 'chr1:11,130,540-11,130,751'],
        # default=['ENST00000330949'],
    )
    # TODO (gdingle): do we even want to convert now?
    # target_fastas = fields.ArrayField(
    #     models.CharField(max_length=65536, validators=[validate_chr_or_seq_or_enst]),
    # )
    # Homology Directed Repair # TODO (gdingle): custom seqence?
    # hdr_seq = models.CharField(max_length=65536, validators=[validate_chr_or_seq_or_enst], blank=True, null=True)
    tag_in = models.CharField(max_length=40, blank=True, null=True, choices=(
        ('FLAG', 'FLAG'),
        ('3XFLAG', '3XFLAG'),
        ('V5', 'V5'),
        ('HA', 'HA'),
        ('MYC', 'MYC'),
        ('TODO', 'TODO: tag used by Manu group'),
    ))

    # TODO (gdingle): extract guide into own model?
    guide_data = JSONField(null=True, default={}, blank=True)
    donor_data = JSONField(null=True, default={}, blank=True)

    def __str__(self):
        return 'GuideDesign({}, {}, {}, ...)'.format(
            self.genome, self.pam, self.targets)


class GuideSelection(models.Model):
    guide_design = models.ForeignKey(GuideDesign, on_delete=models.PROTECT)
    selected_guides_tagin = JSONField(null=True, default={}, blank=True,
                                      help_text='sgRNAs from tagin.stembio.org')
    selected_donors = JSONField(null=True, default={}, blank=True,
                                help_text='ssDNAs from tagin.stembio.org')
    selected_guides = JSONField(null=True, default={}, blank=True,
                                help_text='sgRNAs from crispor.tefor.net')
    # TODO (gdingle): temp for debuggin

    def __str__(self):
        return 'GuideSelection({}, ...)'.format(self.selected_guides)


class GuidePlateLayout(models.Model):
    guide_selection = models.ForeignKey(
        GuideSelection, on_delete=models.PROTECT)
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


class PrimerDesign(models.Model):
    guide_selection = models.ForeignKey(
        GuideSelection, on_delete=models.PROTECT)
    # TODO (gdingle): Addgene plasmid type
    # TODO (gdingle): define range better
    primer_temp = models.IntegerField(default=60)
    max_amplicon_length = models.IntegerField(default=400)
    # TODO (gdingle): extract data into own model?
    primer_data = JSONField(null=True, default={}, blank=True)

    def __str__(self):
        return 'PrimerDesign({}, {}, ...)'.format(self.primer_temp, self.max_amplicon_length)


class PrimerSelection(models.Model):
    primer_design = models.ForeignKey(PrimerDesign, on_delete=models.PROTECT)
    # TODO (gdingle): extract primer into own model
    selected_primers = JSONField(null=True, default={}, blank=True,
                                 help_text='Primers from crispor.tefor.net')

    def __str__(self):
        return 'PrimerSelection({}, ...)'.format(self.selected_primers)


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


class Analysis(models.Model):
    experiment = models.ForeignKey(
        Experiment, on_delete=models.PROTECT)
    researcher = models.ForeignKey(
        Researcher, on_delete=models.PROTECT)
    name = models.CharField(max_length=40)
    # TODO (gdingle): how to upload? Or refer to s3? validate file type?
    pcr_file = models.FileField(blank=True, null=True)
    read_length = models.IntegerField(default=250,
                                      help_text='What is the read length used in the experiment?')
    fragment_length = models.IntegerField(default=300,
                                          help_text='What is average fragment length (before adapters)?')
    stdev_fragment_length = models.IntegerField(default=45,
        help_text='What is the standard deviation of the fragment legnths? If unknown, assume 10% of the average.')  # noqa
