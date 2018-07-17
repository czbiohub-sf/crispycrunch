from django.contrib.postgres.fields import JSONField
from django.db import models

from main.platelayout import Plate384Layout, Plate96Layout
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
    # TODO (gdingle): help_text and naming convention
    name = models.CharField(max_length=40)
    researcher = models.ForeignKey(Researcher, on_delete=models.PROTECT)
    # TODO (gdingle): status field

    def __str__(self):
        return 'Experiment({}, ...)'.format(self.name)


class GuideDesign(models.Model):
    # TODO (gdingle): review guide to experiment cardinality
    experiment = models.ForeignKey(Experiment, on_delete=models.PROTECT)

    genome = models.CharField(max_length=80, choices=[
        ('hg19', 'Homo sapiens - Human - UCSC Feb. 2009 (GRCh37/hg19) + SNPs: 1000Genomes, ExaC'),
        ('todo', 'TODO: more genomes'),
    ], default='hg19')
    pam = models.CharField(max_length=80, choices=[
        ('NGG', '20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
        ('todo', 'TODO: more pams'),
    ], default='NGG')
    targets = models.CharField(max_length=65536, validators=[validate_chr_or_seq])
    # TODO (gdingle): should hdr be one or many per targets?
    # Homology Directed Repair
    hdr_seq = models.CharField(max_length=65536, validators=[validate_chr_or_seq], blank=True, null=True)

    # TODO (gdingle): put in HDR here as well?

    # TODO (gdingle): postgres array... what about relational queries?
    # from django.contrib.postgres import fields
    # targets = fields.ArrayField(
    #     models.CharField(max_length=65536, validators=[validate_seq, validate_chr]))

    # TODO (gdingle): extract guide into own model?
    guide_data = JSONField(null=True, default={}, blank=True)

    def __str__(self):
        return 'GuideDesign({}, {}, {}, ...)'.format(
            self.genome, self.pam, self.targets)


class GuideSelection(models.Model):
    guide_design = models.ForeignKey(GuideDesign, on_delete=models.PROTECT)
    # TODO (gdingle): migrate
    selected_guides = JSONField(null=True, default={}, blank=True)

    def __str__(self):
        return 'GuideSelection({}, ...)'.format(self.selected_guides)

    def order_form(self):
        # TODO (gdingle): implement IDT order form from selection
        raise


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
    # TODO: Should this always be same as pam of GuideSelection?
    pam = models.CharField(max_length=80, choices=[
        ('NGG', '20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1'),
        ('todo', 'TODO: more pams'),
    ], default='NGG')
    # TODO (gdingle): Addgene plasmid type, primers_per_guide
    # TODO (gdingle): define range better
    primer_temp = models.IntegerField(default=60)
    maximum_amplicon_length = models.IntegerField(default=400)
    # TODO (gdingle): this should actually be many pam_ids from the set of selected_guides
    pam_id = models.CharField(max_length=40, blank=True, null=True)
    # TODO (gdingle): extract data into own model?
    primer_data = JSONField(null=True, default={}, blank=True)

    def __str__(self):
        return 'PrimerDesign({}, {}, ...)'.format(self.pam, self.pam_id)


class PrimerSelection(models.Model):
    primer_design = models.ForeignKey(PrimerDesign, on_delete=models.PROTECT)
    # TODO (gdingle): extract primer into own model
    selected_primers = JSONField(null=True, default={}, blank=True)

    def __str__(self):
        return 'PrimerSelection({}, ...)'.format(self.selected_primers)

    def order_form(self):
        # TODO (gdingle): implement IDT order form from selection
        raise


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
        return Plate384Layout(self.primer_selection.selected_primers)


# TODO (gdingle): finish me... see R script
# TODO (gdingle): does this replace the experiment review step?
# TODO (gdingle): is this needed if the analysis pipeline reads straight from the DB?
class SampleSheet(models.Model):
    # TODO (gdingle): derive
    guide_selection = models.ForeignKey(
        GuideSelection, on_delete=models.PROTECT)
    primer_selection = models.ForeignKey(
        PrimerSelection, on_delete=models.PROTECT)
    # TODO (gdingle):
    # check,,Redone
