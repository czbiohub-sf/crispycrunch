# Generated by Django 2.1.1 on 2018-11-12 23:42

import django.contrib.postgres.fields
import django.contrib.postgres.fields.jsonb
import django.core.validators
import django.db.models.deletion
import functools
import main.models
import utils.validators

from django.conf import settings
from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='Analysis',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('update_time', models.DateTimeField(auto_now=True)),
                ('s3_bucket', models.CharField(default='ryan.leenay-bucket', help_text='The Amazon S3 bucket that contains the FastQ files to be analyzed', max_length=80)),
                ('s3_prefix', models.CharField(default='Greg_CXCR4_iPSC', help_text='The S3 directory that contains the FastQ files to be analyzed', max_length=160)),
                ('results_data', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default=list, help_text='Data returned by external service')),
                ('fastq_data', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default=list)),
            ],
            options={
                'ordering': ['-id'],
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Experiment',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('update_time', models.DateTimeField(auto_now=True)),
                ('name', models.CharField(max_length=40, unique=True)),
                ('description', models.CharField(blank=True, max_length=65536)),
                ('owner', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-id'],
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='GuideDesign',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('update_time', models.DateTimeField(auto_now=True)),
                ('genome', models.CharField(choices=[('hg38', 'Homo sapiens - Human - UCSC Dec. 2013 (GRCh38/hg38)')], default='hg38', max_length=80)),
                ('pam', models.CharField(choices=[('NGG', '20bp-NGG (SpCas9, SpCas9-HF1, eSpCas9, ...)')], default='NGG', help_text='Protospacer Adjacent Motif', max_length=80, verbose_name='PAM')),
                ('targets_raw', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=65536, validators=[utils.validators.validate_chr_or_seq_or_enst_or_gene]), default=['ENST00000617316,N', 'ENST00000278840,C', 'ENST00000638572,N', 'ENST00000361781,C', 'ENST00000317551,C', 'ENST00000460006,N', 'ENST00000323646,N', 'ENST00000300737,C', 'ENST00000356978,N', 'ENST00000325110,C', 'ENST00000287936,C', 'ENST00000228510,C', 'ENST00000368467,N', 'ENST00000301012,C', 'ENST00000381344,C', 'ENST00000282841,C', 'ENST00000220584,C', 'ENST00000265896,N', 'ENST00000430767,C', 'ENST00000356396,C', 'ENST00000450723,C', 'ENST00000279263,N', 'ENST00000261507,C', 'ENST00000352397,C', 'ENST00000370274,N', 'ENST00000495186,C', 'ENST00000264027,C', 'ENST00000355527,N', 'ENST00000371269,C', 'ENST00000216484,C', 'ENST00000406396,C', 'ENST00000623882,C', 'ENST00000271688,C', 'ENST00000251363,C', 'ENST00000323699,C', 'ENST00000374279,C', 'ENST00000306851,N', 'ENST00000341156,N', 'ENST00000394684,C', 'ENST00000351288,N', 'ENST00000323374,C', 'ENST00000245222,C', 'ENST00000247225,N', 'ENST00000321276,N', 'ENST00000373202,C', 'ENST00000216264,N', 'ENST00000352035,C', 'ENST00000306749,C', 'ENST00000372458,C', 'ENST00000354666,C', 'ENST00000369816,C', 'ENST00000304434,C', 'ENST00000394607,C', 'ENST00000508821,C', 'ENST00000370355,C', 'ENST00000319540,C', 'ENST00000350997,C', 'ENST00000278829,C', 'ENST00000611707,C', 'ENST00000396987,C', 'ENST00000371696,C', 'ENST00000398063,C', 'ENST00000320285,C', 'ENST00000285518,C', 'ENST00000302182,N', 'ENST00000283415,N', 'ENST00000262134,N', 'ENST00000261407,C', 'ENST00000314891,N', 'ENST00000305997,C', 'ENST00000245615,C', 'ENST00000366997,C', 'ENST00000264775,C', 'ENST00000371250,C', 'ENST00000424479,C', 'ENST00000381883,C', 'ENST00000295887,N', 'ENST00000545121,N', 'ENST00000219789,C', 'ENST00000229266,N', 'ENST00000517309,N', 'ENST00000308020,N', 'ENST00000367466,C', 'ENST00000539276,C', 'ENST00000374316,N', 'ENST00000283006,C', 'ENST00000447750,C', 'ENST00000380191,c', 'ENST00000306480,C', 'ENST00000268607,N', 'ENST00000377190,C', 'ENST00000358901,N', 'ENST00000274140,C'], help_text='Chromosome location, fasta sequence, ENST transcript ID, or\n            gene name. One per line. Append ",N" or ",C" to a line to tag at N or\n            C terminus.', size=None, validators=[utils.validators.validate_unique_set], verbose_name='Target regions')),
                ('target_locs', django.contrib.postgres.fields.ArrayField(base_field=main.models.ChrLocField(max_length=80, validators=[utils.validators.validate_chr]), size=None, validators=[utils.validators.validate_unique_set], verbose_name='Target chromosome locations')),
                ('target_seqs', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=65536, validators=[utils.validators.validate_seq]), size=None, validators=[utils.validators.validate_unique_set], verbose_name='Target sequences')),
                ('target_genes', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, max_length=40, validators=[utils.validators.validate_gene]), size=None, verbose_name='Target gene symbols')),
                ('target_tags', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(blank=True, choices=[('start_codon', 'Within 36bp after start codon (N-terminus)'), ('stop_codon', 'Within 36bp before or after stop codon (C-terminus)'), ('per_target', 'As specified per target ("N" or "C")')], max_length=40), blank=True, size=None, verbose_name='Target HDR tags')),
                ('hdr_tag', models.CharField(blank=True, choices=[('start_codon', 'Within 36bp after start codon (N-terminus)'), ('stop_codon', 'Within 36bp before or after stop codon (C-terminus)'), ('per_target', 'As specified per target ("N" or "C")')], help_text='Insert a sequence by HDR (Homology Directed Repair). Requires ENST transcript IDs.', max_length=40, verbose_name='Insert tag by HDR')),
                ('hdr_start_codon_tag_seq', models.CharField(blank=True, choices=[('ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT', 'Neon Green - ACCGAGCTCAA...'), ('GACTACAAAGACGATGACGACAAG', 'FLAG - GACTACAAAGA...'), ('GACTACAAGGACCACGACGGTGACTACAAGGACCACGACATCGACTACAAGGACGACGACGACAAG', 'XFLAG - GACTACAAGGA...'), ('GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG', 'V5 - GGTAAGCCTAT...'), ('TACCCATACGATGTTCCAGATTACGCT', 'HA - TACCCATACGA...'), ('GAACAAAAACTCATCTCAGAAGAGGATCTG', 'MYC - GAACAAAAACT...'), ('CGTGACCACATGGTCCTTCATGAGTATGTAAATGCTGCTGGGATTACAGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT', 'GFP11 - CGTGACCACAT...')], default=('ACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATGGGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGT', 'Neon Green - ACCGAGCTCAA...'), help_text='Sequence of tag to insert just after start codon', max_length=65536, validators=[utils.validators.validate_seq], verbose_name='Tag sequence for start codon')),
                ('hdr_stop_codon_tag_seq', models.CharField(blank=True, choices=[('GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG', 'Neon Green - GGTGGCGGATT...'), ('CTTGTCGTCATCGTCTTTGTAGTC', 'FLAG - CTTGTCGTCAT...'), ('CTTGTCGTCGTCGTCCTTGTAGTCGATGTCGTGGTCCTTGTAGTCACCGTCGTGGTCCTTGTAGTC', 'XFLAG - CTTGTCGTCGT...'), ('CGTAGAATCGAGACCGAGGAGAGGGTTAGGGATAGGCTTACC', 'V5 - CGTAGAATCGA...'), ('AGCGTAATCTGGAACATCGTATGGGTA', 'HA - AGCGTAATCTG...'), ('CAGATCCTCTTCTGAGATGAGTTTTTGTTC', 'MYC - CAGATCCTCTT...'), ('ACCACTTCCTGGACCTTGAAACAAAACTTCCAATCCGCCACCTGTAATCCCAGCAGCATTTACATACTCATGAAGGACCATGTGGTCACG', 'GFP11 - ACCACTTCCTG...')], default=('GGTGGCGGATTGGAAGTTTTGTTTCAAGGTCCAGGAAGTGGTACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATATGATG', 'Neon Green - GGTGGCGGATT...'), help_text='Sequence of tag to insert just before stop codon', max_length=65536, validators=[utils.validators.validate_seq], verbose_name='Tag sequence for stop codon')),
                ('guide_data', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default=list, help_text='Data returned by external service')),
                ('experiment', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='main.Experiment')),
                ('owner', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-id'],
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='GuideSelection',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('update_time', models.DateTimeField(auto_now=True)),
                ('selected_guides', django.contrib.postgres.fields.jsonb.JSONField(default=dict, help_text='Guides returned by Crispor. Filtered and ranked.', validators=[functools.partial(utils.validators.validate_num_wells, *(), **{'max': 288}), main.models.GuideSelection._validate_selected_guides])),
                ('guide_design', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='main.GuideDesign')),
                ('owner', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-id'],
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='PrimerDesign',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('update_time', models.DateTimeField(auto_now=True)),
                ('primer_temp', models.IntegerField(default=60, validators=[django.core.validators.MinValueValidator(58), django.core.validators.MaxValueValidator(62)], verbose_name='Primer melting temperature')),
                ('max_amplicon_length', models.IntegerField(default=400, help_text='amplicon = primer product', validators=[django.core.validators.MinValueValidator(200), django.core.validators.MaxValueValidator(400)], verbose_name='Maximum amplicon length')),
                ('primer_data', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default=list, help_text='Data returned by external service')),
                ('guide_selection', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='main.GuideSelection')),
                ('owner', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
            options={
                'ordering': ['-id'],
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='PrimerSelection',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('create_time', models.DateTimeField(auto_now_add=True)),
                ('update_time', models.DateTimeField(auto_now=True)),
                ('selected_primers', django.contrib.postgres.fields.jsonb.JSONField(default=dict, help_text='Primers returned by Crispor, grouped by guide, forward primer then reverse primer', validators=[functools.partial(utils.validators.validate_num_wells, *(), **{'max': 384}), main.models.PrimerSelection._validate_selected_primers])),
                ('owner', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
                ('primer_design', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='main.PrimerDesign')),
            ],
            options={
                'ordering': ['-id'],
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='analysis',
            name='experiment',
            field=models.ForeignKey(help_text='The CrispyCrunch experiment to be analyzed', on_delete=django.db.models.deletion.CASCADE, to='main.Experiment'),
        ),
        migrations.AddField(
            model_name='analysis',
            name='owner',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL),
        ),
    ]
