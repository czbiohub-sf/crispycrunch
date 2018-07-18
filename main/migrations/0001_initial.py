# Generated by Django 2.0.7 on 2018-07-18 21:33

import django.contrib.postgres.fields
import django.contrib.postgres.fields.jsonb
from django.db import migrations, models
import django.db.models.deletion
import main.validators


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Experiment',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=40)),
            ],
        ),
        migrations.CreateModel(
            name='GuideDesign',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('genome', models.CharField(choices=[('hg19', 'Homo sapiens - Human - UCSC Feb. 2009 (GRCh37/hg19) + SNPs: 1000Genomes, ExaC'), ('todo', 'TODO: more genomes')], default='hg19', max_length=80)),
                ('pam', models.CharField(choices=[('NGG', '20bp-NGG - Sp Cas9, SpCas9-HF1, eSpCas9 1.1'), ('todo', 'TODO: more pams')], default='NGG', max_length=80)),
                ('targets', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=65536, validators=[main.validators.validate_chr_or_seq]), help_text='Chr location or seq, one per line', size=None)),
                ('hdr_seq', models.CharField(blank=True, max_length=65536, null=True, validators=[main.validators.validate_chr_or_seq])),
                ('guide_data', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default={}, null=True)),
                ('experiment', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.Experiment')),
            ],
        ),
        migrations.CreateModel(
            name='GuidePlateLayout',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('group_by', models.CharField(choices=[('cell_type', 'Cell Type'), ('random', 'Random'), ('todo', 'TODO: more plate groupings')], default='cell_type', max_length=40)),
                ('order_by', models.CharField(choices=[('alphabetical', 'Alphabetical'), ('todo', 'TODO: more plate orderings')], default='alphabetical', max_length=40)),
            ],
        ),
        migrations.CreateModel(
            name='GuideSelection',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('selected_guides', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default={}, null=True)),
                ('guide_design', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.GuideDesign')),
            ],
        ),
        migrations.CreateModel(
            name='PrimerDesign',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('primer_temp', models.IntegerField(default=60)),
                ('maximum_amplicon_length', models.IntegerField(default=400)),
                ('primer_data', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default={}, null=True)),
                ('guide_selection', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.GuideSelection')),
            ],
        ),
        migrations.CreateModel(
            name='PrimerPlateLayout',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('group_by', models.CharField(choices=[('cell_type', 'Cell Type'), ('random', 'Random'), ('todo', 'TODO: more plate groupings')], default='cell_type', max_length=40)),
                ('order_by', models.CharField(choices=[('alphabetical', 'Alphabetical'), ('todo', 'TODO: more plate orderings')], default='alphabetical', max_length=40)),
            ],
        ),
        migrations.CreateModel(
            name='PrimerSelection',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('selected_primers', django.contrib.postgres.fields.jsonb.JSONField(blank=True, default={}, null=True)),
                ('primer_design', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.PrimerDesign')),
            ],
        ),
        migrations.CreateModel(
            name='Researcher',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('first_name', models.CharField(max_length=40)),
                ('last_name', models.CharField(max_length=40)),
            ],
        ),
        migrations.CreateModel(
            name='SampleSheet',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('guide_selection', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.GuideSelection')),
                ('primer_selection', models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.PrimerSelection')),
            ],
        ),
        migrations.AddField(
            model_name='primerplatelayout',
            name='primer_selection',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.PrimerSelection'),
        ),
        migrations.AddField(
            model_name='guideplatelayout',
            name='guide_selection',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.GuideSelection'),
        ),
        migrations.AddField(
            model_name='experiment',
            name='researcher',
            field=models.ForeignKey(on_delete=django.db.models.deletion.PROTECT, to='main.Researcher'),
        ),
    ]
