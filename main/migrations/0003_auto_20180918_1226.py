# Generated by Django 2.1.1 on 2018-09-18 19:26

import django.contrib.postgres.fields
from django.db import migrations, models
import main.validators


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0002_auto_20180917_1640'),
    ]

    operations = [
        migrations.AlterField(
            model_name='guidedesign',
            name='targets',
            field=django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=65536, validators=[main.validators.validate_chr_or_seq_or_enst_or_gene]), default=['chr2:136114360-136114419', 'chr2:136115613-136115672', 'chr2:136116738-136116797', 'chr2:136117544-136117603'], help_text='Chr location, seq, ENST, or gene. One per line. For reverse strand, write chr location right-to-left.', size=None),
        ),
    ]