# Generated by Django 2.0.7 on 2018-07-20 20:32

import django.contrib.postgres.fields.jsonb
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0004_guidedesign_donor_data'),
    ]

    operations = [
        migrations.AddField(
            model_name='guideselection',
            name='selected_donors',
            field=django.contrib.postgres.fields.jsonb.JSONField(blank=True, default={}, null=True),
        ),
        migrations.AddField(
            model_name='guideselection',
            name='selected_guides_tagin',
            field=django.contrib.postgres.fields.jsonb.JSONField(blank=True, default={}, null=True),
        ),
    ]