# Generated by Django 2.0.7 on 2018-07-20 20:23

import django.contrib.postgres.fields.jsonb
from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0003_auto_20180720_1204'),
    ]

    operations = [
        migrations.AddField(
            model_name='guidedesign',
            name='donor_data',
            field=django.contrib.postgres.fields.jsonb.JSONField(blank=True, default={}, null=True),
        ),
    ]