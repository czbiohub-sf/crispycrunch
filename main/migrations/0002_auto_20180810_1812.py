# Generated by Django 2.1 on 2018-08-11 01:12

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0001_initial'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='guideplatelayout',
            name='guide_selection',
        ),
        migrations.RemoveField(
            model_name='primerplatelayout',
            name='primer_selection',
        ),
        migrations.AddField(
            model_name='primerdesign',
            name='primer_temp',
            field=models.IntegerField(default=60),
        ),
        migrations.DeleteModel(
            name='GuidePlateLayout',
        ),
        migrations.DeleteModel(
            name='PrimerPlateLayout',
        ),
    ]
