# Generated by Django 2.0.7 on 2018-07-10 23:32

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('main', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='GuideSet',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('genome', models.CharField(choices=[('homo_sapiens', 'Home sapiens - USCS (GRCh37/hg19)'), ('todo', 'TODO: more genomes')], max_length=80)),
                ('PAM', models.CharField(choices=[('Cas9', '20bp-NGG - Sp Cas9'), ('todo', 'TODO: more PAMs')], max_length=80)),
                ('targets', models.CharField(max_length=65536)),
            ],
        ),
    ]
