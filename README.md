# CrispyCrunch
Web app for CRISPR experiment setup and analysis. Fill in a 96-well plate and analyze it as a whole.

## Background

This app was developed by Greg Dingle starting in July 2018 for the genome engineering group in the CZ Biohub. The intention was to streamline the CRISPR experimentation process of researchers in the Biohub and provide a permanent record of experiments and results.

## Services

CrispyCrunch calls out to several bioinformatics web services.

* [Crispor](http://crispor.tefor.net/)
* [Crispresso2](http://crispresso.pinellolab.partners.org)
* [Togows](http://togows.org)
* [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks)
* [Ensemble](http://rest.ensembl.org/)

Currently, CrispyCrunch calls out to mirrors for Crispor and Crispresso:
* [Crispor mirror](http://ec2-34-219-237-20.us-west-2.compute.amazonaws.com/crispor.py)
* [Crispresso mirror](http://ec2-52-12-22-81.us-west-2.compute.amazonaws.com/)

*NOTE:* The Crispor mirror is modified to support HDR primer design. See [this PR](https://github.com/maximilianh/crisporWebsite/pull/21).

## Tech stack

CrispyCrunch is built with the following technologies and frameworks.

* Postgres
* Python3.6+
* Django
* Bootstrap

# Installation

Clone the repo. You must have czbiohub github access.

```git clone git@github.com:czbiohub/crispycrunch.git```

Install python dependencies

```pip install -r requirements.txt```

Start postgres. The exact command will depend on how you installed it.

```brew services start postgresql```

Configure `DATABASES` in `settings.py` if needed.

Initialize the database.

```docker-compose exec web python manage.py migrate```

Create a superuser for admin.

```python manage.py createsuperuser```

NOTE: a superuser is automatically created by `eb deploy`.

# Usage

Go to http://localhost:8000/. Create a new experiment and follow the subsequent steps.

# Deployment

Currently, CrispyCrunch uses Amazon's Elasticbeanstalk for deployment and hosting. EB manages a web server (EC2) and a database (RDS) in one deployment environment (prod). Open the EB control panel at https://us-west-2.console.aws.amazon.com/elasticbeanstalk/home?region=us-west-2#/environment/dashboard?applicationName=crispycrunch&environmentId=e-vnvbedub4n.

Install EB CLI.

```pip install awsebcli```

Commit any changes you wish to deploy.

```git commit -a -m'fancy stuff'```

Deploy to Amazon Elasticbeanstalk.

```eb deploy```

Wait for deployment to finish, then view the live site at http://crispycrunch.ds.czbiohub.org/.

For more info, see https://realpython.com/deploying-a-django-app-and-postgresql-to-aws-elastic-beanstalk/.

To download fastq files stored in s3, you must configure
```
AWS_ACCESS_KEY_ID
AWS_SECRET_ACCESS_KEY
```
See https://docs.aws.amazon.com/elasticbeanstalk/latest/dg/environments-cfg-softwaresettings.html and https://boto3.amazonaws.com/v1/documentation/api/latest/guide/configuration.html#environment-variable-configuration.

# Admin

To inspect the state of the app, use the built-in Django admin interface at http://localhost:8000/admin/main/.

To inspect the source of fastq files, open https://console.aws.amazon.com/s3/buckets/czb-seqbot/.

To read the error logs in prod:

```
eb ssh
tail /opt/python/log/django.log
```

To login to prod postgres:

```
eb ssh
psql -h aa798nzxm9ji03.cpmmq0ugoybf.us-west-2.rds.amazonaws.com -U ebroot ebdb
```

Hint: Look in `/opt/python/current/env` for secrets.
