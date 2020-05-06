# CrispyCrunch
Web app for CRISPR experiment setup and analysis. Fill in a 96-well plate and analyze it as a whole.

See it live at https://crispycrunch.czbiohub.org.

## Background

This app was developed by Greg Dingle starting in July 2018 for the genome engineering group in the CZ Biohub. The intention was to streamline the CRISPR experimentation process of researchers in the Biohub and provide a permanent record of experiments and results.

## Services

CrispyCrunch calls out to several bioinformatics web services.

* [Crispor](http://crispor.tefor.net/)
* [Crispresso2](http://crispresso.pinellolab.partners.org)
* [Togows](http://togows.org)
* [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks)
* [Ensemble](http://rest.ensembl.org/)
* [GGGenome](http://gggenome.dbcls.jp)

Currently, CrispyCrunch calls out to mirrors for Crispor and Crispresso:
* [Crispor mirror](http://ec2-34-222-186-65.us-west-2.compute.amazonaws.com/crispor.py)
* [Crispresso mirror](http://ec2-34-222-186-65.us-west-2.compute.amazonaws.com:81/)

*NOTE:* The Crispor mirror is modified to support HDR primer design. See [this PR](https://github.com/maximilianh/crisporWebsite/pull/21).

*NOTE:* The Crispresso mirror is running in a Docker container on the same machine as Crispor. Docker service is modified on the machine to store all its data in `/mnt/data/docker`. The code of `CRISPRessoCORE.py` was modified from the public image to fix a bug in HDR stats. The number of Celery workers was increased from the default of 1 to 3 for better parallelism.

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

```python manage.py migrate```

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

Wait for deployment to finish, then view the live site at http://crispycrunch.czbiohub.org/.

For more info, see https://realpython.com/deploying-a-django-app-and-postgresql-to-aws-elastic-beanstalk/ and https://hashedin.com/blog/5-gotchas-with-elastic-beanstalk-and-django/ .

To download fastq files stored in s3, you must configure
```
AWS_ACCESS_KEY_ID
AWS_SECRET_ACCESS_KEY
```
See https://docs.aws.amazon.com/elasticbeanstalk/latest/dg/environments-cfg-softwaresettings.html and https://boto3.amazonaws.com/v1/documentation/api/latest/guide/configuration.html#environment-variable-configuration.

# Admin

To inspect the state of the app, use the built-in Django admin interface at http://localhost:8000/admin/main/ or http://crispycrunch.czbiohub.org/admin/.

A default `admin` user is created on first installation. **Change the password** of that first superuser on first deployment.

You can give any other user admin or superuser status here: https://crispycrunch.czbiohub.org/admin/auth/user/ .

To inspect the source of fastq files, open https://console.aws.amazon.com/s3/buckets/czb-seqbot/.

To read the error logs in prod:

```
eb ssh
tail /opt/python/log/django.log
```

To login to prod postgres:

```
eb ssh
psql -h RDS_HOSTNAME -U RDS_USERNAME RDS_DB_NAME
```

Hint: Look in `/opt/python/current/env` for secrets.

To create pre-authenticaed URLs for sharing using https://github.com/aaugustin/django-sesame:

```
# local
eb ssh
# remote
cd /opt/python/current/
source env
source /opt/python/run/venv/bin/activate
cd app && python manage.py shell
# python shell
>>> from django.contrib.auth.models import User
>>> user = User.objects.get(username='demo')
>>> from sesame import utils
>>> utils.get_query_string(user)
```

# Monitoring

There should be at least one downtime alert setup here:
https://console.aws.amazon.com/cloudwatch/home?region=us-west-2#s=Alarms&alarm=awseb-e-vnvbedub4n-stack-Severe-KUZS17WH3VE9

And you should receive emails of Django errors if you are the owner of `ADMIN_EMAIL`. They are throttled with https://github.com/krisys/django-error-email-throttle .

You can see a table of error counts at https://crispycrunch.czbiohub.org/admin/error_email_throttle/

