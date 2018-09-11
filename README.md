# crispycrunch
Web app for CRISPR experiment setup and analysis. Fill in a 96-well plate and analyze it as a whole.

## Background

This app was developed by Greg Dingle starting in July 2018 for Andy May's genome engineering group in the CZ Biohub. The intention is to streamline the CRISPR experimentation process of researchers in the Biohub and provide a permanent record of experiments and results.

## Services

Crispycrunch calls out to several bioinformatics web services.

* [Crispor](http://crispor.tefor.net/)
* [Crispresso2](http://crispresso.pinellolab.partners.org)
* [Togows](http://togows.org)
* [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgTracks)

## Tech stack

* Postgres
* Python3
* Django
* Bootstrap

# Installation

Clone the repo. You must have czbiohub github access.

```git clone git@github.com:czbiohub/crispycrunch.git```

Build the docker images.

```cd crispycrunch && docker-compose build```

This will also run `pip install -r requirements.txt`.

Start the app services.

```docker-compose up```

Initialize the database.

```docker-compose exec web python manage.py migrate```

# Usage
<!-- TODO: better homepage -->

Go to http://localhost:8000/main/experiment. Create a new experiment and follow the subsequent steps.

## Admin

To inspect the state of the app, use the built-in Django admin interface at http://localhost:8000/admin/main/ . But first you must create an admin account with `python manage.py createsuperuser`.