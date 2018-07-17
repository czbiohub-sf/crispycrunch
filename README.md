# crispycrunch
Web app for CRISPR experiment setup and analysis

## Background

This app was developed by Greg Dingle starting in July 2018 for Andy May's genome engineering group in the CZ Biohub. The intention is to streamline the CRISPR experimentation process of researchers in the Biohub and provide a permanent record of experiments and results. 

## Installation

Clone the repo. You must have czbiohub github access.

```git clone git@github.com:czbiohub/crispycrunch.git```

## Technologies

* Python3 
* Django
* Postgres
* Bootstrap

Build and run the app. Crispycrunch is packaged in a docker container. 

```cd crispycrunch && docker-compose up --build```

Check that it builds and runs successfully. The app is currently served by the django development server. You should see something like the following.

```
Building web
Step 1/7 : FROM python:3
 ---> 17453243214e
Step 2/7 : ENV PYTHONUNBUFFERED 1
 ---> Using cache
 ---> 8490e4f9d7c8
Step 3/7 : RUN mkdir /code
 ---> Using cache
 ---> be01f152335d
Step 4/7 : WORKDIR /code
 ---> Using cache
 ---> 53b563c47442
Step 5/7 : ADD requirements.txt /code/
 ---> Using cache
 ---> e304535703df
Step 6/7 : RUN pip install -r requirements.txt
 ---> Using cache
 ---> 8dc0f2c6e308
Step 7/7 : ADD . /code/
 ---> f0415afb3591
Successfully built f0415afb3591
Successfully tagged crispycrunch_web:latest
Starting crispycrunch_db_1 ... done
Recreating crispycrunch_web_1 ... done
Attaching to crispycrunch_db_1, crispycrunch_web_1
db_1   | 2018-07-17 16:54:13.352 UTC [1] LOG:  listening on IPv4 address "0.0.0.0", port 5432
db_1   | 2018-07-17 16:54:13.352 UTC [1] LOG:  listening on IPv6 address "::", port 5432
db_1   | 2018-07-17 16:54:13.359 UTC [1] LOG:  listening on Unix socket "/var/run/postgresql/.s.PGSQL.5432"
db_1   | 2018-07-17 16:54:13.384 UTC [23] LOG:  database system was shut down at 2018-07-17 16:08:13 UTC
db_1   | 2018-07-17 16:54:13.409 UTC [1] LOG:  database system is ready to accept connections
web_1  | Performing system checks...
web_1  | 
web_1  | System check identified no issues (0 silenced).
web_1  | July 17, 2018 - 09:54:15
web_1  | Django version 2.0.7, using settings 'crispycrunch.settings'
web_1  | Starting development server at http://0.0.0.0:8000/
web_1  | Quit the server with CONTROL-C.
```
