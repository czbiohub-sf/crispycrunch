set -o errexit
find . -path "*/migrations/*.py" -not -name "__init__.py" -delete
psql postgres -c 'drop schema if exists public cascade; create schema public;'
python manage.py makemigrations
python manage.py migrate

# For manual of prod, do:
# eb ssh
# psql -h aa798nzxm9ji03.cpmmq0ugoybf.us-west-2.rds.amazonaws.com -U ebroot ebdb
