set -o errexit
find . -path "*/migrations/*.py" -not -name "__init__.py" -delete
psql postgres -c 'drop schema if exists public cascade; create schema public;'
python manage.py makemigrations
python manage.py migrate
