DEV_DB_CONTAINER_NAME = crispycrunch-dev
TEST_DB_CONTAINER_NAME = crispycrunch-test

.PHONY: init
init:
	pip install -r requirements.txt

.PHONY: setup-develop
setup-develop:
	pip install -e .'[dev]'
	pre-commit install

.PHONY: pre-commit
pre-commit:
	pre-commit run --all-files

.PHONY: lint
lint:
	flake8 . --count --statistics --exit-zero
	python -m pylint ./opencell

.PHONY: test
test:
	pytest -v --ignore ./client/node_modules

.PHONY: start-dev-db
start-dev-db: drop-dev-db
	docker create \
	--name $(DEV_DB_CONTAINER_NAME) \
	-e POSTGRES_USER=crispycrunch \
	-e POSTGRES_PASSWORD=password \
	-e POSTGRES_DB=crispycrunchdb \
	-p 5435:5432 \
	postgres
	docker start $(DEV_DB_CONTAINER_NAME) && sleep 2

.PHONY: drop-dev-db
drop-dev-db:
	-docker rm --force $(DEV_DB_CONTAINER_NAME);

start-app:
	python manage.py runserver