.PHONY: all lint test test-cov test-docker install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	python -m pytest --cov=q2_ntcard -n 4 && coverage xml -o coverage.xml

test-docker: all
	qiime info
	qiime ntcard --help

install: all
	$(PYTHON) -m pip install -v .

dev: all
	pip install pre-commit
	pip install -e .
	pre-commit install

clean: distclean

distclean: ;
