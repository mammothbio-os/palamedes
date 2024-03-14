install-deps:
	pip install -r requirements_dev.txt --upgrade

install:
	pip install -e .

test:
	py.test tests/ \
	--cov-branch \
	--cov-report term-missing \
	--cov palamedes

lint: 
	isort --profile black --check-only palamedes/ tests/
	flake8 palamedes/ tests/
	black --check --skip-string-normalization --line-length 120 palamedes/ tests/

clean:
	isort --profile black palamedes/ tests/
	black --skip-string-normalization --line-length 120 palamedes/ tests/
