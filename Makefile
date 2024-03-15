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
	ruff check palamedes/ tests/
	mypy palamedes/ tests/

clean:
	ruff format palamedes/ tests/
