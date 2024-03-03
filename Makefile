CODE := src/objectnat

build-and-publish: clean build publish

lint:
	poetry run pylint $(CODE)

format:
	poetry run isort $(CODE)
	poetry run black $(CODE)

install:
	pip install .

install-dev:
	poetry install --with dev

install-dev-pip:
	pip install -e . --config-settings editable_mode=strict

clean:
	rm -rf ./dist

build:
	poetry build

publish:
	poetry publish

install-from-build:
	python -m wheel install dist/graph_lib-*.whl
