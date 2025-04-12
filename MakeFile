ENVIRONMENT = chem_analyzer_env

environment:
	conda env create -f environment.yaml

tests:
	pytest -v test_molecule_final.py
