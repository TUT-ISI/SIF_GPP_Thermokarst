#!/bin/bash -l

for year in {2005..2020}
do
	python convertAllDataToNC.py ${year}
done
