"""
Script: This script plots the beeswarm diagrams. 
Author: Lauri Hatakka
Version: 1.0
Date: Jun, 2023

"""
import numpy as np
import pandas as pd
import dask.dataframe as dd
import pickle
import seaborn as sns
import random
import ipywidgets as widgets
import shap 
import matplotlib.pyplot as plt
import xarray as xr
import sys
import time
import itertools

variables = ['GPP', 'SIF']
dropped_variable_groups = [[], ['lai_hv', 'lai_lv', 'fal']]
month_groups = [[4,5],[6,7],[8,9]]

sample_size=10000000

listoflists = [variables, dropped_variable_groups, month_groups]
listofoptions = list(itertools.product(*listoflists))


index = int(sys.argv[1])
options = listofoptions[index]

variable = options[0]
dropped_variables = options[1]
months_to_analyse = options[2]

print(f'{variable}/{dropped_variables}/{months_to_analyse}', flush=True)

starttime = time.time()

if variable == 'GPP':
    location_modifier = f'ERAGPP'
if variable == 'SIF':
    location_modifier = f'ERAFLUXSAT'

location = f'/fmi/scratch/project_2005798/data/Modelfiles/{location_modifier}/{months_to_analyse}/Dropped_{dropped_variables}'
figsavelocation = f'/fmi/scratch/project_2005798/figures/beeswarm'

explainer = pickle.load( open(f'{location}/explainer_{sample_size}.pickle', "rb") )
plot = shap.plots.beeswarm(explainer, show=False)

#To modify beeswarm plots, use normal pyplot syntax here!
plt.suptitle(f'SHAP values for {variable}')
plt.title(f'Months: {months_to_analyse}, Dropped: {dropped_variables}')

plt.savefig(f'{figsavelocation}/beeswarm_{variable}_{dropped_variables}_{months_to_analyse}.png')
plt.close()