"""
Script: This script generates the beeswarm plots. 
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


# variables = ['GPP', 'SIF']
# dropped_variable_groups = [[], ['lai_hv', 'lai_lv', 'fal']]
# month_groups = [[4,5],[6,7],[8,9]]
variables = ['GPP']
dropped_variable_groups = [[]]
month_groups = [[6,7],[8,9]]
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
    
model = pickle.load(open(f'{location}/XGBoost.pickle', 'rb'))
x_test = pd.read_parquet(f'{location}/x_test.parquet')
x_train = pd.read_parquet(f'{location}/x_train.parquet') 
y_test = pd.read_parquet(f'{location}/y_test.parquet') 
y_train = pd.read_parquet(f'{location}/y_train.parquet') 
locales = pd.read_parquet(f'{location}/latlontime.parquet')

x = pd.concat([x_train, x_test])
y = pd.concat([y_train, y_test])


x_testsample = x.sample(sample_size, random_state=300)


explainer = shap.TreeExplainer(model=model)
shap_values = explainer.shap_values(x_testsample)

explainer_for_sample = explainer(x_testsample)


explanation_values = xr.DataArray(explainer_for_sample.values,dims=['index', 'column'], coords = [x_testsample.index, x_testsample.columns] ).to_dataset(name = 'shap_values')
explanation_base_values = xr.DataArray(explainer_for_sample.base_values, dims = ['index'], coords = [x_testsample.index]).to_dataset(name = 'base_values')
explanation_data = xr.DataArray(explainer_for_sample.data, dims =['index','column'], coords = [x_testsample.index, x_testsample.columns]).to_dataset(name = 'original_data')
explanation_locales = locales.to_xarray().rename({'__null_dask_index__':'index'})
shap_ds = xr.merge([explanation_locales, explanation_values,explanation_base_values, explanation_data])
shap_ds.to_netcdf(f'{location}/SHAP_{sample_size}_point_sample_{variable}_{dropped_variables}_{months_to_analyse}')
print(f'Dumping! {location}/explainer_{sample_size}.pickle')
pickle.dump(explainer_for_sample, open(f'{location}/explainer_{sample_size}.pickle', "wb") )
pickle.dump(x_testsample, open(f'{location}/testsample_{sample_size}.pickle', "wb"))
# np.set_printoptions(threshold=sys.maxsize)


print(f'For dropped: {dropped_variables} and months: {months_to_analyse} and {variable}, the beeswarm plot for {variable} is:')
# summary_plot = shap.summary_plot(shap_values, x_testsample, show=False)
plot = shap.plots.beeswarm(explainer(x_testsample), show=False)
plt.suptitle(f'SHAP values for {variable}')
plt.title(f'Months: {months_to_analyse}, Dropped: {dropped_variables}')

plt.savefig(f'{figsavelocation}/beeswarm_{variable}_{dropped_variables}_{months_to_analyse}.png')
plt.close()

print(f'For sample size of {sample_size}: Running time: {time.time() - starttime}')