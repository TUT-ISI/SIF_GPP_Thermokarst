"""
Script: This script trains the model for the GPP analysis and store the results in a parquet output format. 
Author: Lauri Hatakka
Version: 1.0
Date: Jun, 2023

"""
import numpy as np
import pandas as pd
import dask.dataframe as dd
import dask.array as da
import pickle
import xgboost as xgb
import xarray as xr
# from sklearn.model_selection import train_test_split
from dask_ml.model_selection import train_test_split
from dask_ml.wrappers import ParallelPostFit
import random
import os

variable_groups = [[], ['lai_hv', 'lai_lv', 'fal']]
month_groups = [[4,5],[6,7],[8,9]]

for variable_to_drop in variable_groups:
    for months_to_analyse in month_groups:
        datapath = '/fmi/scratch/project_2005798/data/'
        years = list(range(2000,2021))
        dataset = xr.Dataset()
        for year in years:
            try:
                Yeardata = xr.open_dataset(os.path.join(datapath, f'PreprocessedData/ERAGPP/{months_to_analyse}/{year}.nc'), chunks = 'auto' )
                dataset = xr.merge([dataset, Yeardata])
            except FileNotFoundError:
                print(f'{year} not found')
        print(dataset)
        dataset = dataset.drop_indexes(['latitude', 'longitude', 'time'])
        dataset = dataset.chunk('auto').unify_chunks()
        dataset = dataset.to_dask_dataframe()
        dataset = dataset.drop('ssr', axis =1)
        dataset = dataset.dropna(how = 'all',)
        dataset = dataset.dropna(how = 'any', subset = ['GPPtotal'])
        dataset = dataset.drop(variable_to_drop, axis = 1)



        outputdir = f'/fmi/scratch/project_2005798/data/Modelfiles/ERAGPP/{months_to_analyse}/Dropped_{variable_to_drop}'
        os.makedirs(outputdir, exist_ok = True)

        dataset.to_parquet(f'{outputdir}/dataset.parquet')
        print('conversion to dataframe done!', flush = True)
        orig_indexes=dataset[['time', 'latitude', 'longitude']]
        dataset = dataset.drop(['latitude','longitude', 'time'], axis = 1)
        
        X = dataset.drop(['GPPtotal', 'GPPShade', 'GPPSun'], axis = 1)
        y = dataset['GPPtotal']

        X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42, test_size=0.2, shuffle=True)

        print(type(X_train))
        print(type(X_test))
        print(type(y_train))
        print(type(y_test))


        model = xgb.XGBRegressor(objective="reg:squarederror", 
                                n_estimators=100,
                                max_depth=10,
                                learning_rate=0.1,
                                verbose=1,
                                missing=float('nan'))



        print('Fitting!', flush = True)

        model.fit(X_train, y_train)
        filename = f'{outputdir}/XGBoost.pickle'
        pickle.dump(model, open(filename, "wb"))

        print(X_test)

        X_test.to_parquet(f'{outputdir}/x_test.parquet')
        X_train.to_parquet(f'{outputdir}/x_train.parquet')
        y_test.to_frame().to_parquet(f'{outputdir}/y_test.parquet')
        y_train.to_frame().to_parquet(f'{outputdir}/y_train.parquet')

        orig_indexes.to_parquet(f'{outputdir}/latlontime.parquet')


