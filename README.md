# SIF_GPP_Thermokarst
This repository contains the most relevant functions generated to assess and plot the Arctic SIF and GPP sensitivity analysis to air temperature in different thermokarst-affected areas.

The functions contained in this repository were used to compute, plot and perform the analysis presented in the ArcticSIF project article:
contact author: neus.sabater@fmi.fi

### Instructions
1. Open Matlab or python
2. Run the script `script1.m` to start the analysis.
3. Ensure that the input files are placed in the appropriate folder. The input and output paths can be defined in the environment section in the conf.json file.


### Environment functions 
Most of the scripts above work by detecting in which environment the code is running (local machine or server) and then reading the corresponding paths where the input data is stored.

- `config.json`: Specifies the list of paths for input/output data for each environment. Edit this to your own paths if you want to run the codes.
- `detect_environment.m`: This function detects your local environment or server. Edit with your local environment name.
- `detect_environment.py`: Script to detect your local environment or server. Edit with your local environment name.
- `example_read_env.m`: Example script demonstrating how to read the environment in a Matlab function.
- `example_read_env.py`: Example script demonstrating how to read the environment in a Python function.
- `load_conf.m`: Function to load the configuration file.
- `load_conf.py`: Script to load the configuration file.



## Matlab

Here we list the Matlab scripts and functions used in the project. 
Make sure to have Matlab installed with the necessary toolboxes: Statistical and machine learning toolbox:

- `computeThermokarstPeatlandCStorage.m:`: This script reads the Thermokarst map and the permafrost-affected Peatland C storage dataset to check the degree of consistency between these two databases.
- `computeNDVI_EVI_MODIS.m`: This script reads the EVI and NDVI products from MODIS Terra and Aqua satellites and re-grids them to the Arctic grid used in this project.
- `computeNDVI_EVI_MODISAQUATERRA_Merging.m`: This script reads the EVI and NDVI NetCDF products from MODIS Terra and Aqua satellites and regrids them to the ArcticSIF 0.1-degree resolution with an 8-day interval.
- `computeSIF_OCO2.m`: This script reads the OCO-2 data from NetCDF files for the years 2015, 2016, and 2017 to validate whether the observed increase in SIF during 2016 (seen with GOSIF) is corroborated and not an artifact due to training with meteorological datasets.
- `computeBAWLDPercents.m`: This script reads the BAWLD database, along with percent information, to compute a stacked bar chart and estimate the percentage of different elements within grid cells classified as each class.
- `computeSIFAno_CSIF.m`: This script reads the cSIF database generated with OCO-2 SIF data and MODIS surface reflectance every 4 days at 0.05-degree resolution.
- `readBAWLD.m`: This script reads the `.shp` and `.prj` files for the land and lake Arctic database, published at [https://arcticdata.io/catalog/view/doi:10.18739/A2C824F9X](https://arcticdata.io/catalog/view/doi:10.18739/A2C824F9X).


Plotting functions

- `plotAnoAirTempThermkarst.m`: Plots air temperature anomalies saved and estimated during the period of interest for the study and evaluates its range per class of BAWLD.
- `plotConfMatrixAnomERA5_ARCLIM.m`: Generates a confusion matrix between the anomalies of the different ERA5 and ARCLIM variables and the anomalies of GPP and SIF to detect if the response is in sync.
- `plotPolarPeatland.m`: Reads peatland information (C-store stocks in peatlands) data and plots the corresponding maps.
- `plotPolarPlotAnomAirTempARCLIM.m`: Generates polar plots for the ARCLIM 8-day interval data used in the SHAP analyses across different years.
- `plotPolarPlotAnomCSIFAll_Thermokast.m`: Reads the cSIF database anomalies from 2003-2016 (or 2003-2018) and computes the OLS regression with the corresponding air temperature anomalies for the same period.
- `plotPolarPlotAnomERA5.m`: Generates polar plots for the ECMWF 8-day interval data used in the SHAP analyses across different years.
- `plotPolarNDVI_EVI.m`: Reads the EVI and NDVI average values in June/July for each year generated using the MODIS Terra and Aqua data, then generates polar plots.
- `plotPolarPlotClassSIFGPPTrend.m`: Generates polar plots for the SIF & GPP classification trends. There are four classes combining negative and positive GPP and SIF trends.
- `plotAnomAirTempBAWLD.m`: Plots air temperature anomalies saved and estimated during the period of interest for the study and evaluates its range per class of BAWLD.
- `plotConfusionMatrixThermoskarstBAWLD.m`: Reads the Thermokarst classification and the BAWLD class, generating a matrix that relates the percent of Thermokarst coverage for each class.
- `plotPolarPlotAnomCSIFAll_Thermokast.m`: Reads the cSIF database anomalies from 2003-2016 (or 2003-2018) and computes the OLS regression with the corresponding air temperature anomalies for the same period.
- `plotPolarPlotAnomCSIFAll_class.m`: Reads the cSIF database anomalies from 2003-2016 (or 2003-2018) and computes the OLS regression with the corresponding air temperature anomalies for the same period.
- `plotPolarThermokast.m`: Reads the maps of thermokarst of permafrost available at [https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332](https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1332). The original maps (in `.shp` and `.prj` formats) have been re-gridded to the lat-lon grid used in the ArcticSIF study. The function to generate the regridded maps and store them in NetCDF is `readThermokast.py`.
- `plotPolarPlotAnomAirTemp.m`: Generates polar plots for air temperature anomalies data from ECMWF 8-day intervals over different years.


## Functions_Matlab Generic
This section provides information about the generic functions used across various Matlab scripts.
- `redblue.m` define a redblue scale


## Shap Analysis

The Shap analysis is critical for explaining machine learning models' predictions. This section focuses on SHapley Additive exPlanations (SHAP) functions.
This section is all coded in Python.

- `Calculate_2t.py`: This script reads from the ARCLIM database the 2m temperature, converts it into Celsius, and crops the latitude region of interest.
- `Calculate_vpd.py`: This script reads from the ARCLIM database the dewpoint temperature and 2m temperature to compute the vapor pressure deficit (VPD).
- `convertAllDataToNC.py`: This script converts the SIF and GPP data into NetCDF format at the required resolution for the project.
- `ERA5Downloader.py`: This script downloads the required meteorological variables from the ECMWF ERA-5 products. The data is aggregated per year and month and stored in NetCDF format.
- `ERA5GPP_preprocesser.py`: This script prepares and merges all the input data required to perform the SHAP analysis for GPP.
- `ERA5SIF_preprocesser.py`: This script prepares and merges all the input data required to perform the SHAP analysis for SIF.
- `GPPEvaluator.py`: This script evaluates the model trained for GPP to be used in the SHAP analysis.
- `SIFEvaluator.py`: This script evaluates the model trained for SIF to be used in the SHAP analysis.
- `Model_trainer_GPP.py`: This script trains the model for the GPP analysis and stores the results in a Parquet output format.
- `Model_trainer_SIF.py`: This script trains the model for the SIF analysis and stores the results in a Parquet output format.
- `SHAP_analysis_ERA_mass.py`: This script generates the beeswarm plots for SHAP analysis.
- `SHAP_beeswarm.py`: This script plots the beeswarm diagrams for SHAP analysis.



### Python 
General scripts coded in Python 
- `convertPeatlandTifToNC.py`: Converts information from the Histel and Histosol datasets from `.tiff` to `.nc`.
- `readThermokast.py`: Extracts and projects data from the shapefile of the Thermokarst database.
- `convertDataToNC.py`: Reads and stores the GPP data from the Fluxsat_v2 database and SIF from the GOSIF databases in a NetCDF output format.
- `convertSifSunShadeDataToNc.py`: Reads and stores the 8-day combined GPPsun, GPPShade, and GPPtotal data in a NetCDF output format. Dataset source: [https://doi.org/10.1038/s41597-022-01309-2](https://doi.org/10.1038/s41597-022-01309-2).
