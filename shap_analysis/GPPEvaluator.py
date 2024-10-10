"""
Script: This script evaluates the model trained for GPP to be used in the SHAP analysis. 
Author: Lauri Hatakka
Version: 1.0
Date: Jun, 2023

"""
import numpy as np
import pandas as pd
import dask.dataframe as dd
import dask.array as da
import pickle
from sklearn.metrics import mean_squared_error
from dask_ml.model_selection import train_test_split
import matplotlib.pyplot as plt
import sys
import warnings
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LogNorm



def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    bit_rgb = np.linspace(0, 1, 256)
    if position is None:
        position = np.linspace(0, 1, len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red': [], 'green': [], 'blue': []}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap', cdict, 256)
    return cmap


cmap2Dhist = make_cmap([(255, 255, 255), (127, 188, 227), (82, 170, 115), (225, 189, 74), (222, 64, 39), (149, 21, 25)], bit=True)

def computeStats(true, predicted, computeEEratio=False):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='invalid value encountered in true_divide')  # temporarily ignore this warning (would give warnings if predicted only contains a single value)
        N = len(true)
        R = np.corrcoef(true, predicted)[0, 1]
        R2 = R**2  # r2_score(true, predicted)
        BIAS = np.median(predicted) - np.median(true)
        
    return R, R2, BIAS

def calculate_errorbars_in_range(x,y,start_val,stop_val,step):
    xx = np.arange(start_val,stop_val,step)
    
    x_loc = np.zeros(xx.size-1)
    y_mean = np.zeros(xx.size-1)
    y_std = np.zeros(xx.size-1)
    
    for i in range(np.size(xx)-1):
        ind = np.where(np.logical_and(x >= xx[i], x < xx[i+1]))
        
        vals = y[ind]
        
        if vals.size == 0:
            y_mean[i] = np.nan
            y_std[i] = np.nan
        else:
            y_mean[i] = np.nanmean(vals)
            y_std[i] = np.nanstd(vals)
            
        x_loc[i] = np.nanmean((xx[i],xx[i+1]))

    return x_loc, y_mean, y_std

def plot_data(xx, yy, y_true, title_name, RMSE, x_str, y_str, save_name, error_bar_on=True, plt_linear=1, path="/home/vaisanea/Desktop/"):
    print("Starting plotting: " + title_name)
    
    yy = np.array(yy)
    xx = np.array(xx)
    
    x_min = np.min(xx) - 0.5
    x_max = np.max(xx) + 0.5
    y_min = np.min(yy) - 0.5
    y_max = np.max(yy) + 0.5
    
    N = int(np.size(xx))
    
    fig = plt.figure()
    ax = fig.subplots()
  
        
    norm2Dhist = LogNorm(vmin=1, vmax=int(0.025 * N))
    
    ax.hist2d(xx, yy, bins=(64,64), range=[[x_min,x_max],[y_min,y_max]], cmap=cmap2Dhist, norm=norm2Dhist, alpha=0.8, zorder=2)
    ax.grid(True)
    
    ax.set_ylabel(y_str)
    ax.set_xlabel(x_str)
    ax.set_title(title_name)
    ax.set_ylim([y_min,y_max])
    ax.set_xlim([x_min,x_max])
    
    if plt_linear == 1:
        ax.plot([-10, 20], [-10, 20], 'k-', linewidth=2.0, alpha=0.8, zorder=3)
    elif plt_linear == 2:
        ax.plot([x_min, y_max], [x_min, y_max], 'k-', linewidth=2.0, alpha=0.8, zorder=3)
    else:
        ax.plot([-10, 20], [0, 0], 'k-', linewidth=2.0, alpha=0.8, zorder=3)
        
    axCB = fig.add_axes([0.935, 0.1, 0.02, 0.85])
    cb1 = mpl.colorbar.ColorbarBase(axCB, cmap=cmap2Dhist, norm=norm2Dhist, orientation='vertical')
    cb1.set_label('N', fontsize=24)
    
    if error_bar_on:
        x_error_loc, y_mean, y_std = calculate_errorbars_in_range(xx,yy,np.min(xx),np.max(xx),1)
        ax.errorbar(x_error_loc, y_mean, y_std,marker='o',ecolor="m",markerfacecolor='m',linestyle='None',capsize=5,zorder=4)
    
    R, R2, BIAS = computeStats(y_true,yy)
    
    if plt_linear:
        ax.text(0.99, 0.08 + 3 * 0.08, 'N={}'.format(N), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, fontsize=18, color='k', zorder=99)#, bbox={'facecolor': 'white', 'alpha': 0.7, 'boxstyle': 'square,pad=0.03'})
        ax.text(0.99, 0.08 + 2 * 0.08, 'R2={:.3f}'.format(R2), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, fontsize=18, color='k', zorder=99)#, bbox={'facecolor': 'white', 'alpha': 0.7, 'boxstyle': 'square,pad=0.03'})
        ax.text(0.99, 0.08 + 1 * 0.08, 'BIAS={:.3f}'.format(BIAS), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, fontsize=18, color='k', zorder=99)#, bbox={'facecolor': 'white', 'alpha': 0.7, 'boxstyle': 'square,pad=0.03'})
        ax.text(0.99, 0.08 + 0 * 0.08, 'RMSE={:.3f}'.format(RMSE), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, fontsize=18, color='k', zorder=99)#, bbox={'facecolor': 'white', 'alpha': 0.7, 'boxstyle': 'square,pad=0.03'})
    else:
        ax.text(0.99, 0.08 + 0 * 0.08, 'N={}'.format(N), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, fontsize=18, color='k', zorder=99)#, bbox={'facecolor': 'white', 'alpha': 0.7, 'boxstyle': 'square,pad=0.03'})

    fig.savefig(path + save_name + ".png", bbox_inches="tight")

    print("Done!\n")

variable_groups = [[],['lai_hv', 'lai_lv', 'fal']]
month_groups = [[4,5],[6,7],[8,9]]

figoutput = f'/fmi/scratch/project_2005798/figures/eval'

# dropped_variables = ['lai_hv', 'lai_lv', 'fal']
# dropped_variables = []
# months_to_analyse = [8,9]
for dropped_variables in variable_groups:
    for months_to_analyse in month_groups:
        inputpath = f'/fmi/scratch/project_2005798/data/Modelfiles/ERAGPP/{months_to_analyse}/Dropped_{dropped_variables}'
        x_test = dd.read_parquet(f'{inputpath}/x_test.parquet')
        y_test = dd.read_parquet(f'{inputpath}/y_test.parquet').compute()

        model = pickle.load(open(f'{inputpath}/XGBoost.pickle', 'rb'))

        y_pred = model.predict(x_test)

        mse = mean_squared_error(y_test, y_pred)
        model_name = 'xgboost'
        target = 'GPP'
        print(y_test)
        y_test = y_test['GPPtotal']
        plot_data( y_test, y_pred,
                #y_true = y_test[target],
                y_true = y_test,
                y_str = 'Prediction',
                x_str = 'Ground Truth',
                error_bar_on=False, 
                save_name=f'eval_{model_name}_{target}_{dropped_variables}_{months_to_analyse}', 
                title_name=f'{model_name} {target}',
                RMSE= mse**(1/2.0),
                path = f'{figoutput}/')

