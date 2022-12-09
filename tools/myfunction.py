#!/usr/bin/env python
# coding: utf-8

# %matplotlib inline

from netCDF4 import Dataset
import numpy as np
import pandas as pd
#import mpld3
import xarray as xr
import datetime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#import matplotlib.colors as mcolors
import seaborn as sns
from datetime import datetime
import math
import scipy
import scipy.cluster.hierarchy as sch
from math import fsum
#from varname import nameof
from datetime import date
from scipy import stats
import os
from scipy.stats import pearsonr, spearmanr
# sns.set_theme(style="ticks")

# Define some lists/arrays
depth = np.array([0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40,
                 0.44, 0.54, 0.64, 0.74, 0.84, 0.94, 1.04, 1.14, 2.39, 4.676, 7.635, 11.140, 15.115])
node_z = [0.01, 0.04, 0.09, 0.16, 0.26, 0.40, 0.587, 0.80, 1.06, 1.36, 1.70, 2.08,
          2.50, 2.99, 3.58, 4.27, 5.06, 5.95, 6.94, 8.03, 9.795, 13.328, 19.483, 28.871, 41.998]
sec_pr_day = 60*60*24
# layer depths to multiply with concentrations
depth_30 = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
depth_50 = np.array([0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.02])  # [meters]

pool_names = ['LITm', 'LITs', 'SAPb', 'SAPf', 'EcM',
              'AM', 'SOMp', 'SOMa', 'SOMc']  # 'ErM', 'AM',
N_pool_names = ['N_LITm', 'N_LITs', 'N_SAPb', 'N_SAPf', 'N_EcM',
                'N_AM', 'N_SOMp', 'N_SOMa', 'N_SOMc', 'NH4_sol', 'NH4_sorp', 'NO3']

PFT_names = ["not_vegetated                           ",
             "needleleaf_evergreen_temperate_tree     ",
             "needleleaf_evergreen_boreal_tree        ",
             "needleleaf_deciduous_boreal_tree        ",
             "broadleaf_evergreen_tropical_tree       ",
             "broadleaf_evergreen_temperate_tree      ",
             "broadleaf_deciduous_tropical_tree       ",
             "broadleaf_deciduous_temperate_tree      ",
             "broadleaf_deciduous_boreal_tree         ",
             "broadleaf_evergreen_shrub               ",
             "broadleaf_deciduous_temperate_shrub     ",
             "broadleaf_deciduous_boreal_shrub        ",
             "c3_arctic_grass                         ",
             "c3_non-arctic_grass                     ",
             "c4_grass                                ", ]

PFT_names = [i.strip() for i in PFT_names]

# DEFINING FUNCTIONS:


def get_observations(list_of_site_numbers):

    df = pd.read_excel(
        "/home/ecaas/Desktop/jordprofiler/CNstocks2016_HAdW_updated_LTS9august21.xlsx", sheet_name="Sheet1")
    komm = pd.read_excel(
        "/home/ecaas/Desktop/jordprofiler/Kommuneklassifisering_1994.xlsx", sheet_name="Sheet1")
    alle = df[1:]  # All lines containing values (all but first line)
    info = df.iloc[0]  # First line contains more info

    # Change dtype to approproate types
    for c in alle.columns:
        if alle.loc[:, c].dtype == object:
            if type(alle.loc[:, c].iloc[0]) != str:
                alle.loc[:, c] = pd.to_numeric(alle.loc[:, c])
            else:
                alle.loc[:, c] = alle.loc[:, c].astype('category')

    # Add kommune as a column to the data frame "alle"
    kommune = []
    for knr in alle.KOMM:
        if knr not in komm.loc[:, 'Kommune-nr.'].values:
            kommune.append('NO MATCH')
        else:
            for snr in komm.loc[:, 'Kommune-nr.']:
                if knr == snr:
                    kommune.append(
                        komm[komm['Kommune-nr.'] == snr]['Kommunenavn'].values[0])
    alle.loc[:, 'Kommunenavn'] = kommune

    # Chosen locations for decomp development (Kongsvinger and Saltdal)
    sites = list_of_site_numbers
    sites_df = alle.loc[alle['PROFNR'].isin(sites)].reset_index()

    for i in range(len(sites_df)):
        site = ('NR'+str(sites_df.loc[i, 'PROFNR']) +
                '_'+sites_df.loc[i, 'Kommunenavn']).strip()
        sites_df.loc[i, 'SITE'] = site
        obs = sites_df.loc[:, ['SITE', 'CLFH', 'C30', 'C50', 'cs100', 'cstock', 'c100min', 'COMratio', 'LFHDYP', 'PROFDYPC', 'NLFH', 'nstock', 'ns100', 'n100min', 'NOMratio', 'CNFOREST', 'KOORDINA',
                               'MAP', 'MAT', 'TREG', 'VEGETASJ', 'ALTITUDE', 'PBONITET', 'VEKSTSES', 'DRAINAGE', 'FAOTex', 'CanTex', "CANADIAN", "VEGTYP", "FID", "BEST_TRE"]].sort_values(by='SITE')

    obs.loc[:, "origin"] = "obs"

    for col in obs[['LFHDYP', 'PROFDYPC']].columns:  # convert observed depths to meters
        obs.loc[:, col] = obs.loc[:, col]/100.0

    for col in obs[['cstock', 'nstock', 'C50', 'C30', 'MAP', 'MAT', "cs100", "ns100", "c100min", "n100min", "CLFH", "NLFH"]].columns:  # Round to two decimals
        obs.loc[:, col] = round(obs.loc[:, col], 2)

    return obs, info

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def make_decomp_df(xarray_list):
    decomp_Cpools_df = pd.DataFrame(columns=['SITE'] + pool_names)
    decomp_Npools_df = pd.DataFrame(columns=['SITE'] + N_pool_names)
    flux_list = ['C_PlantLITm', 'C_PlantLITs', 'C_PlantEcM',
                 'C_PlantAM', 'C_PlantSOMc', 'C_PlantSOMp']

    for i in range(len(xarray_list)):
        layers = no_of_layers(xarray_list[i])
        # add site name to val (which is a dataframe row)
        decomp_Cpools_df.loc[i, 'SITE'] = str(xarray_list[i].site_name)
        # add site name to val (which is a dataframe row)
        decomp_Npools_df.loc[i, 'SITE'] = str(xarray_list[i].site_name)
        # Put the time averaged C content of each soil pool column (total) in dataframe decomp_Cpools_df
        for pool in pool_names:
            mean = float(
                fsum(xarray_list[i][pool].mean(dim="time")*depth[0:layers]))
            decomp_Cpools_df.loc[i, pool] = mean

        # Put the time averaged N content of each soil pool column (total) in dataframe decomp_Npools_df
        for pool in N_pool_names:
            mean = float(
                sum(xarray_list[i][pool].mean(dim="time")*depth[0:layers]))
            decomp_Npools_df.loc[i, pool] = mean

        decomp_Cpools_df.loc[i, "TOTC"] = decomp_Cpools_df[pool_names].sum(axis=1)[
            i]
        # NOTE: excluding inorganic N (NO3 and NH4sol & NH4sorp)
        decomp_Npools_df.loc[i,
                             "TOTN"] = decomp_Npools_df[N_pool_names[:-3]].sum(axis=1)[i]

        # Put mean of total C input into dataframe decomp_Cpools_df
        for f in flux_list:
            tot = float(fsum(xarray_list[i][f].mean(
                dim="time")*depth[0:layers]))
            decomp_Cpools_df.loc[i, f] = tot

        decomp_Cpools_df.loc[i, "TOTC_input"] = decomp_Cpools_df[flux_list].sum(axis=1)[
            i]

        for f in ['HR_flux']:
            tot = float(fsum(xarray_list[i][f].mean(
                dim="time")*depth[0:layers]))
            decomp_Cpools_df.loc[i, f] = tot

        col = xarray_list[i][pool_names]
        col = col.mean(dim="time")
        col = col.to_dataframe().reset_index(drop=True)

        for var in col.columns:
            col30 = col[var][0:len(depth_30)]*depth_30
            sum30 = col30.sum()
            decomp_Cpools_df.loc[i, 'c30_%s' % var] = round(sum30, 2)

            col50 = col[var][0:len(depth_50)]*depth_50
            sum50 = col50.sum()
            decomp_Cpools_df.loc[i, 'c50_%s' % var] = round(sum50, 2)

        sum_SOM30 = decomp_Cpools_df.loc[i, 'c30_SOMa'] + \
            decomp_Cpools_df.loc[i, 'c30_SOMp'] + \
            decomp_Cpools_df.loc[i, 'c30_SOMc']
        decomp_Cpools_df.loc[i, 'c30_SOM_tot'] = sum_SOM30

        sum30 = decomp_Cpools_df.loc[i, 'c30_SOMa']+decomp_Cpools_df.loc[i, 'c30_SOMp']+decomp_Cpools_df.loc[i, 'c30_SOMc'] + \
            decomp_Cpools_df.loc[i, 'c30_LITm']+decomp_Cpools_df.loc[i, 'c30_LITs']+decomp_Cpools_df.loc[i, 'c30_SAPb'] + decomp_Cpools_df.loc[i, 'c30_SAPf']\
            + decomp_Cpools_df.loc[i, 'c30_EcM']
        decomp_Cpools_df.loc[i, 'c30_tot'] = sum30

        sum_SOM50 = decomp_Cpools_df.loc[i, 'c50_SOMa'] + \
            decomp_Cpools_df.loc[i, 'c50_SOMp'] + \
            decomp_Cpools_df.loc[i, 'c50_SOMc']
        decomp_Cpools_df.loc[i, 'c50_SOM_tot'] = sum_SOM50

        sum50 = decomp_Cpools_df.loc[i, 'c50_SOMa']+decomp_Cpools_df.loc[i, 'c50_SOMp']+decomp_Cpools_df.loc[i, 'c50_SOMc'] + \
            decomp_Cpools_df.loc[i, 'c50_LITm']+decomp_Cpools_df.loc[i, 'c50_LITs']+decomp_Cpools_df.loc[i, 'c50_SAPb'] + decomp_Cpools_df.loc[i, 'c50_SAPf']\
            + decomp_Cpools_df.loc[i, 'c50_EcM']
        decomp_Cpools_df.loc[i, 'c50_tot'] = sum50

        decomp_Cpools_df.loc[i, 'depth'] = xarray_list[i]['depth'].values
        decomp_Cpools_df.loc[i, 'f_clay'] = xarray_list[i]['f_clay'].values

    decomp_Cpools_df['SAP'] = decomp_Cpools_df.SAPb+decomp_Cpools_df.SAPf
    decomp_Cpools_df['MYC'] = decomp_Cpools_df.EcM+decomp_Cpools_df.AM

    decomp_Npools_df['N_SAP'] = decomp_Npools_df.N_SAPb+decomp_Npools_df.N_SAPf
    decomp_Npools_df['N_MYC'] = decomp_Npools_df.N_EcM+decomp_Npools_df.N_AM

    decomp_Cpools_df['FBratio'] = decomp_Cpools_df.SAPf.astype(
        float)/decomp_Cpools_df.SAPb.astype(float)
    decomp_Cpools_df['pct_microbes'] = (
        (decomp_Cpools_df.SAP+decomp_Cpools_df.EcM+decomp_Cpools_df.AM)/decomp_Cpools_df.TOTC)*100.

    decomp_Npools_df['TOTLITN'] = decomp_Npools_df.N_LITm + \
        decomp_Npools_df.N_LITs
    decomp_Npools_df['TOTSOMN'] = decomp_Npools_df.N_SOMa + \
        decomp_Npools_df.N_SOMp+decomp_Npools_df.N_SOMc

    decomp_Cpools_df['TOTLITC'] = decomp_Cpools_df.LITm+decomp_Cpools_df.LITs
    decomp_Cpools_df['TOTSOMC'] = decomp_Cpools_df.SOMa + \
        decomp_Cpools_df.SOMp+decomp_Cpools_df.SOMc

    df_tmp = decomp_Cpools_df[['LITm', 'LITs', 'SAPb', 'SAPf', 'EcM', 'AM', 'SOMp', 'SOMa', 'SOMc', 'c50_LITm', 'c50_LITs', 'c50_SAPb', 'c50_SAPf', 'c50_EcM', 'c50_AM', 'c50_SOMp', 'c50_SOMa', 'c50_SOMc',
                               'c30_LITm', 'c30_LITs', 'c30_SAPb', 'c30_SAPf', 'c30_EcM', 'c30_AM', 'c30_SOMp', 'c30_SOMa', 'c30_SOMc', 'C_PlantLITm', 'C_PlantLITs', 'C_PlantEcM', 'C_PlantAM', 'C_PlantSOMc', 'C_PlantSOMp']]
    decomp_Cpools_df.drop(['LITm', 'LITs', 'SAPb', 'SAPf', 'EcM', 'AM', 'SOMp', 'SOMa', 'SOMc', 'c30_LITm', 'c30_LITs', 'c30_SAPb', 'c30_SAPf', 'c30_EcM', 'c30_AM', 'c30_SOMp', 'c30_SOMa', 'c30_SOMc',
                           'c50_LITm', 'c50_LITs', 'c50_SAPb', 'c50_SAPf', 'c50_EcM', 'c50_AM', 'c50_SOMp', 'c50_SOMa', 'c50_SOMc', 'C_PlantLITm', 'C_PlantLITs', 'C_PlantEcM',
                           'C_PlantAM', 'C_PlantSOMc', 'C_PlantSOMp'], axis=1, inplace=True)
    #decomp_Cpools_df = decomp_Cpools_df[decomp_Cpools_df.columns[1:]].astype(float)
    #decomp_Npools_df = decomp_Npools_df[decomp_Npools_df.columns[1:]].astype(float)

    decomp_data = pd.concat([decomp_Cpools_df, decomp_Npools_df[[
                            'TOTN', "TOTLITN", "TOTSOMN", 'NH4_sol', 'NH4_sorp', 'NO3']]], axis=1)
    # decomp_data["TSA"]=CLM_data.TSA
    # decomp_data["TSOI_10CM"]=CLM_data.TSOI_10CM

    # Unit conversion:
    decomp_data.loc[:, ['TOTC', 'c30_SOM_tot', 'c30_tot', 'c50_SOM_tot', 'c50_tot', 'SAP', 'MYC', 'TOTLITC', 'TOTSOMC', 'TOTN']] = decomp_data.loc[:, [
        'TOTC', 'c30_SOM_tot', 'c30_tot', 'c50_SOM_tot', 'c50_tot',  'SAP', 'MYC', 'TOTLITC', 'TOTSOMC', 'TOTN']].astype(float)*1e-3

    decomp_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    decomp_data.loc[:, "origin"] = "decomp"
    decomp_data = decomp_data.sort_values(by="SITE")
    decomp_data.reset_index(inplace=True, drop=True)
    return decomp_data, df_tmp

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def make_CLM_df(xarray_list):

    column_names = ['RAIN', 'SNOW', 'TSA', 'TOTSOMC',
                    'TOTSOMN', 'TOTLITC', 'TOTLITN', 'SITE']
    mean_df = pd.DataFrame(columns=column_names)
    mean_df['MAP'] = np.nan
    hr = []
    # layer depths to multiply with concentrations
    depth_30 = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
    depth_50 = np.array([0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.02])  # [meters]
    sec_pr_day = 60*60*24

    for i in range(len(xarray_list)):  # loop through all xarrays (one for each site)
        val = xarray_list[i][['RAIN', 'SNOW', 'TSA', 'TOTSOMC',
                              'TOTSOMN', 'TOTLITC', 'TOTLITN', "HR", "NPP", "LITFALL"]]
        # mean over total time period for the vars in the line above
        val = val.mean(dim="time")
        # convert from xarray to pd dataframe (row)
        val = val.to_dataframe().reset_index(drop=True)
        # add site name to val (which is a dataframe row)
        val['SITE'] = [xarray_list[i].site_name]
        # add row to main dataframe
        mean_df = mean_df.append(val, ignore_index=True)

#        hr.append(val_HR.HR.values*(60*60)) #convert to /hr and add to initially empty list "hr"

        col = xarray_list[i][['SOIL1C_vr', 'SOIL2C_vr',
                              'SOIL3C_vr', 'LITR1C_vr', 'LITR2C_vr', 'LITR3C_vr']]
        col = col.mean(dim="time")  # Mean concentration for each layer
        col = col.to_dataframe().reset_index(drop=True)

        # get correct row in dataframe mean_df
        t = mean_df.loc[mean_df['SITE'] == xarray_list[i].site_name]

        for var in col.columns:  # loop variables
            # multiply with relevant depths to get mass (gC/m2)
            col30 = col[var][0:len(depth_30)]*depth_30
            col50 = col[var][0:len(depth_50)]*depth_50
            sum30 = col30.sum()  # Sum of mass in the relevant layers
            sum50 = col50.sum()
            mean_df.loc[t.index, 'c50_%s' % var[0:5]] = round(
                sum50, 2)  # add to the row where site name matches
            mean_df.loc[t.index, 'c30_%s' % var[0:5]] = round(sum30, 2)

        # Put shallowest bedrock layer (nbedrock) in dataframe mean_df
        nbedrock = xarray_list[i].nbedrock.values
        mean_df.loc[t.index, 'n_bedrock'] = nbedrock

        # Put mean annual precipitation from simulations in mean_df
        # Add RAIN and SNOW to get total precip
        daily_precip = (xarray_list[i].SNOW + xarray_list[i].RAIN)*sec_pr_day
        # Find annual precip for each year
        yearly_precip = daily_precip.groupby('time.year').sum()
        # Find mean annual precip for the years 1988-1992 (drop 1993)
        precip = yearly_precip[:-1].mean().values
        mean_df.loc[t.index, 'MAP'] = round(float(precip), 1)

        sum_soil50 = mean_df.loc[t.index, 'c50_SOIL1'] + \
            mean_df.loc[t.index, 'c50_SOIL2']+mean_df.loc[t.index, 'c50_SOIL3']
        mean_df.loc[t.index, 'c50_SOIL_tot'] = sum_soil50

        sum_soil30 = mean_df.loc[t.index, 'c30_SOIL1'] + \
            mean_df.loc[t.index, 'c30_SOIL2']+mean_df.loc[t.index, 'c30_SOIL3']
        mean_df.loc[t.index, 'c30_SOIL_tot'] = sum_soil30

        sum_tot50 = mean_df.loc[t.index, 'c50_SOIL1']+mean_df.loc[t.index, 'c50_SOIL2']+mean_df.loc[t.index, 'c50_SOIL3'] + \
            mean_df.loc[t.index, 'c50_LITR1']+mean_df.loc[t.index,
                                                          'c50_LITR2']+mean_df.loc[t.index, 'c50_LITR3']
        mean_df.loc[t.index, 'c50_tot'] = sum_tot50

        sum_tot30 = mean_df.loc[t.index, 'c30_SOIL1']+mean_df.loc[t.index, 'c30_SOIL2']+mean_df.loc[t.index, 'c30_SOIL3'] + \
            mean_df.loc[t.index, 'c30_LITR1']+mean_df.loc[t.index,
                                                          'c30_LITR2']+mean_df.loc[t.index, 'c30_LITR3']
        mean_df.loc[t.index, 'c30_tot'] = sum_tot30

    # mean_df.loc[:,'HR']=hr #Add "hr" as column to main df
    mean_df.loc[:, "origin"] = "CLM"
    # convert mean_df temperature to celsius
    mean_df.loc[:, "TSA_degC"] = round(mean_df.TSA-273.15, 1)

    mean_df["TOTN"] = mean_df.TOTSOMN+mean_df.TOTLITN
    mean_df["TOTC"] = mean_df.TOTSOMC+mean_df.TOTLITC
    mean_df["SITE"] = (mean_df.SITE).astype("category")

    for col in mean_df[['TOTSOMC', 'TOTSOMN', 'TOTLITC', 'c50_tot', 'c30_tot', 'TOTC', 'TOTN',
                        'c50_SOIL1', 'c30_SOIL1', 'c50_SOIL2', 'c30_SOIL2', 'c50_SOIL3', 'c30_SOIL3', 'c50_LITR1', 'c30_LITR1', 'c50_LITR2', 'c30_LITR2', 'c50_LITR3', 'c30_LITR3', 'c50_SOIL_tot', 'c30_SOIL_tot']].columns:  # convert mean_df content to kg/m2
        mean_df.loc[:, col] = round(mean_df.loc[:, col]*1e-3, 2)

    mean_df.loc[:, "TOTLITN"] = round(mean_df.loc[:, "TOTLITN"]*1e-3, 4)

    mean_df = mean_df.sort_values(by="SITE")
    mean_df.reset_index(inplace=True)
    mean_tmp = mean_df[['c50_SOIL1', 'c30_SOIL1', 'c50_SOIL2', 'c30_SOIL2', 'c50_SOIL3', 'c30_SOIL3', 'c50_LITR1',
                        'c30_LITR1', 'c50_LITR2', 'c30_LITR2', 'c50_LITR3', 'c30_LITR3', 'c50_SOIL_tot', 'c30_SOIL_tot', 'SITE']]
    mean_df.drop(['c50_SOIL1', 'c30_SOIL1', 'c50_SOIL2', 'c30_SOIL2', 'c50_SOIL3', 'c30_SOIL3', 'c50_LITR1', 'c30_LITR1',
                 'c50_LITR2', 'c30_LITR2', 'c50_LITR3', 'c30_LITR3', 'c50_SOIL_tot', 'c30_SOIL_tot'], axis=1, inplace=True)

    return mean_df, mean_tmp

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def no_of_layers(data):
    return len(data.levsoi)
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Convert mcdate to datetime objects from, and use this as coordinate variable instead of hours since start.


def add_date_dim(data):  # data: xarray dataset

    to_int = data.mcdate.astype('int')

    to_str = to_int.astype('str')

    dates_list = [datetime.strptime(date, '%Y%m%d')
                  for date in list(to_str.values)]

    data.coords["time"] = dates_list

    return data

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def significance_ttest(data_series1,name1, data_series2,name2, alpha):

    A = data_series1.values.tolist()
    B = data_series2.values.tolist()

    t_check = stats.ttest_ind(A, B)
    t_check
    if (t_check[1] < alpha):
        print('%s IS significantly different from %s with alpha = %1.3f' %(name1,name2,alpha))
    else:
        print('%s IS NOT significantly different from %s with alpha = %1.3f' %(name1,name2, alpha))
    print('t-statistic = %6.3f pvalue = %6.4f' % (t_check[0], t_check[1]))

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Read texture data from CLM surface file


def get_texture(surfdata):
    tex = {}
    mean_tex = {}
    tex['SAND'] = list(surfdata.PCT_SAND.values.squeeze())
    tex['CLAY'] = list(surfdata.PCT_CLAY.values.squeeze())
    silt = 100.0-surfdata.PCT_CLAY-surfdata.PCT_SAND
    tex['SILT'] = list(silt.values.squeeze())

    mean_tex['SAND'] = surfdata.PCT_SAND.values.squeeze().mean()
    mean_tex['CLAY'] = surfdata.PCT_CLAY.values.squeeze().mean()
    mean_tex['SILT'] = silt.mean()

    return tex, mean_tex
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Create a dict that contains the PFT distribution from the CLM surface files


def create_PFT_dict(surface_file):
    PFT = {}
    PFT["Site"] = surface_file.site_name
    for i in range(15):
        PFT[PFT_names[i]] = list(surface_file.PCT_NAT_PFT.values.squeeze())[0][i]
    return PFT

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Read netcdf file from decomp directory into an xarray dataset


def read(filename):
    f = xr.open_dataset('/home/ecaas/soil_decomp/results/'+filename)
    return f
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# close an xarray dataset:


def close(dataset):
    dataset.close()
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def normalized(data_array):
    norm = (data_array-data_array.min())/(data_array.max()-data_array.min())
    return norm
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# COPY THIS TO NOTEBOOK, THE globals() THING DOES NOT WORK FROM HERE.
# def import_decomp(dir_list,file_ext):
#     site_data=[]
#     decomp_names=[]
#     dirs = dir_list
#     path_of_the_directory = '/home/ecaas/soil_decomp/results/'

#     for folder in dirs:
#         ext = (folder+file_ext)
#         for files in os.listdir(path_of_the_directory+folder):
#             if files.endswith(ext):
#                 filename = os.fsdecode(files)
#                 name = filename[:-len(ext)-1]
#                 if name[0:2] != 'NR':
#                     name='NR'+name#+folder
#                 else:
#                     name = name #+folder
#                 decomp_names.append(name)
#                 globals()[name+folder] = xr.load_dataset(path_of_the_directory+folder+"/"+filename)

#                 globals()[name+folder].attrs['site_name'] = name
#                 globals()[name+folder].attrs['run_name'] = folder

#                 add_date_dim(globals()[name+folder])
#                 site_data.append(globals()[name+folder])
#             else:
#                 continue
#         site_data =  sorted(site_data,key=lambda x: x.site_name)
#         decomp_names = sorted(decomp_names)
#         #decomp_names: list of name strings
#         #site_data: list of xarrays
#     return decomp_names,site_data

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# Convert mean of C pools in xarray to dataframe
def decomp_to_dataframe(dataset):
    df = dataset[pool_names].mean(dim='time').to_dataframe()
    for pool in df:
        df[pool] = df[pool]*depth[0:len(df)]
    return df
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# sums up total mass (g/m2) in pool. returns dataArray (data = xarray)
def total_mass_pool(data, pool, depth):
    # empty arrays to fill:
    pool_mass_sum = np.zeros(len(data.time))

    for i in range(no_of_layers(data)):  # Compute mass
        # array of masses in layer i of pool for each timestep
        pool_mass = data[pool][i]*depth[i]
        pool_mass_sum += pool_mass  # sum mass
    pool_mass_sum.attrs['units'] = 'g/(m2)'
    return pool_mass_sum
    

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# sums up total mass (g/m2). 
def total_mass(data, pool_list, depth):
    # empty arrays to fill:
    pool_mass_sum = np.zeros(len(data.time))
    tot_mass_sum = np.zeros(len(data.time))
    for pool in pool_list:
        for i in range(no_of_layers(data)):  # Compute mass
            # array of masses in layer i of pool for each timestep
            pool_mass = data[pool][i]*depth[i]
            pool_mass_sum += pool_mass  # sum mass
        tot_mass_sum += pool_mass_sum  # sum mass
        pool_mass_sum = np.zeros(len(data.time))
        
    tot_mass_sum.attrs['units'] = 'g/(m2)'
    return tot_mass_sum

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# returns a time series of the total deposited N in g/m2 h


def total_depN_mass(data, depth):
    # empty array to fill:
    dep_mass_sum = np.zeros(len(data.time))

    for i in range(no_of_layers(data)):  # Compute mass
        dep_mass_layer = data['N_Deposition'][i]*depth[i]
        dep_mass_sum += dep_mass_layer

    dep_mass_sum.attrs['units'] = 'g/(m2 h)'
    return dep_mass_sum
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# returns a time series of the total leached N in g/m2 h


def total_leachN_mass(data, depth):
    # empty array to fill:
    leach_mass_sum = np.zeros(len(data.time))

    for i in range(no_of_layers(data)):  # Compute mass
        leach_mass_layer = data['N_Leaching'][i]*depth[i]
        leach_mass_sum += leach_mass_layer

    leach_mass_sum.attrs['units'] = 'g/(m2 h)'
    return leach_mass_sum


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# returns a time series of the total respired C mass from the whole soil column
def total_HR_mass(data, depth):
    # empty array to fill:
    HR_mass_sum = np.zeros(len(data.time))

    for i in range(no_of_layers(data)):  # Compute mass
        HR_mass_layer = data['HR_flux'][i]*depth[i]
        HR_mass_sum += HR_mass_layer

    HR_mass_sum.attrs['units'] = 'g/(m2 h)'
    return HR_mass_sum
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_depths_vr(data, pool_list, depth, plot_title):

    fig = plt.figure(figsize=(25, 5*no_of_layers(data)))
    fig.suptitle(plot_title, fontsize=18, y=.90)
    for i in range(no_of_layers(data)):
        ax = fig.add_subplot(no_of_layers(data), 1, i+1)
        ax.set_title("Layer nr: %i" % i, fontsize=14)
        for p in range(len(pool_list)):
            pool_mass = data[pool_list[p]][i]*depth[i]
            plt.plot(data.time, pool_mass, label='%s' % pool_list[p])
            plt.legend()
            plt.grid(True)
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# This one does not plot the total mass!


def plot_pools_vr(data, pool_list, depth, plot_title):

    fig = plt.figure(figsize=(25, 5*len(pool_list)))
    fig.suptitle(plot_title, fontsize=18, y=.92)

    for p in range(len(pool_list)):
        ax = fig.add_subplot(len(pool_list), 1, p+1)
        ax.set_title(pool_list[p], fontsize=14)
        for i in range(no_of_layers(data)):
            pool_mass = data[pool_list[p]][i]*depth[i]
            plt.plot(data.time, pool_mass, label='Layer nr: %i' % (i+1))
            plt.legend()
            plt.grid(True)


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# This one plots total mass + mass in each layer.
# Plots total mass in pools, divided by layers (gC/m2) against time (data = xarray)
def plot_pools(data, name_of_simulation, pool_list, depth):

    for pool in pool_list:
        pool_mass_sum = np.zeros(len(data.time))

        plt.figure(figsize=[20, 7])
        plt.grid(True, alpha=0.5)

        for i in range(no_of_layers(data)):
            plt.title("%s, %s " % (pool, name_of_simulation))
            pool_mass = data[pool][i]*depth[i]
            plt.plot(data.time, pool_mass, label='Layer nr: %i' % (i+1))
            plt.xlabel("Time")
            pool_mass_sum += pool_mass
            total = total_mass_pool(data, pool, depth)
        total.plot(ls='--', lw=2)

        plt.legend()
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_SAP_CUE_vr(data, depth, plot_title):

    fig = plt.figure(figsize=(25, 10))
    fig.suptitle(plot_title, fontsize=18, y=.95)
    CUE = ["CUEb", "CUEf"]
    for p in range(len(CUE)):
        ax = fig.add_subplot(len(CUE), 1, p+1)
        ax.set_title(CUE[p], fontsize=14)
        for i in range(no_of_layers(data)):
            CUE_plot = data[CUE[p]][i]
            plt.plot(data.time, CUE_plot, label='Layer nr: %i' % (i+1))
            plt.legend()
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_cons_fluxes_vr(data, flux_list, plot_title):

    fig = plt.figure(figsize=(25, 5*len(flux_list)))
    fig.suptitle(plot_title, fontsize=18, y=.99)
    for f in range(len(flux_list)):
        ax = fig.add_subplot(len(flux_list), 1, f+1)
        ax.set_title(flux_list[f], fontsize=14)

        for i in range(10,no_of_layers(data)):
            flux = data[flux_list[f]][i]
            plt.plot(data.time, flux, label='Layer nr: %i' % (i+1))
            plt.legend()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# Plot CN ration of all layers in pools against time: (data = xarray)
def plot_CN_ratio(data, name_of_simulation, pool_list):

    for pool in pool_list:

        plt.figure(figsize=[20, 7])
        plt.grid(True, alpha=0.5)

        for i in range(no_of_layers(data)):
            plt.title("C:N %s, %s" % (pool, name_of_simulation))
            pool_C = data[pool][i]
            pool_N = data['N_%s' % pool][i]
            ratio = pool_C/pool_N
            plt.plot(data.time, ratio, label='Layer nr: %i' % (i+1))
            plt.xlabel("Time")
            plt.ylabel("C:N ratio")

        plt.legend()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_HR(clm_data, decomp_data1, plot_title):
    no_of_levels = len(decomp_data1.levsoi)
    fig = plt.figure(figsize=(25, 50))
    fig.suptitle(plot_title, fontsize=18, y=.89)

    for i in range(no_of_levels):
        ax = fig.add_subplot(no_of_levels, 1, i+1)
        ax.plot(clm_data.time,
                (decomp_data1.HR_flux[i, 1:]/(60*60)), label="decomp model")
       # ax.plot(clm_data.time,(decomp_data2.HR_flux[i,1:]/(60*60)),label="HighN")

        clm_data.HR_vr[:, i].plot(label="CLM")
        plt.legend()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_pools_clm(data, name_of_simulation, pool_list, depth):

    for pool in pool_list:
        pool_mass_sum = np.zeros(len(data.time))

        plt.figure(figsize=[12, 7])
        plt.grid(True, alpha=0.5)

        for i in range(int(data.nbedrock)):
            plt.title("%s, %s " % (pool, name_of_simulation))
            t = data[pool].squeeze()
            pool_mass = t[:, i]*depth[i]
            plt.plot(data.time, pool_mass, label='Layer nr: %i' % (i+1))
            plt.xlabel("Time")
            #plt.ylabel("g /m2"%element)

            pool_mass_sum += pool_mass
            #total = total_mass_pool(data, pool,depth)
        # data.TOTSOMC.plot(label='TOTSOM')
        # data.TOTLITC.plot(label='TOTLIT')
      #  plt.plot(data.time,pool_mass_sum,'--', label="Total")
        #plt.plot(data.time,total,'--', label="Total check")

        plt.legend()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_input_and_HR_decomp(dataset):
    data = dataset.copy()

    IN_LITm = data.C_PlantLITm.copy()
    IN_LITm = IN_LITm.assign_coords(coords={'levsoi': data.levsoi})
    IN_LITs = data.C_PlantLITs.copy()
    IN_LITs = IN_LITs.assign_coords(coords={'levsoi': data.levsoi})
    IN_EcM = data.C_PlantEcM.copy()
    IN_EcM = IN_EcM.assign_coords(coords={'levsoi': data.levsoi})

    IN_AM = data.C_PlantAM.copy()
    IN_AM = IN_AM.assign_coords(coords={'levsoi': data.levsoi})
    IN_SOMp = data.C_PlantSOMp.copy()
    IN_SOMp = IN_SOMp.assign_coords(coords={'levsoi': data.levsoi})
    IN_SOMa = data.C_PlantSOMa.copy()
    IN_SOMa = IN_SOMa.assign_coords(coords={'levsoi': data.levsoi})
    IN_SOMc = data.C_PlantSOMc.copy()
    IN_SOMc = IN_SOMc.assign_coords(coords={'levsoi': data.levsoi})
    for i in range(len(data.time)):
        IN_LITm[:, i] = IN_LITm[:, i]*depth[0:len(data.levsoi)]
        IN_LITs[:, i] = IN_LITs[:, i]*depth[0:len(data.levsoi)]
        IN_EcM[:, i] = IN_EcM[:, i]*depth[0:len(data.levsoi)]
        IN_AM[:, i] = IN_AM[:, i]*depth[0:len(data.levsoi)]
        IN_SOMa[:, i] = IN_SOMa[:, i]*depth[0:len(data.levsoi)]
        IN_SOMp[:, i] = IN_SOMp[:, i]*depth[0:len(data.levsoi)]
        IN_SOMc[:, i] = IN_SOMc[:, i]*depth[0:len(data.levsoi)]

    HR = total_HR_mass(data, depth)
    tot_in = IN_LITm.sum(axis=0)+IN_LITs.sum(axis=0) + IN_EcM.sum(axis=0)+IN_AM.sum(
        axis=0)+IN_SOMa.sum(axis=0)+IN_SOMp.sum(axis=0)+IN_SOMc.sum(axis=0)
    plt.figure(figsize=[15, 7])
    (tot_in).plot(label="input")
    HR.plot(label="HR")
    plt.legend()
    print("IN  : ", round(tot_in.mean().values*(24*365*5)/1000, 2), "kgC/m2")
    print("HR  : ", round(HR.mean().values*(24*365*5)/1000, 2), "kgC/m2")
    print("Diff: ", round((tot_in.mean().values-HR.mean().values)*(24*365*5)/1000, 2))
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_input_and_HR(data):  # For CLM data
    CWD_TO_LIT3 = (data.CWDC_TO_LITR3C_vr.squeeze()*depth).sum(axis=1)
    CWD_TO_LIT2 = (data.CWDC_TO_LITR2C_vr.squeeze()*depth).sum(axis=1)
    LEAFC_TO_LITTER = data.LEAFC_TO_LITTER
    FROOTC_TO_LITTER = data.FROOTC_TO_LITTER

    plt.figure(figsize=[15, 7])
    tot_in = CWD_TO_LIT3 + CWD_TO_LIT2 + LEAFC_TO_LITTER + \
        FROOTC_TO_LITTER  # +data.LITFALL# CWD_TO_LIT3 + CWD_TO_LIT2
    # +data.LITFALL# CWD_TO_LIT3 + CWD_TO_LIT2
    tot_in_LITFALL = CWD_TO_LIT3 + CWD_TO_LIT2 + data.LITFALL
    tot_in_LITFALL.plot(label="LITFALL")

    tot_in.plot(label="input")

    data.HR.plot(label="HR")
    print("IN LITFALL   : ", round(
        tot_in_LITFALL.mean().values*(60*60*24*365*5)/1000, 3), "kgC/m2")
    print("- HR         : ", round(data.HR.mean().values *
          (60*60*24*365*5)/1000, 3), "kgC/m2")
    print("= Diff       : ", round((tot_in_LITFALL.mean().values -
          data.HR.mean().values)*(60*60*24*365*5)/1000, 3))
    print("---------------------------------------")

    print("IN LEAF&FROOT: ", round(
        tot_in.mean().values*(60*60*24*365*5)/1000, 3), "kgC/m2")
    print("- HR         : ", round(data.HR.mean().values *
          (60*60*24*365*5)/1000, 3), "kgC/m2")
    print("= Diff       : ", round((tot_in.mean().values -
          data.HR.mean().values)*(60*60*24*365*5)/1000, 3))
    plt.legend()


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def plot_net_mineralization(data):
    net_mineralization = -(data.N_INSAPb+data.N_INSAPf)
    fig, ax = plt.subplots(figsize=(30, 7))
    net_mineralization.plot(yincrease=False)
    plt.title("Net mineralization (SAP release of N) when positive, immobilization (SAP uptake of N) when negative")

    # Use to convert from gC/m3 to gC/m2
    depths = xr.DataArray(depth[:no_of_layers(data)], dims="levsoi")
    t_mass = net_mineralization*depths
    t_mass_tot = t_mass.sum(axis=0)  # Sum over all layers
    fig, ax = plt.subplots(figsize=(30, 7))
    t_mass_tot.plot()
    (xr.zeros_like(t_mass_tot)).plot(color="black")
    t_mass.attrs['units'] = 'g/(m2 h)'
    plt.title("Total net mineralization")
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_correlation_matrix(dataframe, cut_off, method):
    dataframe = dataframe.corr(method=method)
    mask = np.triu(np.ones_like(dataframe, dtype=bool))

    mask |= np.abs(dataframe) < cut_off
    dataframe_copy = dataframe[~mask]  # fill in NaN in the non-desired cells
    nr_of_vars = len(dataframe)

    if method == "spearman":
        pval = dataframe.corr(method=lambda x, y: spearmanr(x, y)[
                              1]) - np.eye(*dataframe.shape)

    elif method == "pearson":
        pval = dataframe.corr(method=lambda x, y: pearsonr(x, y)[
                              1]) - np.eye(*dataframe.shape)

    else:
        print("Invalid method")

    p = pval.applymap(lambda x: ''.join(
        ['*' for t in [0.01, 0.05, 0.1] if x <= t]))

    with sns.axes_style("dark"):
        f, ax = plt.subplots(figsize=(0.9*nr_of_vars, 0.7*nr_of_vars))
        ax = sns.heatmap(dataframe_copy, vmin=-1, vmax=1, annot=dataframe.round(2).astype(
            str) + p, fmt='', annot_kws={"size": 10}, mask=mask, linewidth=.3, cmap="coolwarm")
        plt.yticks(rotation=0)
    ax.set_title("%s correlation heatmap" % method)
    textstr = '\n'.join((
        r'***  : $\alpha$ = 0.01',
        r'**   : $\alpha$ = 0.05',
        r'*    : $\alpha$ = 0.1'))

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(0.75, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_mass_and_conc_carbon(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")

    # Use to convert from gC/m3 to gC/m2
    depths = xr.DataArray(depth[:no_of_layers(data)], dims="levsoi")
    weighted_season_mean_mass = weighted_season_mean*depths

    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("Mass of C, "+title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(pool_names):
        ax = plt.subplot(5, 3, n + 1)
        ax.invert_yaxis()
        plt.xlabel("[gC/m2]")
        plt.ylabel("Depth [m]")
        plt.grid()
        for i in range(4):
            ax.plot(weighted_season_mean_mass[pool_names[n]][i, :], node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (pool_names[n]))
            ax.legend()

    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("Concentration of C, "+title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(pool_names):
        ax = plt.subplot(5, 3, n + 1)
        ax.invert_yaxis()
        plt.xlabel("[gC/m3]")
        plt.ylabel(" [m]")
        plt.grid()
        for i in range(4):
            ax.plot(weighted_season_mean[pool_names[n]][i, :], node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (pool_names[n]))
            ax.legend()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_CN(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")

    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("C:N ratio, "+title, y=0.92, fontsize="x-large")
    plt.grid()

    for n, ticker in enumerate(pool_names):
        ax = plt.subplot(5, 3, n + 1)
        ax.invert_yaxis()
        plt.ylabel(" [m]")

        for i in range(4):

            ax.plot((weighted_season_mean[pool_names[n]][i, :]/weighted_season_mean["N_"+pool_names[n]][i, :]),
                    node_z[0:no_of_layers(weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (pool_names[n]))
            ax.legend()
            plt.grid()
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_mass_and_conc_nitrogen(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")

    # Use to convert from gC/m3 to gC/m2
    depths = xr.DataArray(depth[:no_of_layers(data)], dims="levsoi")
    weighted_season_mean_mass = weighted_season_mean*depths

    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("Mass of N, "+title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(N_pool_names):
        ax = plt.subplot(5, 3, n + 1)
        ax.invert_yaxis()
        plt.xlabel("[gN/m2]")
        plt.ylabel("Depth [m]")
        plt.grid(True)

        for i in range(4):
            ax.plot(weighted_season_mean_mass[N_pool_names[n]][i, :], node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (N_pool_names[n]))
            ax.legend()
            plt.grid(True)

    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("Concentration of N, "+title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(N_pool_names):
        ax = plt.subplot(5, 3, n + 1)
        ax.invert_yaxis()
        plt.xlabel("[gN/m3]")
        plt.ylabel("Depth [m]")
        plt.grid(True)
        plt.grid(True)

        for i in range(4):

            ax.plot(weighted_season_mean[N_pool_names[n]][i, :], node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (N_pool_names[n]))
            ax.legend()
            plt.grid(True)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_CUE(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")
    var = ["CUEb", "CUEf", "CUE_ecm", "CUE_am"]
    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("CUE for " + title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(var):
        ax = plt.subplot(5, 2, n + 1)
        ax.invert_yaxis()
        for i in range(4):

            ax.plot((weighted_season_mean[var[n]][i, :]), node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (var[n]))
            ax.legend()
            plt.grid(True)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_HR(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")
    var = ["HR_flux", "HRb", "HRf", "HRe", "HRa"]
    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("HR for "+title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(var):
        ax = plt.subplot(5, 3, n + 1)
        ax.invert_yaxis()
        plt.xlabel("[gC/(m3 s)]")
        plt.ylabel(" [m]")
        for i in range(4):

            ax.plot((weighted_season_mean[var[n]][i, :]/(60*60)), node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (var[n]))
            # ax.set_xscale('log')
            ax.legend()
            ax.grid(True)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_ROI(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")
    var = ["ROI_ecm", "ROI_am", "f_ecm", "f_am"]
    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("ROI for "+title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(var):
        ax = plt.subplot(5, 2, n + 1)
        ax.invert_yaxis()
        plt.xlabel("Return Of Investment function")
        plt.ylabel("Depth [m]")
        for i in range(4):

            ax.plot((weighted_season_mean[var[n]][i, :]), node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (var[n]))
            ax.grid(True)
            ax.legend()
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_immobilization(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")
    var = ["N_INSAPf", "N_INSAPb"]
    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("Immobilization/Mineralization" +
                 title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(var):
        ax = plt.subplot(5, 2, n + 1)
        ax.invert_yaxis()
        plt.xlabel("mean gN/m3h")
        plt.ylabel("Depth [m]")
        for i in range(4):

            ax.plot((weighted_season_mean[var[n]][i, :]), node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (var[n]))
            ax.legend()
            plt.grid(True)
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_seasonal_climate(data, title):
    month_length = data.time.dt.days_in_month
    weights = (month_length.groupby("time.season") /
               month_length.groupby("time.season").sum())
    weighted_season_mean = (
        weights*data).groupby("time.season").sum(dim="time")
    var = ["Temp", "Moisture"]  # , "N_Deposition","N_Leaching"]
    fig = plt.figure(figsize=[20, 25])
    fig.suptitle(title, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(var):
        ax = plt.subplot(5, 2, n + 1)
        ax.invert_yaxis()
        for i in range(4):

            ax.plot((weighted_season_mean[var[n]][i, :]), node_z[0:no_of_layers(
                weighted_season_mean)], '-*', label="%s" % (weighted_season_mean.season[i].values))
            ax.set_title("%s" % (var[n]))
            ax.grid(True)
            ax.legend()

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_vertical_change(data):
    var = []
    for pool in pool_names:
        var.append("vert_change"+pool)
    fig = plt.figure(figsize=[20, 30])
    fig.suptitle("Vertical changes", y=0.92, fontsize="x-large")
    for n, ticker in enumerate(var):
        ax = plt.subplot(5, 3, n + 1)
        ax.invert_yaxis()

        # (data[var[n]][i,:]),node_z[0:no_of_layers(data)],'-*')
        data[var[n]].plot(label=var[n])

        # ax.pcolormesh(data.time,data[var[n]],data)
        ax.set_title("%s" % (var[n]))
        # plt.legend()


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Note: needs conc_ variables in dataframes!
def plot_profiles_clim_mineral(obs, decomp, CLM, title, save):
    fig, ax = plt.subplots(figsize=(10, 10))

    d = [0.3, 0.5]
    d3 = [0.3 for i in range(10)]
    d5 = [0.5 for i in range(10)]
    #x = np.linspace(0,1,10)
    plt.ylabel("Depth below surface [m]")
    plt.xlabel("Concentration [kgC/m3]")
    plt.title("Mineral soil only, " + title)
    plt.ylim(0, 1)
    plt.grid(True)

    ax.invert_yaxis()

    ax.axhspan(0.0, 0.3, facecolor='gray', alpha=0.1)
    ax.axhspan(0.3, 0.5, facecolor='gray', alpha=0.2)
    ax.axhspan(0.5, 1.0, facecolor='gray', alpha=0.3)

    ax.errorbar(obs.conc_min30.mean(), [
                0.13], xerr=obs.conc_min30.sem(), marker="D", color="green")
    ax.errorbar(obs.conc_min30.median(), [0.13], marker="*", color="green")

    ax.errorbar(CLM.conc_min30.mean(), [
                0.15], xerr=CLM.conc_min30.sem(), marker="D", color="blue")
    ax.errorbar(CLM.conc_min30.median(), [0.15], marker="*", color="blue")

    ax.errorbar(decomp.conc_min30.mean(), [
                0.17], xerr=decomp.conc_min30.sem(), marker="D", color="orange")
    ax.errorbar(decomp.conc_min30.median(), [0.17], marker="*", color="orange")

    ax.errorbar(obs.conc_min30_50.mean(), [
                0.38], xerr=obs.conc_min30_50.sem(), marker="D", color="green")
    ax.errorbar(obs.conc_min30_50.median(), [0.38], marker="*", color="green")

    ax.errorbar(CLM.conc_min30_50.mean(), [
                0.4], xerr=CLM.conc_min30_50.sem(), marker="D", color="blue")
    ax.errorbar(CLM.conc_min30_50.median(), [0.4], marker="*", color="blue")

    ax.errorbar(decomp.conc_min30_50.mean(), [
                0.42], xerr=decomp.conc_min30_50.sem(), marker="D", color="orange")
    ax.errorbar(decomp.conc_min30_50.median(), [
                0.42], marker="*", color="orange")

    ax.errorbar(obs.conc_min50_100.mean(), [
                0.78], xerr=obs.conc_min50_100.sem(), marker="D", color="green")
    ax.errorbar(obs.conc_min50_100.median(), [
                0.78], marker="*", color="green", label="obs, %i" % len(obs))

    ax.errorbar(CLM.conc_min50_100.mean(), [
                0.82], xerr=CLM.conc_min50_100.sem(), marker="D", color="blue")
    ax.errorbar(CLM.conc_min50_100.median(), [
                0.82], marker="*", color="blue", label="CLM, %i" % len(CLM))

    ax.errorbar(decomp.conc_min50_100.mean(), [
                0.86], xerr=decomp.conc_min50_100.sem(), marker="D", color="orange")
    ax.errorbar(decomp.conc_min50_100.median(), [
                0.86], marker="*", color="orange", label="decomp, %i" % len(decomp))
    ax.legend()
    if save:
        plt.savefig("profile_mineralsoil"+title+(date.today()
                                                 ).strftime("%d%m%y")+"_sem.png", bbox_inches="tight")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# Note: needs conc_ variables in dataframes!
def plot_profiles_clim(obs, decomp, CLM, title, save):
    fig, ax = plt.subplots(figsize=(10, 10))

    d = [0.3, 0.5]
    d3 = [0.3 for i in range(10)]
    d5 = [0.5 for i in range(10)]
    #x = np.linspace(0,1,10)
    plt.ylabel("Depth below surface [m]")
    plt.xlabel("Concentration [kgC/m3]")
    plt.title(title)
    plt.ylim(0, 1)
    plt.grid(True)

    ax.invert_yaxis()

    ax.axhspan(0.0, 0.3, facecolor='gray', alpha=0.1)
    ax.axhspan(0.3, 0.5, facecolor='gray', alpha=0.2)
    ax.axhspan(0.5, 1.0, facecolor='gray', alpha=0.3)

    ax.errorbar(obs.conc30.mean(), [
                0.13], xerr=obs.conc30.sem(), marker="D", color="green")
    ax.errorbar(obs.conc30.median(), [0.13], marker="*", color="green")

    ax.errorbar(CLM.conc30.mean(), [
                0.15], xerr=CLM.conc30.sem(), marker="D", color="blue")
    ax.errorbar(CLM.conc30.median(), [0.15], marker="*", color="blue")

    ax.errorbar(decomp.conc30.mean(), [
                0.17], xerr=decomp.conc30.sem(), marker="D", color="orange")
    ax.errorbar(decomp.conc30.median(), [0.17], marker="*", color="orange")

    ax.errorbar(obs.conc30_50.mean(), [
                0.38], xerr=obs.conc30_50.sem(), marker="D", color="green")
    ax.errorbar(obs.conc30_50.median(), [0.38], marker="*", color="green")

    ax.errorbar(CLM.conc30_50.mean(), [
                0.4], xerr=CLM.conc30_50.sem(), marker="D", color="blue")
    ax.errorbar(CLM.conc30_50.median(), [0.4], marker="*", color="blue")

    ax.errorbar(decomp.conc30_50.mean(), [
                0.42], xerr=decomp.conc30_50.sem(), marker="D", color="orange")
    ax.errorbar(decomp.conc30_50.median(), [0.42], marker="*", color="orange")

    ax.errorbar(obs.conc50_100.mean(), [
                0.78], xerr=obs.conc50_100.sem(), marker="D", color="green")
    ax.errorbar(obs.conc50_100.median(), [
                0.78], marker="*", color="green", label="obs, %i" % len(obs))

    ax.errorbar(CLM.conc50_100.mean(), [
                0.82], xerr=CLM.conc50_100.sem(), marker="D", color="blue")
    ax.errorbar(CLM.conc50_100.median(), [
                0.82], marker="*", color="blue", label="CLM, %i" % len(CLM))

    ax.errorbar(decomp.conc50_100.mean(), [
                0.86], xerr=decomp.conc50_100.sem(), marker="D", color="orange")
    ax.errorbar(decomp.conc50_100.median(), [
                0.86], marker="*", color="orange", label="decomp, %i" % len(decomp))
    ax.legend()
    if save:
        plt.savefig("profile_"+title+(date.today()
                                      ).strftime("%d%m%y")+"_sem.png", bbox_inches="tight")
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_profiles_clim_same(min_df, middle_df, max_df, clim, dataset, color_list, save):
    fig, ax = plt.subplots(figsize=(10, 10))

    d = [0.3, 0.5]
    d3 = [0.3 for i in range(10)]
    d5 = [0.5 for i in range(10)]
    #x = np.linspace(0,1,10)
    plt.ylabel("Depth below surface [m]")
    plt.xlabel("Concentration [kgC/m3]")
    plt.title(clim+dataset)
    plt.ylim(0, 1)
    plt.grid(True)

    ax.invert_yaxis()

    ax.axhspan(0.0, 0.3, facecolor='gray', alpha=0.1)
    ax.axhspan(0.3, 0.5, facecolor='gray', alpha=0.2)
    ax.axhspan(0.5, 1.0, facecolor='gray', alpha=0.3)

    ax.errorbar(min_df.conc30.mean(), [
                0.13], xerr=min_df.conc30.sem(), marker="D", color=color_list[0])
    ax.errorbar(min_df.conc30.median(), [
                0.13], marker="*", color=color_list[0])

    ax.errorbar(middle_df.conc30.mean(), [
                0.15], xerr=middle_df.conc30.sem(), marker="D", color=color_list[1])
    ax.errorbar(middle_df.conc30.median(), [
                0.15], marker="*", color=color_list[1])

    ax.errorbar(max_df.conc30.mean(), [
                0.17], xerr=max_df.conc30.sem(), marker="D", color=color_list[2])
    ax.errorbar(max_df.conc30.median(), [
                0.17], marker="*", color=color_list[2])

    ax.errorbar(min_df.conc30_50.mean(), [
                0.38], xerr=min_df.conc30_50.sem(), marker="D", color=color_list[0])
    ax.errorbar(min_df.conc30_50.median(), [
                0.38], marker="*", color=color_list[0])

    ax.errorbar(middle_df.conc30_50.mean(), [
                0.4], xerr=middle_df.conc30_50.sem(), marker="D", color=color_list[1])
    ax.errorbar(middle_df.conc30_50.median(), [
                0.4], marker="*", color=color_list[1])

    ax.errorbar(max_df.conc30_50.mean(), [
                0.42], xerr=max_df.conc30_50.sem(), marker="D", color=color_list[2])
    ax.errorbar(max_df.conc30_50.median(), [
                0.42], marker="*", color=color_list[2])

    ax.errorbar(min_df.conc50_100.mean(), [
                0.78], xerr=min_df.conc50_100.sem(), marker="D", color=color_list[0])
    ax.errorbar(min_df.conc50_100.median(), [0.78], marker="*", color=color_list[0],
                label="%s, %i" % ((min_df.reset_index()).loc[0, clim], len(min_df)))

    ax.errorbar(middle_df.conc50_100.mean(), [
                0.80], xerr=middle_df.conc50_100.sem(), marker="D", color=color_list[1])
    ax.errorbar(middle_df.conc50_100.median(), [0.80], marker="*", color=color_list[1],
                label="%s, %i" % ((middle_df.reset_index()).loc[0, clim], len(middle_df)))

    ax.errorbar(max_df.conc50_100.mean(), [
                0.82], xerr=max_df.conc50_100.sem(), marker="D", color=color_list[2])
    ax.errorbar(max_df.conc50_100.median(), [0.82], marker="*", color=color_list[2],
                label="%s, %i" % ((max_df.reset_index()).loc[0, clim], len(max_df)))
    ax.legend()
    if save:
        plt.savefig(dataset+"_"+clim+(date.today()
                                      ).strftime("%d%m%y")+"_sem.png", bbox_inches="tight")
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_spunup_mean(data, name_of_sim, log_bool):

    fig = plt.figure(figsize=[20, 28])
    fig.suptitle("Concentrations "+name_of_sim, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(pool_names):
        # Mean of last 20 years of spinup
        data20 = (data[pool_names[n]][:, -20:]).mean(axis=1)

        ax = plt.subplot(5, 4, n + 1)
        ax.invert_yaxis()
        plt.xlabel("[gC/m3]", fontsize="x-small")
        plt.ylabel("Depth [m]", fontsize="x-small")
        plt.grid()
        ax.plot(data20, node_z[0:no_of_layers(data20)], '-*')

        if log_bool:
            ax.set_xscale('log')
        ax.set_title("%s" % (pool_names[n]), fontsize="x-small")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_spunup_mean_CN(data, name_of_sim):
    fig = plt.figure(figsize=[20, 28])
    fig.suptitle("CN "+name_of_sim, y=0.92, fontsize="x-large")
    for n, ticker in enumerate(pool_names):
        # Mean of last 20 years of spinup
        data20_C = (data[pool_names[n]][:, -20:]).mean(axis=1)
        # Mean of last 20 years of spinup
        data20_N = (data[N_pool_names[n]][:, -20:]).mean(axis=1)
        CN = data20_C/data20_N
        ax = plt.subplot(5, 4, n + 1)
        ax.invert_yaxis()
        plt.xlabel("[-]", fontsize="x-small")
        plt.ylabel("Depth [m]")
        plt.grid()
        ax.plot(CN, node_z[0:no_of_layers(data)], '-*')
        ax.set_title("%s" % (pool_names[n]), fontsize="x-small")

    # Mean of last 20 years of spinup
    data20_C_LITm = (data[pool_names[0]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_LITs = (data[pool_names[1]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMp = (data[pool_names[6]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMa = (data[pool_names[7]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMc = (data[pool_names[8]][:, -20:]).mean(axis=1)

    # Mean of last 20 years of spinup
    data20_N_LITm = (data[N_pool_names[0]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_LITs = (data[N_pool_names[1]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMp = (data[N_pool_names[6]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMa = (data[N_pool_names[7]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMc = (data[N_pool_names[8]][:, -20:]).mean(axis=1)

    CN_SOM = (data20_C_SOMp+data20_C_SOMc+data20_C_SOMa) / \
        (data20_N_SOMp+data20_N_SOMc+data20_N_SOMa)
    CN_LIT = (data20_C_LITm + data20_C_LITs)/(data20_N_LITm+data20_N_LITs)
    CN_SOM_LIT = (data20_C_SOMp+data20_C_SOMc+data20_C_SOMa+data20_C_LITm+data20_C_LITs) / \
        (data20_N_SOMp+data20_N_SOMc+data20_N_SOMa+data20_N_LITm+data20_N_LITs)

    ax = plt.subplot(5, 4, n + 2)
    ax.invert_yaxis()
    plt.xlabel("[-]", fontsize="x-small")
    plt.ylabel("Depth [m]")
    plt.grid()
    ax.plot(CN_LIT, node_z[0:no_of_layers(data)], '-*')
    ax.set_title("%s" % (pool_names[0:2]), fontsize="x-small")

    ax = plt.subplot(5, 4, n + 3)
    ax.invert_yaxis()
    plt.xlabel("[-]", fontsize="x-small")
    plt.ylabel("Depth [m]")
    plt.grid()
    ax.plot(CN_SOM, node_z[0:no_of_layers(data)], '-*')
    ax.set_title("%s" % (pool_names[6:9]), fontsize="x-small")

    ax = plt.subplot(5, 4, n + 4)
    ax.invert_yaxis()
    plt.xlabel("[-]", fontsize="x-small")
    plt.ylabel("Depth [m]")
    plt.grid()
    ax.plot(CN_SOM_LIT, node_z[0:no_of_layers(data)], '-*')
    ax.set_title("%s" % (pool_names[0:2]+pool_names[6:9]), fontsize="x-small")
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_cons_vr(data, flux_list, plot_title):

    fig = plt.figure(figsize=(25, 5*len(flux_list)))
    fig.suptitle(plot_title, fontsize=18, y=.99)
    for f in range(len(flux_list)):
        ax = fig.add_subplot(len(flux_list), 1, f+1)
        ax.set_title(flux_list[f], fontsize=14)
        # plt.xlim(4500,5000)
        for i in range(7, 9):
            flux = data[flux_list[f]][i]
            plt.plot(data.time, flux, label='Layer nr: %i' % (i+1))
            plt.legend()


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def plot_pool_fractions(data, name, fig, ax):

    # Mean of last 20 years of spinup
    data20_C_LITm = (data["LITm"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_LITs = (data["LITs"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMp = (data["SOMp"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMa = (data["SOMa"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMc = (data["SOMc"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SAPb = (data["SAPb"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SAPf = (data["SAPf"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_EcM = (data["EcM"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_AM = (data["AM"][:, -20:]).mean(axis=1)

    mass_sum_LITm = (data20_C_LITm*depth[0:no_of_layers(data)]).sum()
    mass_sum_LITs = (data20_C_LITs*depth[0:no_of_layers(data)]).sum()
    mass_sum_SOMp = (data20_C_SOMp*depth[0:no_of_layers(data)]).sum()
    mass_sum_SOMa = (data20_C_SOMa*depth[0:no_of_layers(data)]).sum()
    mass_sum_SOMc = (data20_C_SOMc*depth[0:no_of_layers(data)]).sum()
    mass_sum_SAPb = (data20_C_SAPb*depth[0:no_of_layers(data)]).sum()
    mass_sum_SAPf = (data20_C_SAPf*depth[0:no_of_layers(data)]).sum()
    mass_sum_EcM = (data20_C_EcM*depth[0:no_of_layers(data)]).sum()
    mass_sum_AM = (data20_C_AM*depth[0:no_of_layers(data)]).sum()

    total_mass = mass_sum_LITm + mass_sum_LITs + mass_sum_SOMp + mass_sum_SOMc + \
        mass_sum_SOMa + mass_sum_EcM+mass_sum_AM+mass_sum_SAPf+mass_sum_SAPb

    # Mean of last 20 years of spinup
    data20_N_LITm = (data[N_pool_names[0]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_LITs = (data[N_pool_names[1]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMp = (data[N_pool_names[6]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMa = (data[N_pool_names[7]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMc = (data[N_pool_names[8]][:, -20:]).mean(axis=1)

    data20_sumC = data20_C_SOMp+data20_C_SOMc+data20_C_SOMa+data20_C_LITm + \
        data20_C_LITs+data20_C_SAPb+data20_C_SAPf+data20_C_AM+data20_C_EcM
    ax.invert_yaxis()
    ax.plot((data20_C_LITm/data20_sumC),
            node_z[0:len(data.levsoi)], label="LITm fraction", marker="o")
    ax.plot((data20_C_LITs/data20_sumC),
            node_z[0:len(data.levsoi)], label="LITs fraction", marker="o")
    ax.plot((data20_C_SOMa/data20_sumC),
            node_z[0:len(data.levsoi)], label="SOMa fraction", marker="o")
    ax.plot((data20_C_SOMc/data20_sumC),
            node_z[0:len(data.levsoi)], label="SOMc fraction", marker="o")
    ax.plot((data20_C_SOMp/data20_sumC),
            node_z[0:len(data.levsoi)], label="SOMp fraction", marker="o")
    ax.plot((data20_C_SAPb/data20_sumC),
            node_z[0:len(data.levsoi)], label="SAPb fraction", marker="o")
    ax.plot((data20_C_SAPf/data20_sumC),
            node_z[0:len(data.levsoi)], label="SAPf fraction", marker="o")
    ax.plot((data20_C_EcM/data20_sumC),
            node_z[0:len(data.levsoi)], label="EcM fraction", marker="o")
    ax.plot((data20_C_AM/data20_sumC),
            node_z[0:len(data.levsoi)], label="AM fraction", marker="o")
    ax.legend()
    ax.set_title(data.site_name)
    ax.set_xlabel("Fraction of total C content")
    ax.set_ylabel("Depth [m]")


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def plot_check_sites(data, name, fig, ax_list):

    # Mean of last 20 years of spinup
    data20_C_LITm = (data["LITm"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_LITs = (data["LITs"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMp = (data["SOMp"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMa = (data["SOMa"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SOMc = (data["SOMc"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SAPb = (data["SAPb"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_SAPf = (data["SAPf"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_EcM = (data["EcM"][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_C_AM = (data["AM"][:, -20:]).mean(axis=1)

    mass_sum_LITm = (data20_C_LITm*depth[0:no_of_layers(data)]).sum()
    mass_sum_LITs = (data20_C_LITs*depth[0:no_of_layers(data)]).sum()
    mass_sum_SOMp = (data20_C_SOMp*depth[0:no_of_layers(data)]).sum()
    mass_sum_SOMa = (data20_C_SOMa*depth[0:no_of_layers(data)]).sum()
    mass_sum_SOMc = (data20_C_SOMc*depth[0:no_of_layers(data)]).sum()
    mass_sum_SAPb = (data20_C_SAPb*depth[0:no_of_layers(data)]).sum()
    mass_sum_SAPf = (data20_C_SAPf*depth[0:no_of_layers(data)]).sum()
    mass_sum_EcM = (data20_C_EcM*depth[0:no_of_layers(data)]).sum()
    mass_sum_AM = (data20_C_AM*depth[0:no_of_layers(data)]).sum()

    total_mass = mass_sum_LITm + mass_sum_LITs + mass_sum_SOMp + mass_sum_SOMc + \
        mass_sum_SOMa + mass_sum_EcM+mass_sum_AM+mass_sum_SAPf+mass_sum_SAPb

    # Mean of last 20 years of spinup
    data20_N_LITm = (data[N_pool_names[0]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_LITs = (data[N_pool_names[1]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMp = (data[N_pool_names[6]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMa = (data[N_pool_names[7]][:, -20:]).mean(axis=1)
    # Mean of last 20 years of spinup
    data20_N_SOMc = (data[N_pool_names[8]][:, -20:]).mean(axis=1)

    data20_sumC = data20_C_SOMp+data20_C_SOMc+data20_C_SOMa+data20_C_LITm + \
        data20_C_LITs+data20_C_SAPb+data20_C_SAPf+data20_C_AM+data20_C_EcM

    # plotting:

    ax_list[0][0].plot((data20_C_SAPf/data20_C_SAPb),
                       node_z[0:len(data.levsoi)], label=name, marker="o")
    (ax_list[0][0]).set_title('F:B ratio (based on concentration)')
    (ax_list[0][0]).legend()

    ax_list[1][0].plot(((data20_C_LITm+data20_C_LITs)/data20_sumC),
                       node_z[0:len(data.levsoi)], label=name, marker="o")
    (ax_list[1][0]).set_title("Litter fraction")
    (ax_list[1][0]).legend()  # ;(ax_list[1][0]).invert_yaxis()

    ax_list[2][0].plot((data20_C_SOMp+data20_C_SOMc+data20_C_SOMa) /
                       data20_sumC, node_z[0:len(data.levsoi)], label=name, marker="o")
    (ax_list[2][0]).set_title("SOM fraction")
    (ax_list[2][0]).legend()  # ;(ax_list[2][0]).invert_yaxis()

    ax_list[3][0].plot(((data20_C_SAPb+data20_C_SAPf+data20_C_EcM+data20_C_AM) /
                       data20_sumC), node_z[0:len(data.levsoi)], label=name, marker="o")
    (ax_list[3][0]).set_title("Microbial fraction")
    (ax_list[3][0]).legend()  # ;(ax_list[3][0]).invert_yaxis()

    ax_list[4][0].plot(((data20_C_SAPb+data20_C_SAPf)/data20_sumC),
                       node_z[0:len(data.levsoi)], label=name, marker="o")
    (ax_list[4][0]).set_title("Saprotroph fraction")
    (ax_list[4][0]).legend()  # ;(ax_list[4][0]).invert_yaxis()

    ax_list[5][0].plot(((data20_C_EcM+data20_C_AM)/data20_sumC),
                       node_z[0:len(data.levsoi)], label=name, marker="o")
    (ax_list[5][0]).set_title("Mycorrhizal fraction")
    (ax_list[5][0]).legend()  # ;(ax_list[5][0]).invert_yaxis()

# use example:'
# fig,  axes = plt.subplots(nrows=6, ncols=1,figsize = [12,50],squeeze=False)
# for i in range(6):
#     (axes[i][0]).invert_yaxis()

# check_sites(Byg_obs_spinup,"Byg_obs_spinup",fig,axes)


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def C_content_bar_plot(data, placement, color, label, fig, ax):

    mass_dict = Ctotal_pools(data)
    # Shape :
    # {'LITm': 849.26,
    # 'LITs': 1864.81,
    # 'SAPb': 190.69,
    # 'SAPf': 190.69,
    # 'EcM': 16.24,
    # 'SOMp': 1708.2,
    # 'SOMa': 2624.64,
    # 'SOMc': 3243.11}

    labels = list(mass_dict.keys())
    x = np.arange(len(labels))  # the label locations
    width = 0.1

    rects1 = ax.bar(x+placement*width, mass_dict.values(),
                    width, color=color, label=label)

    plt.title("Total carbon content by pool", size="x-large")
    plt.legend(fontsize="x-large")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, size="x-large")

    ax.bar_label(rects1, padding=3, size="small", fmt='%.1f', rotation=45)

    ax.set_ylabel('gC/m2', size="x-large")
    plt.rcParams['xtick.bottom'] = True

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def C30_content_bar_plot(data, placement, color, label, fig, ax):

    mass_dict = C30_pools(data)
    # Shape :
    # {'LITm': 849.26,
    # 'LITs': 1864.81,
    # 'SAPb': 190.69,
    # 'SAPf': 190.69,
    # 'EcM': 16.24,
    # 'SOMp': 1708.2,
    # 'SOMa': 2624.64,
    # 'SOMc': 3243.11}

    labels = list(mass_dict.keys())
    x = np.arange(len(labels))  # the label locations
    width = 0.15

    rects1 = ax.bar(x+placement*width, mass_dict.values(),
                    width, color=color, label=label)

    plt.title("C content down to 30cm by pool", size="x-large")
    plt.legend(fontsize="x-large")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, size="x-large")

    ax.bar_label(rects1, padding=3, size="small", fmt='%.1f', rotation=45)

    ax.set_ylabel('gC/m2', size="x-large")
    plt.rcParams['xtick.bottom'] = True
    # ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def Ninorg_content_bar_plot(data, placement, color, label, fig, ax):

    mass_dict = Ntotal_pools(data)

    inorg = {k: mass_dict[k] for k in ('NH4_sorp', 'NO3', 'NH4_sol')}
    labels = list(inorg.keys())
    x = np.arange(len(labels))  # the label locations
    width = 0.1

    rects1 = ax.bar(x+placement*width, inorg.values(),
                    width, color=color, label=label)

    plt.title("Total N content by pool", size="x-large")
    plt.legend(fontsize="x-large")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, size="x-large")

    ax.bar_label(rects1, padding=3, size="small", fmt='%.1f', rotation=45)

    ax.set_ylabel('gN/m2', size="x-large")
    plt.rcParams['xtick.bottom'] = True

    # ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def C_content_bar_plot_split(data, placement, color, label, deep, fig, ax, ax1,spinup=False):
    if spinup:
        if deep == "tot":
            mass_dict = {}
            d = depth[0:len(data.levsoi)]   
            col = (data[pool_names]).isel(time=slice(-20, None)).mean(dim="time").to_dataframe().reset_index(drop=True)

            for var in col.columns:
                col_tot = col[var][0:len(depth)]*d
                sum_tot = round(col_tot.sum(), 2)
                mass_dict[var] = sum_tot

            fig.suptitle("Total C content by pool", size="x-large")
        elif deep == "30cm":

            fig.suptitle("C content down to 30cm by pool", size="x-large")
            mass_dict = {}
            depth_30 = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
            col = (data[pool_names]).isel(time=slice(-20, None)).mean(dim="time").to_dataframe().reset_index(drop=True)

            for var in col.columns:
                col30 = col[var][0:len(depth_30)]*depth_30
                sum_tot = round(col30.sum(), 2)
                mass_dict[var] = sum_tot
        else:
            print("deep either \"30cm\" or \"tot\"")
    else:
        if deep == "30cm":
            mass_dict = C30_pools(data)
            fig.suptitle("C content down to 30cm by pool", size="x-large")
        elif deep == "tot":
            mass_dict = Ctotal_pools(data)
            fig.suptitle("Total C content by pool", size="x-large")

        else:
            print("deep either \"30cm\" or \"tot\"")

    mic = {k: mass_dict[k] for k in ('SAPb', 'SAPf', 'EcM', "AM")}
    somlit = {k: mass_dict[k]
              for k in ("LITm", "LITs", "SOMp", "SOMc", "SOMa")}
    labels_mic = list(mic.keys())
    labels_som = list(somlit.keys())

    x_mic = np.arange(len(labels_mic))  # the label locations
    x_som = np.arange(len(labels_som))  # the label locations

    width = 0.06
    rects = ax.bar(x_som+placement*width, somlit.values(),
                   width, color=color, label=label)
    rects1 = ax1.bar(x_mic+placement*width, mic.values(),
                     width, color=color, label=label)
    ax.legend(fontsize="large")
    ax1.legend(fontsize="large")
    ax1.set_xticks(x_mic)
    ax1.set_xticklabels(labels_mic, size="x-large")
    ax.set_xticks(x_som)
    ax.set_xticklabels(labels_som, size="x-large")

    ax1.bar_label(rects1, padding=3, size="small", fmt='%.1f', rotation=45)
    ax.bar_label(rects, padding=3, size="small", fmt='%.1f', rotation=45)

    ax1.set_ylabel('gC/m2', size="x-large")
    ax.set_ylabel('gC/m2', size="x-large")

    plt.rcParams['xtick.bottom'] = True

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def N_content_bar_plot(data, placement, color, label, fig, ax):

    mass_dict = Ntotal_pools(data)

    labels = list(mass_dict.keys())
    x = np.arange(len(labels))  # the label locations
    width = 0.1

    rects1 = ax.bar(x+placement*width, mass_dict.values(),
                    width, color=color, label=label)

    plt.title("Total N content by pool", size="x-large")
    plt.legend(fontsize="x-large")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, size="x-large")

    ax.bar_label(rects1, padding=3, size="small", fmt='%.1f', rotation=45)

    ax.set_ylabel('gN/m2', size="x-large")
    plt.rcParams['xtick.bottom'] = True

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# returns the final and initial values of the xarray dataset


def final_init(data):  # (data = xarray)
    final = data.sel(time=data.time.values[-1])
    init = data.sel(time=data.time.values[0])
    return final, init
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


# To calulate carbon content down to 30cm and 50 cm to compare with observations (C30 C50) for decomp results
# NB: Uses mean of time period, so should be used mainly for comparing 1988-1992 period with observations!
def C30_C50_decomp(data):  # data = xarray
    sumC30 = 0
    sumC50 = 0

    # layer depths to multiply with concentrations
    depth_30 = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
    depth_50 = np.array([0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.02])  # [meters]

    col = data[pool_names]
    col = col.mean(dim="time")  # mean over time period
    col = col.to_dataframe().reset_index(drop=True)

    for var in col.columns:
        col30 = col[var][0:len(depth_30)]*depth_30
        col50 = col[var][0:len(depth_50)]*depth_50
        sum30 = round(col30.sum(), 3)
        sum50 = round(col50.sum(), 3)
        sumC30 += sum30
        sumC50 += sum50

    return sumC30, sumC50  # g/m2
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# C content down to 30 cm NB: Uses mean of time period, so should be used mainly for comparing 1988-1992 period with observations!


def C30_pools(data):
    pools_30 = {}
    # layer depths to multiply with concentrations
    depth_30 = np.array([0.02, 0.04, 0.06, 0.08, 0.10])

    col = data[pool_names]
    col = col.mean(dim="time")  # mean over timeperiod 1988-1992
    col = col.to_dataframe().reset_index(drop=True)

    for var in col.columns:
        col30 = col[var][0:len(depth_30)]*depth_30
        sum30 = round(col30.sum(), 2)
        pools_30[var] = sum30
    return pools_30
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Total C content in column. NB: Uses mean of time period, so should be used mainly for comparing 1988-1992 period with observations!


def Ctotal_pools(data):
    pools = {}
    d = depth[0:len(data.levsoi)]
    # layer depths to multiply with concentrations

    col = data[pool_names]
    col = col.mean(dim="time")  # mean over timeperiod 1988-1992
    col = col.to_dataframe().reset_index(drop=True)

    for var in col.columns:
        col_tot = col[var][0:len(depth)]*d
        sum_tot = round(col_tot.sum(), 2)
        pools[var] = sum_tot

    return pools
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Total N content in column. NB: Uses mean of time period, so should be used mainly for comparing 1988-1992 period with observations!


def Ntotal_pools(data):
    pools = {}
    d = depth[0:len(data.levsoi)]
    # layer depths to multiply with concentrations

    col = data[N_pool_names]
    col = col.mean(dim="time")  # mean over timeperiod 1988-1992
    col = col.to_dataframe().reset_index(drop=True)

    for var in col.columns:
        col_tot = col[var][0:len(depth)]*d
        sum_tot = round(col_tot.sum(), 2)
        pools[var] = sum_tot
    return pools
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# To calulate carbon content down to 30cm and 50 cm to compare with observations (C30 C50) for CLM results


def C30_C50_clm(data):
    sumC30 = 0
    sumC50 = 0

    # layer depths to multiply with concentrations
    depth_30 = np.array([0.02, 0.04, 0.06, 0.08, 0.10])
    depth_50 = np.array([0.02, 0.04, 0.06, 0.08, 0.12, 0.16, 0.02])  # [meters]

    col = data[['SOIL1C_vr', 'SOIL2C_vr', 'SOIL3C_vr',
                'LITR1C_vr', 'LITR2C_vr', 'LITR3C_vr']]
    col = col.mean(dim="time")
    col = col.to_dataframe().reset_index(drop=True)

    for var in col.columns:
        col30 = col[var][0:len(depth_30)]*depth_30
        col50 = col[var][0:len(depth_50)]*depth_50
        sum30 = round(col30.sum(), 3)
        sum50 = round(col50.sum(), 3)
        sumC30 += sum30
        sumC50 += sum50

    return sumC30, sumC50


# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def for_barplots(pd_data, pool_list, depth):
    final_pools = pd.DataFrame(columns=pool_list)
    total_mass = 0
    for p in pool_list:
        final_pools[p] = np.array(
            pd_data[p])*depth[0:len(pd_data.index)]  # From g/m3 to g/m2
        total_mass += math.fsum(final_pools[p])
    return final_pools, total_mass
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def bar_plot(df, placement, legend_on, ax):
    columns = df.columns
    rows = df.index
    # Get some shades for the colors
    colors = plt.cm.copper(np.linspace(0, 1, len(rows)))
    n_rows = len(df)
    index = np.arange(len(columns))
    bar_width = 0.2
    # Initialize the vertical-offset for the stacked bar chart.
    y_offset = np.zeros(len(columns))
    y_offset = df.sum().values

    # Plot bars and create text labels for the table
    for row in range(n_rows):
        y_offset = y_offset - df.iloc[row].values
        plt.bar(index+placement*bar_width,
                df.iloc[row].values, bar_width, bottom=y_offset, label="%i" % row, color=colors[row])
    if legend_on == True:
        plt.legend()
    ax.set_xticks(index)
    ax.set_xticklabels(columns, size="large")

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def compare_two_variables(data1, label1, data2, label2, variable):
    for i in range(min(no_of_layers(data1), no_of_layers(data2))):
        plt.figure(figsize=[20, 7])
        data1[variable][i].plot(label=label1)
        data2[variable][i].plot(label=label2)
        plt.legend()
        plt.title("Layer %s" % (i+1))

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def calculate_fungi_decomp_fluxes(data):
    # Decomposed by fungi:
    LITm_to_SAPf = data.C_LITmSAPf*data.CUEf
    LITs_to_SAPf = data.C_LITsSAPf*data.CUEf
    SOMa_to_SAPf = data.C_SOMaSAPf*data.CUEf

    return LITm_to_SAPf, LITs_to_SAPf, SOMa_to_SAPf


def calculate_bacteria_decomp_fluxes(data):
    # Decomposed by bacteria:
    LITm_to_SAPb = data.C_LITmSAPb*data.CUEb
    LITs_to_SAPb = data.C_LITsSAPb*data.CUEb
    SOMa_to_SAPb = data.C_SOMaSAPb*data.CUEb

    return LITm_to_SAPb, LITs_to_SAPb, SOMa_to_SAPb
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Copied from: https://wil.yegelwel.com/cluster-correlation-matrix/


def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 

    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 

    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold,
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)

    if not inplace:
        corr_array = corr_array.copy()

    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]
    return corr_array[idx, :][:, idx]

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


def mean_turnover_time(data, pool_name, list_of_loss_terms, return_bool=False):
    # Use to convert from gC/m3 to gC/m2
    depths = xr.DataArray(depth[:no_of_layers(data)], dims="levsoi")
    # mean value over the whole time period for each layer
    mean_layer_mass = depths*data[pool_name].mean(axis=1)
    t = []
    losses = []
    placement = [['Left', 'TopRight'],['Left','CenterRight'],['Left','BottomRight'],['Left','VeryBottomRight']]
    fig = plt.figure(constrained_layout=True,figsize=[15, 8])
    axes = fig.subplot_mosaic(placement,gridspec_kw={"height_ratios":[1,1,1,1]})


    # fig1, axes = plt.subplots(nrows=2, ncols=len(
    #     # list_of_loss_terms), figsize=[20, 12])

    fig.suptitle(data.site_name+data.run_name)

    for l in range(len(list_of_loss_terms)):
        p = placement[l][1]
        # tau_l = m/l, converted to years from hrs
        tau = (mean_layer_mass /
               (depths*data[list_of_loss_terms[l]]).mean(axis=1))/(24*365)
        t.append(tau)
        losses.append((depths*data[list_of_loss_terms[l]]).mean(axis=1))
        axes[p].plot(tau, tau.levsoi, '-*')
        axes[p].invert_yaxis()
        axes[p].set_title("Turnover time to loss term: %s for pool %s" % (
            list_of_loss_terms[l], pool_name))
        axes[p].set_ylabel("layer number")
        axes[p].set_xlabel("Mean turnover time in years")

        axes[p].grid(True)

    tau_tot = ((mean_layer_mass/sum(losses))/(24*365)).round(2)
    #plt.figure(figsize=[15, 8])
    axes["Left"].plot(tau_tot, tau_tot.levsoi, '-*')
    axes["Left"].grid()
    axes["Left"].invert_yaxis()
    axes["Left"].set_title("Mean turnover time for each layer of %s" % pool_name)
    axes["Left"].set_ylabel("layer number")
    axes["Left"].set_xlabel("Mean turnover time in years")

    #tau_tot: mass/(L1+L2..Ln) for each layer
    #t: list: [mass/L1, mass/L2, ... , mass/Ln]
    #mean_layer_mass: mean content in each layer [g/m2]


    if return_bool:
        return (tau_tot, t, mean_layer_mass)
