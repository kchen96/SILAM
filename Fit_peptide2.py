#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fit peptides for degradation constant and half-life.

Project: SILAM

Step #2

Goals: (1) Filter for peptides with R^2 passed the threshold;
(2) Scale the value for each peptide from 0 to 1 based on the max value.

Input data: (1) Peptide mass spec readout, normalized based on long-lived proteins; 
(2) Fitting results from Step #1, which also contains the max value of each peptide that will be used for scaling.

Output data: Filtered normalized peptide abundance. Calculation of scaled values for each peptide.

Created on 2022.04.24

Author: Kirsten Chen

Last changed: 2022.06.09
"""

import pandas as pd
import numpy as np
import dask.dataframe as dd

import warnings
warnings.filterwarnings('ignore')

## Pre-processing
df = pd.read_csv('peptide_level_filtered.csv',index_col = 0)
df = df[(df['normalized']==True) & (df['quant']=='MS2.Area')]

df.fillna(0,inplace=True)
ddf = dd.from_pandas(df, npartitions=32)


fit1 = pd.read_csv('peptide_fit_rsquareddata.csv')
fit1 = dd.from_pandas(fit1[['Modified.Sequence','Oxygen','R2','ymax']], npartitions=32) ## only these 4 columns are useful



########## Part II -- Filter 'good' peptides based on the R2
R2_cutoff = 0.6  # Set an arbitrary cutoff on R2


## Split datasets into 3 groups
ddf1 = ddf[ddf.Oxygen==8]
ddf2 = ddf[ddf.Oxygen==21]
ddf3 = ddf[ddf.Oxygen==60]

ff1 = fit1[fit1.Oxygen==8].drop('Oxygen',axis=1)
ff2 = fit1[fit1.Oxygen==21].drop('Oxygen',axis=1)
ff3 = fit1[fit1.Oxygen==60].drop('Oxygen',axis=1)

## merge 
ddf1 = ddf1.merge(ff1, how = 'left', on =['Modified.Sequence'])
ddf2 = ddf2.merge(ff2, how = 'left', on =['Modified.Sequence'])
ddf3 = ddf3.merge(ff3, how = 'left', on =['Modified.Sequence'])

## Filter
ddf1 = ddf1[ddf1['R2'] > R2_cutoff]
ddf2 = ddf2[ddf2['R2'] > R2_cutoff]
ddf3 = ddf3[ddf3['R2'] > R2_cutoff]

## Calculate scaled y
ddf1['scaled_value'] = ddf1['value']/ddf1['ymax']
ddf2['scaled_value'] = ddf2['value']/ddf2['ymax']
ddf3['scaled_value'] = ddf3['value']/ddf3['ymax']

ddf = ddf1.append(ddf2)
ddf = ddf.append(ddf3)

filtered_df = ddf.compute()

filtered_df.to_csv('peptide_level_filtered.csv')