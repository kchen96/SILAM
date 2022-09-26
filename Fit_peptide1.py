#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fit peptides for degradation constant and half-life.

Project: SILAM

Step #1

Input data: Peptide mass spec readout, normalized based on long-lived proteins.

Output data: Results of single peptide fitting.

Created on 2022.04.24

Author: Kirsten Chen

Last changed: 2022.09.24
"""

import pandas as pd
import numpy as np
import dask.dataframe as dd
import statsmodels.api as sm

import warnings
warnings.filterwarnings('ignore')



########## Part I -- Fit single peptide and find R2
# Function for peptide fitting
def fit_single_peptide(d, protein_level = True):
    subdf = d[d.value != 0]
    
    ## Extract x (days)
    x = np.array(subdf['Day'], dtype='float').reshape(-1, 1)
      
    if subdf.shape[0] > 0 and len(np.unique(x)) > 1:

        ## Extract mass spec peak intensities
        y = np.array(subdf['value'], dtype='float')
        
        ## Scale the y values
        ymax = np.median(subdf[subdf.Day == 0].value.values)
        y = y  / ymax 
        ## Natural log 
        y = np.log(y) 

    else:
        return pd.Series({'R2': -1,'Kd': np.nan, 'Kd_SE':np.nan,
                          'Half_life': np.nan,'ymax': np.nan}) 


    ## fit
    result = sm.OLS(y, x).fit()
    ## Calculate the slope (degradation rate), standard error of the slope
    kd = -result.params[0]
    bse = result.bse[0]
    ## Half life based on the degradation rate
    hl = (np.log(2)/(kd)) if kd > 0 else np.nan
    
    
    return pd.Series({'R2': result.rsquared,'Kd': kd, 'Kd_SE': bse,
                      'Half_life': hl, 'ymax': ymax})

# Import normalized peptide datatable
df = pd.read_csv('peptide_level.csv',index_col = 0)
df = df[(df['normalized']==True) & (df['quant']=='MS2.Area')]

df.fillna(0,inplace=True)
ddf = dd.from_pandas(df, npartitions=8)
ddf = ddf[(ddf['Protein.Group'].str.contains(';')== False)]
ddf.compute().to_csv('peptide_level_filtered.csv')

out = ddf.groupby(['Protein.Group','Genes','Modified.Sequence','Oxygen']
                  ).apply(lambda x: fit_single_peptide(x, protein_level = False),
                          meta = {'R2': 'float64', 'Kd': 'float64',
                                  'Kd_SE': 'float64',
                                  'Half_life': 'float64',
                                  'ymax': 'float64'})

reduced_df = out.compute()
reduced_df.to_csv('peptide_fit_rsquareddata.csv')
