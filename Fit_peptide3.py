#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fit peptides for degradation constant and half-life.

Project: SILAM

Step #3

Input data: Filtered normalized peptide abundance.

Output data: Results of fitting across all replicates.

Created on 2022.04.24

Author: Kirsten Chen

Last changed: 2022.09.14
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

import warnings
warnings.filterwarnings('ignore')


## Load data
df = pd.read_csv('peptide_level_filtered.csv')



########## Part III -- Compute degradation rate constant and half lives based on protein-level data
def fit_single_peptide(d):
    subdf = d[d.value != 0]
    
    ## Extract x (days)
    x = np.array(subdf['Day'], dtype='float').reshape(-1, 1)

    if subdf.shape[0] > 1 and len(np.unique(x)) > 1: 
        ## Extract y (scaled values), take the natural log
        y = np.array(subdf['scaled_value'], dtype='float')
        subdf['ln_value'] = np.log(y) 
        
        try:
            ## fit
            md = smf.mixedlm("ln_value ~ 0+ Day", subdf, groups=subdf['Modified.Sequence'])
            mdf = md.fit(method=["lbfgs"])
            
        except:
            return pd.DataFrame({'Kd': np.nan, 'Kd_pvalue':np.nan, 
                                 'Half_life': np.nan}, index=[0]) 

    
    else:
        return pd.DataFrame({'Kd': np.nan, 'Kd_pvalue':np.nan, 
                             'Half_life': np.nan}, index=[0]) 


    ## Rate constant
    kd = -mdf.params['Day']
    kd_pval = mdf.pvalues['Day']
    
    
    ## Half life based on the degradation rate
    hl = (np.log(2)/(kd)) if kd > 0 else np.nan
    
    
    return pd.DataFrame({'Kd':kd, 'Kd_pvalue': kd_pval,
                         'Half_life': hl}, index=[0])



out = df.groupby(['Protein.Group','Genes','Oxygen']).apply(lambda x: fit_single_peptide(x))
out = out.reset_index().drop(columns = ['level_3'])                          

out.dropna(inplace = True)
out.to_csv('peptide_fit.csv')
