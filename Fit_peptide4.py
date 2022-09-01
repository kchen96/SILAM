#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analysis of degradation rates across oxygen.

Project: SILAM

Step #4

Input data: Filtered normalized peptide abundance.

Output data: Mixed linear model effects analysis.

Created on 2022.06.09

Author: Kirsten Chen

Last changed: 2022.07.22
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats import *

import warnings
warnings.filterwarnings('ignore')


## Load data
df = pd.read_csv('peptide_level_filtered.csv')
fitRes = pd.read_csv('peptide_fit.csv')


########## Part IV -- Mixed linear effects 

## Fold changes
fitRes_fc = pd.DataFrame(columns = ['Protein.Group', 'Genes',
                                    'Oxy8vs21_Log2FC', 'Oxy60vs21_Log2FC'])

for i in np.arange(fitRes.shape[0]):
    p = fitRes.Genes.values[i]
    pg = fitRes['Protein.Group'].values[i]
    kd_21 = fitRes[np.logical_and(fitRes.Genes == p,fitRes.Oxygen==21)]['Kd'].values
    kd_8 = fitRes[np.logical_and(fitRes.Genes == p,fitRes.Oxygen==8)]['Kd'].values
    kd_60 = fitRes[np.logical_and(fitRes.Genes == p,fitRes.Oxygen==60)]['Kd'].values
    
    if len(kd_21)>0:
        kd_21 = kd_21[0]
        
        if len(kd_8) >0: 
            kd_8 = kd_8[0]
            lfc_8 = np.log2(kd_8/kd_21)
        else:
            lfc_8 = np.nan
    
        if len(kd_60) >0: 
            kd_60 = kd_60[0]
            lfc_60 = np.log2(kd_60/kd_21)
            
        else:
            lfc_60 = np.nan
            
    else: 
        lfc_8 = np.nan
        lfc_60 = np.nan
        
    fitRes_fc.loc[len(fitRes_fc.index)] = [pg, p, lfc_8, lfc_60]

fitRes_fc = fitRes_fc.set_index('Protein.Group')


## Mixed linear effects
def MixedML_peptides(d):
    subdf = d[d.value != 0]
      

    if subdf.shape[0] > 1 : 
        ## Extract y (scaled values), take the natural log
        y = np.array(subdf['scaled_value'], dtype='float')
        subdf['ln_value'] = np.log(y) 
        
        subdf1 = subdf[subdf.Oxygen != 60]
        subdf2 = subdf[subdf.Oxygen != 8]
        
        ## Extract x (days)
        x1 = np.array(subdf1['Day'], dtype='float').reshape(-1, 1)
        x2 = np.array(subdf2['Day'], dtype='float').reshape(-1, 1)

        if subdf1.shape[0] > 0 and len(np.unique(x1))>1:
            ## fit
            try:
                md1 = smf.mixedlm("ln_value ~ 0+ Day + Oxygen + Day*Oxygen", subdf1, groups=subdf1['Modified.Sequence'])
                mdf1 = md1.fit(method=["lbfgs"])
                mdf1_p = mdf1.pvalues['Day:Oxygen']
            except:
                mdf1_p = np.nan
        else:
            mdf1_p = np.nan
            
        if subdf2.shape[0] > 0 and len(np.unique(x2))>1:
            ## fit
            try:
                md2 = smf.mixedlm("ln_value ~ 0+ Day + Oxygen + Day*Oxygen", subdf2, groups=subdf2['Modified.Sequence'])
                mdf2 = md2.fit(method=["lbfgs"])
                mdf2_p = mdf2.pvalues['Day:Oxygen']
            except:
                mdf2_p = np.nan
        else:
            mdf2_p = np.nan
        
    else:
        mdf1_p = np.nan
        mdf2_p = np.nan
        
    return pd.DataFrame({'Oxy8vs21_pval': mdf1_p,
                         'Oxy60vs21_pval': mdf2_p}, index = [0])




out = df.groupby(['Protein.Group','Genes']
                 ).apply(lambda x: MixedML_peptides(x))
out = out.reset_index()                      

out1 = out[['Protein.Group','Oxy8vs21_pval']].dropna().set_index('Protein.Group')
out1['Oxy8vs21_FDR'] = multitest.multipletests(pvals = out1.Oxy8vs21_pval.values,
                                          alpha = 0.05, method = 'fdr_bh')[1]

out2 = out[['Protein.Group','Oxy60vs21_pval']].dropna().set_index('Protein.Group')
out2['Oxy60vs21_FDR'] = multitest.multipletests(pvals = out2.Oxy60vs21_pval.values,
                                          alpha = 0.05, method = 'fdr_bh')[1]

result = pd.merge(fitRes_fc['Oxy8vs21_Log2FC'].dropna(), out1, how='left',left_index=True, right_index=True).drop_duplicates()
result2 = pd.merge(fitRes_fc['Oxy60vs21_Log2FC'].dropna(), out2, how='left',left_index=True, right_index=True).drop_duplicates()
result = pd.merge(result, result2,how='outer',left_index=True, right_index=True).drop_duplicates()
result = pd.merge(fitRes_fc['Genes'], result, how='right', left_index=True, right_index=True).drop_duplicates()

result.drop(result[np.logical_and(np.isnan(result['Oxy8vs21_Log2FC']), 
                                  np.isnan(result['Oxy60vs21_Log2FC']))].index, inplace = True)

result.to_csv('peptide_fit_MixedLM.csv')
