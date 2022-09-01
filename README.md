# SILAM

Fit peptides for degradation constant and half-life.

Project: SILAM

Pipeline:
	1	Fit_peptide1.py:  Scale peptides from 0 to 1, fit the decay of each peptide R^2
	2	Fit_peptide2.py:  Based on the R^2, filter for “good” peptides (I set the cutoff as 0.6) and scale each peptide from 0 to 1
	3	Fit_peptide3.py:  Use mixed linear model to calculate the degradation rate constant from the peptide-level data
	4	Fit_peptide4.py:  Use mixed linear model to analyze the effects of oxygen on the degradation rate constant from the peptide-level data


Input data: Peptide mass spec readout, normalized based on long-lived proteins.

Output data:
    1. Results of single peptide fitting.
    2. Filtered peptides based on the R2.
    3. Results of linear fitting based on the "good" peptides.
    4. Mixed linear model analysis on the oxygen effects. 

Updates in this version:
    1. Use the mixed linear model (statsmodels mixedlm)


Created on 2022.04.24

Author: Kirsten Chen

Last changed: 2022.06.13
