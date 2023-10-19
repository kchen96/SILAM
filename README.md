# SILAM

Fit peptides for degradation constant and half-life.

Project: SILAM

Flow chart:

<img width="827" alt="Screen Shot 2023-10-18 at 5 41 51 PM" src="https://github.com/kchen96/SILAM/assets/53800187/33882515-335a-4ac1-9eaa-972e436e2e13">

Analysis pipeline:

1. **Mass spectrometry data processing**: The DIA-PASEF data was searched with DIA-NN v.1.7.1 using a library-centric approach. Identified spectrums with MS1 precursors within 10 ppm and MS2 precursor within 15 ppm were selected and a second library was generated (double-pass mode). Quantification was set to robust (high-accuracy) while signal was re-normalized as function of the spray stability (RT-dependent). Protein inference was disabled, and library generation was set to smart profiling. The transition level data was filtered at 1% library Q-value (PG) and transitions were summed into precursor MS2 abundances and precursors were averaged to a single peptide-abundance. Only the abundance of all the unlabeled peptides (no 15N) is used for downstream analysis.

2. **Normalization**: Normalize data internally based on the median abundance of all peptides of the 10 long-lived proteins (half-lives > 2 months in mouse heart), under the assumption that fixed abundance between different conditions and labelling time.

3. **Scaling data**: The normalized abundance was divided by the median abundance at time 0 of the labeling to calculate the fraction of unlabeled peptides at each time point. This way, the median value for time 0 will be 1, and the scaled value of peptides should be between 0 and 1.

4. **Initial fitting**: Next, a preliminary fitting was performed to filter peptides with good linear correlations using the OLS function from the python (v3.11) statsmodels package (v0.13.5). Each peptide for a given condition was fit into the first-order kinetic model using the ordinary linear model (x = time. y = ln(scaled_value)). 
 
5. **Filter peptides**: The goodness of fitting was assessed by r2. (1) r2 > 0.65. (2) Peptides are detected in more than 5 samples in a given condition. (3) Peptides are detected in at least 1 sample at day 0
	
6. **Calculate the Kd for each protein at each oxygen tension**: To calculate the Kd for each protein, a linear mixed effects model was applied using the mixedlm function from the python (v3.11) statsmodels package (v0.13.5). For each oxygen condition and each tissue, time was set as the fixed variable, and to account for the variability of different peptides for a given protein, the peptides were set as the random variable. Kd values with p-values less than 0.20 were reported. Protein half-lives (t1/2) were calculated based on 
t_(1/2)=(ln⁡(2))/K_d 

7. **Compare Kd between oxygen tensions**: To compare Kd between oxygen tensions, a linear mixed effects model using mixedlm with oxygen, time, and oxygen:time as fixed effects, and peptides as random variables was applied. The p values of the interaction between oxygen and time (oxygen:time) indicated the differences in the slopes, which assessed the differences of the degradation rates. Adjusted p values were calculated using the Benjamini-Hochberg correction.



Scripts: 

	0. Normalization:  Normalize the MS area based on the long-lived proteins
	1. Fit_peptide1.py:  Scale peptides from 0 to 1, fit the decay of each peptide R^2
	2. Fit_peptide2.py:  Based on the R^2, filter for “good” peptides (I set the cutoff as 0.6) and scale each peptide from 0 to 1
	3. Fit_peptide3.py:  Use mixed linear model to calculate the degradation rate constant from the peptide-level data
	4. Fit_peptide4.py:  Use mixed linear model to analyze the effects of oxygen on the degradation rate constant from the peptide-level data


**Input data:** Peptide mass spec readout (DIA-NN output).


Output data:

    0. Peptide mass spec readout (DIA-NN output), normalized based on long-lived proteins.
    1. Results of single peptide fitting.
    2. Filtered peptides based on the R2.
    3. Results of linear fitting based on the "good" peptides.
    4. Mixed linear model analysis on the oxygen effects. 
    



Created on 2022.04.24

Author: Kirsten Chen

Last changed: 2023.10.18
