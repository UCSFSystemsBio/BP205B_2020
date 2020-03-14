#### Suvival Analysis results

Survival analysis carried out in the METABRIC cohort - breast cancer patients with normalized gene expression (z-scored) + metadata

Signatures obtained from HCC_MK and MDA data-sets are scored for each patient by summing the normalized expression value for genes in a given signature. 
Patients are then sorted into two groups - high score or low score - by setting a threshold. The threshold is identified as the value which maximizes the difference between two groups.
This can sometimes result in unbalanced groups (i.e 5 in group, 1995 in the other). 

For genes that are 'downregulated' --> "DOWN", or 'upregulated' --> 'UP, expression is summed without multiplying by -1. For aggregate gene signatures with both UP and DOWN components (lacking a DOWN or UP label), score is calculated as score = UP + -1*DOWN. 

Kaplan-Meir survival curves (km_plots) as well as hazard ratio plots or forrest plots (hazard_plots) are provided for each signature. Factors are leveled so that the 'normal' condition is the reference, where it makes sense. 

The CSV file survival_results.csv consists of the following: (columns renamed slightly)
- KM_pval : The Kaplan-Meir p value - if signficant indicates signficant difference in # of people surviving in each group
- KM_HZ : the hazard ratio for being in one group vs othe other -- the larger the more signficant the effect, in either direction (> 1 or < 1)
- KM_HZ_0.95_lower : the lower end of the 95% confidence interval for the KM hazard ratio
- KM_HZ_0.95_upper : the upper end of the 95% confidence interval for the KM hazard ratio
- Cox_HZ_... : the hazard ratio for being in one group vs another generated from the Cox model, 95% conf.int is not included at the moment
- Cox_pval_... : the p value for being in one group vs another generated from the Cox model
---- Specifically, Breast cancer subtype status (Subtype), chemotherapy treatment (Chemotherapy), and tumor cellularity (Cellularity) are provided