# TPL-H1_Mechanism
A collection of scripts used to perform computational analysis related to the TPL-H1 mechanism

Details of code package:

For APEX2 Data analysis

GO_APEX2_graphing.R
Plotting GO term enrichment on proximity labeling results. 

For R-SGA analysis pipelines
SGA_Analysis.R:
1. Read colony size files and append rows/rbind.
2. Add conditions and replicate vectors.
3. Cbind colony sizes with replicate and condition.
4. Read 1536 FLEX map and cbind. 
5. Name file and write to computer.

SGA_analysis_TS.R:
This code is the variation (from above) for temperature sensitive arrays.
Array_colony_size_measurement.R is used to extract colony size for all array images in a given folder.

Plot_DA/TS_Zscores.R:
This set of code was used to plot heatmaps of the combined Z-scores from all R-SGA experiments (HT1), the APEX2 proximity labelling scores (HT2), the validated ratio by cytometry for the SPARC-H1-H5 (HT4), and SPARC-TPLN188 (HT5) liquid growth experiments. This code was used to create the . The DA refers to the Deletion Array R-SGA experiments, and the TS refers to the Temperature Sensitive Array R-SGA experiments. These are kept seperate due to the use of an alternative array ID in the TS array, which has multiple alleles of the same gene, rendering use of the normal yeast gene locus impossible. 


RSGA_GO_Analysis.R
Plotting GO term enrichment on R-SGA upregulated hits results. 
