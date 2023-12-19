# ResistRecoverDrought 
[![DOI](https://zenodo.org/badge/620832781.svg)](https://zenodo.org/badge/latestdoi/620832781)

**The repo is going to be updated soon by a third release!**

This repo includes (R) scripts and data to replicate the results presented in the manuscript **"NEW TITLE"**.

Data used in the analyses can be found in the 'Data' folder. Note that SPEI data (v 2.6), which we used to quantify the intensity of drought, are openly available at: https://digital.csic.es/handle/10261/202305

Brief description of data found in the 'Data' folder:

- biomass.csv: csv file with (yearly) values of vegetation biomass and date of biomass cut (the latter reported in a column named 'day_of_year');
- plots.csv: csv file with (yearly) soil moisture data for each vegetation plot;
- rich+lui: .RData file with (yearly) values of species richness (i.e., number of plant species recorded each year, in each vegetation plot) and of the land-use index;
- fundata: .RData file with (yearly) values of community weighted means (computed for above-ground plant traits), and of different indices of functional diversity (including measures 'coupling' functional and phylogenetic diversity). 

Brief description of scripts to replicate the results presented in the manuscript: 

Please, run the scripts following this order: Functions -> Data_prep -> SPEI_data -> Time_series_log_ratios -> Recovery_resistance -> ModelTables

- Functions: this script includes the R functions needed for (i) computing the annual log ratio (in Data_prep); and (ii) performing a principal component analysis on the CWMs of the plant traits (also happening in Data_prep);
- Data_prep: essentially, this script is for generating and processing all data needed to run the analysis in scripts 'Time_series_log_ratios' and 'Recovery_resistance';
- Time_series_log_ratios: this script reports the time-series analysis of the log response ratios, which are modelled as a function of biodiversity and land-use intensity;
- Recovery_resistance: this script reports the analysis on the effect of the interaction between drought intensity and biodiversity on biomass resistance and recovery;
- ModelTables: this script allows reproducing the tables of model summaries (presented in the supplementary material of the manuscript).   

**All data made available in this repo are released under CC-BY license 4.0**
