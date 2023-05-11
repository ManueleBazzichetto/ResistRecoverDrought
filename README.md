# ResistRecoverDrought
This repo includes (R) scripts and data to replicate the results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought"

Data used in the analyses can be found in the 'Data' folder. Note that SPEI data (v 2.6), which we used to quantify the intensity of drought, are openly available at: https://digital.csic.es/handle/10261/202305

Brief description of data found in the 'Data' folder:

- biomass.csv: csv file with (yearly) values of vegetation biomass and date of biomass cut (the latter reported in a column named 'day_of_year');
- plots.csv: csv file with (yearly) soil moisture data for each vegetation plot;
- rich+lui: .RData file with (yearly) values of species richness (i.e., number of plant species recorded each year, in each vegetation plot) and of the land-use index;
- fundata: .RData file with (yearly) values of community weighted means (computed for above-ground plant traits), and of different indices of functional diversity (including measures 'coupling' functional and phylogenetic diversity). 

Brief description of scripts to replicate the results presented in the manuscript: 

Please, run the scripts following this order: Functions -> Data_prep -> SPEI_data -> Analysis -> Analysis2 -> Diagnostics -> ModelTables

- Functions: this script reports several functions, including those for: (i) computing the annual log ratio; (ii) performing a principal component analysis on the CWMs of the plant traits; and (iii) fit linear regression for year-specific analysis (see script 'Analysis');
- Data_prep: this script provides the code for generating the data needed to run the analysis in scripts 'Analysis' and 'Analysis2';
- Analysis: this script reports year-specific analysis of the log response ratios, which we modelled as a function of biodiversity and land-use intensity;
- Analysis2: this script reports the analysis of the effect of the interaction between drought intensity and biodiversity on biomass resistance and recovery;
- Diagnostics: diagnostics of the models fitted in 'Analysis';
- ModelTables: this script allows reproducing the tables of model summaries.   

**All data made available in this repo are released under CC-BY license**
