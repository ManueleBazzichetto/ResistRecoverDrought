##Code to replicate results presented in the manuscript "Biodiversity promotes resistance but dominant species shape recovery of grasslands under extreme drought" by Bazzichetto et al.

##Model summary tables for analyses in the Recovery_resistance script

library(modelsummary)
library(performance)

##----------------------------------------------------Get tables of models' summary

##----------------------------------------------------------------Recovery

##Functional components

##SPEI-3: Recovery_SPEI3
##SPEI-12: Recovery_SPEI12
##SPEI-24: Recovery_SPEI24

##Taxonomic component

##SPEI-3: Recovery_SPEI3.tx
##SPEI-12: Recovery_SPEI12.tx
##SPEI-24: Recovery_SPEI24.tx

Coef_names_models <- c("(Intercept)" = "Intercept",
                     "SPEI_final_catModerateWet" = "ModerateWet",
                     "SPEI_final_catExtremeWet" = "ExtremeWet",
                     "SPEI_final_catExtremeDry" = "ExtremeDry",
                     "ExploHAI" = "Central",
                     "ExploSCH" = "North-East",
                     "Across_PC1_CWM" = "Slow-fast continuum",
                     "Yr_Comb_rao_FP" = "Functional diversity",
                     "Yr_species_rich" = "Species richness",
                     "Yr_lui" = "Land-use intensity",
                     "day_of_year" = "Day of year",
                     "SoilMoist10" = "Soil moisture",
                     "SPEI_final_catModerateWet:Across_PC1_CWM" = "ModerateWet:Slow-fast c.",
                     "SPEI_final_catExtremeWet:Across_PC1_CWM" = "ExtremeWet:Slow-fast c.",
                     "SPEI_final_catModerateWet:Yr_Comb_rao_FP" = "ModerateWet:Funct. div.",
                     "SPEI_final_catExtremeWet:Yr_Comb_rao_FP" = "ExtremeWet:Funct. div.",
                     "SPEI_final_catExtremeDry:Across_PC1_CWM" = "ExtremeDry:Slow-fast c.",
                     "SPEI_final_catExtremeDry:Yr_Comb_rao_FP" = "ExtremeDry:Funct. div.",
                     "SPEI_final_catModerateWet:Yr_species_rich" = "ModerateWet:Sp. rich.",
                     "SPEI_final_catExtremeWet:Yr_species_rich" = "ExtremeWet:Sp. rich.", 
                     "SPEI_final_catExtremeDry:Yr_species_rich" = "ExtremeDry:Sp. rich.")

#models for recovery functional components
RecoFun.list <- list(SPEI3 = Recovery_SPEI3, SPEI12 = Recovery_SPEI12, SPEI24 = Recovery_SPEI24)

#R2 and adj-R2 computed by modelsummary are same of those computed by performance::r2
lapply(RecoFun.list, r2)

#docx
modelsummary(models = RecoFun.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Reco_Fun.docx")

#models for recovery taxonomic component
RecoTx.list <- list(SPEI3 = Recovery_SPEI3.tx, SPEI12 = Recovery_SPEI12.tx,
                    SPEI24 = Recovery_SPEI24.tx)

#R2 and adj-R2 computed by modelsummary are same of those computed by performance::r2
lapply(RecoTx.list, r2)

#docx
modelsummary(models = RecoTx.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Reco_Tx.docx")

#Anova tables for functional components

#Incremental F-test
car::Anova(Recovery_SPEI3)
car::Anova(Recovery_SPEI12)
car::Anova(Recovery_SPEI24)

#Anova tables for taxonomic component

#Incremental F-test
car::Anova(Recovery_SPEI3.tx)
car::Anova(Recovery_SPEI12.tx)
car::Anova(Recovery_SPEI24.tx)


##----------------------------------------------------------------Resistance

##Functional components

##SPEI-3: Resistance_SPEI3
##SPEI-12: Resistance_SPEI12
##SPEI-24: Resistance_SPEI24

##Taxonomic component

##SPEI-3: Resistance_SPEI3.tx
##SPEI-12: Resistance_SPEI12.tx
##SPEI-24: Resistance_SPEI24.tx

#models for recovery functional components
ResiFun.list <- list(SPEI3 = Resistance_SPEI3, SPEI12 = Resistance_SPEI12, SPEI24 = Resistance_SPEI24)

#R2 and adj-R2 computed by modelsummary are same of those computed by performance::r2
lapply(ResiFun.list, r2)

#docx
modelsummary(models = ResiFun.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Resi_Fun.docx")

#models for recovery taxonomic component
ResiTx.list <- list(SPEI3 = Resistance_SPEI3.tx, SPEI12 = Resistance_SPEI12.tx,
                    SPEI24 = Resistance_SPEI24.tx)

#R2 and adj-R2 computed by modelsummary are same of those computed by performance::r2
lapply(ResiTx.list, r2)

#docx
modelsummary(models = ResiTx.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Resi_Tx.docx")

#Anova tables for functional components

#Incremental F-test
Anova(Resistance_SPEI3)
Anova(Resistance_SPEI12)
Anova(Resistance_SPEI24)

#Anova tables for taxonomic component

#Incremental F-test
Anova(Resistance_SPEI3.tx)
Anova(Resistance_SPEI12.tx)
Anova(Resistance_SPEI24.tx)

