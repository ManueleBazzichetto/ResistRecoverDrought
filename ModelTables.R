##Code to replicate results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought" by Bazzichetto et al.
##Model summary tables

library(modelsummary)
library(performance)

##----------------------------------------------------Get tables of models' summary

##----------------------------------------------------------------Recovery

##Functional components

##SPEI-3: Recovery_SPEI3.nlme
##SPEI-12: Recovery_SPEI12.nlme
##SPEI-24: Recovery_SPEI24 (lm)

##Taxonomic component

##SPEI-3: Recovery_SPEI3.tx.nlme
##SPEI-12: Recovery_SPEI12.tx.nlme
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
                     "Yr_Myco_intensity" = "Mycorrhizae int.",
                     "day_of_year" = "Day of harvest",
                     "SPEI_final_catModerateWet:Across_PC1_CWM" = "ModerateWet:Slow-fast c.",
                     "SPEI_final_catExtremeWet:Across_PC1_CWM" = "ExtremeWet:Slow-fast c.",
                     "SPEI_final_catModerateWet:Yr_Comb_rao_FP" = "ModerateWet:Funct. div.",
                     "SPEI_final_catExtremeWet:Yr_Comb_rao_FP" = "ExtremeWet:Funct. div.",
                     "SPEI_final_catExtremeDry:Across_PC1_CWM" = "ExtremeDry:Slow-fast c.",
                     "SPEI_final_catExtremeDry:Yr_Comb_rao_FP" = "ExtremeDry:Funct. div.",
                     "SPEI_final_catModerateWet:Yr_species_rich" = "ModerateWet:Sp. rich.",
                     "SPEI_final_catExtremeWet:Yr_species_rich" = "ExtremeWet:Sp. rich.", 
                     "SPEI_final_catExtremeDry:Yr_species_rich" = "ExtremeDry:Sp. rich.",
                     "SD (Intercept Useful_EP_PlotID)" = "SD PlotID",
                     "SD (Observations)" = "SD Residual")

#models for recovery functional components
RecoFun.list <- list(SPEI3 = Recovery_SPEI3.nlme, SPEI12 = Recovery_SPEI12.nlme, SPEI24 = Recovery_SPEI24)

#marg and cond R2 computed by modelsummary are same of those computed by performance::r2
lapply(RecoFun.list, r2)

#docx
modelsummary(models = RecoFun.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Reco_Fun.docx")

#models for recovery taxonomic component
RecoTx.list <- list(SPEI3 = Recovery_SPEI3.tx.nlme, SPEI12 = Recovery_SPEI12.tx.nlme,
                    SPEI24 = Recovery_SPEI24.tx)

#marg and cond R2 computed by modelsummary are same of those computed by performance::r2
lapply(RecoTx.list, r2)

#docx
modelsummary(models = RecoTx.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Reco_Tx.docx")

#Deviance/Anova tables for functional components

#Kenward-Roger
car::Anova(Recovery_SPEI3.lmer, test = "F")
car::Anova(Recovery_SPEI12.lmer, test = "F")
#Incremental F-test
car::Anova(Recovery_SPEI24)

#Deviance/Anova tables for taxonomic component

#Kenward-Roger
car::Anova(Recovery_SPEI3.tx.lmer, test = "F")
car::Anova(Recovery_SPEI12.tx.lmer, test = "F")
#Incremental F-test
car::Anova(Recovery_SPEI24.tx)


##----------------------------------------------------------------Resistance

##Functional components

##SPEI-3: Resistance_SPEI3.nlme
##SPEI-12: Resistance_SPEI12.nlme
##SPEI-24: Resistance_SPEI24.nlme

##Taxonomic component

##SPEI-3: Resistance_SPEI3.tx.nlme
##SPEI-12: Resistance_SPEI12.tx.nlme
##SPEI-24: Resistance_SPEI24.tx.nlme

#models for recovery functional components
ResiFun.list <- list(SPEI3 = Resistance_SPEI3.nlme, SPEI12 = Resistance_SPEI12.nlme, SPEI24 = Resistance_SPEI24.nlme)

#marg and cond R2 computed by modelsummary are same of those computed by performance::r2
lapply(ResiFun.list, r2)

#docx
modelsummary(models = ResiFun.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Resi_Fun.docx")

#models for recovery taxonomic component
ResiTx.list <- list(SPEI3 = Resistance_SPEI3.tx.nlme, SPEI12 = Resistance_SPEI12.tx.nlme,
                    SPEI24 = Resistance_SPEI24.tx.nlme)

#marg and cond R2 computed by modelsummary are same of those computed by performance::r2
lapply(ResiTx.list, r2)

#docx
modelsummary(models = ResiTx.list, coef_map = Coef_names_models, statistic = "{std.error} ({p.value})",
             gof_omit = 'AIC|BIC|ICC|RMSE',
             output = "~/Documents/ClimateExtremes/CleanCode/ModelSummary/Resi_Tx.docx")

#Deviance tables for functional components

#Kenward-Roger
car::Anova(Resistance_SPEI3.lmer, test = "F")
car::Anova(Resistance_SPEI12.lmer, test = "F")
car::Anova(Resistance_SPEI24.lmer, test = "F")

#Deviance tables for taxonomic component

#Kenward-Roger
car::Anova(Resistance_SPEI3.tx.lmer, test = "F")
car::Anova(Resistance_SPEI12.tx.lmer, test = "F")
car::Anova(Resistance_SPEI24.tx.lmer, test = "F")

