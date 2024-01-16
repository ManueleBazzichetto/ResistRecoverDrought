##Code to replicate results presented in the manuscript "Biodiversity promotes resistance but dominant species shape recovery of grasslands under extreme drought" by Bazzichetto et al.

##Statistical analyses - models testing the statistical interaction between functional composition and diversity, taxonomic diversity and drought intensity

library(dplyr)
#for plotting
library(ggplot2)
library(ggpubr)
library(patchwork)
library(cowplot)
library(magick)
#for modelling
library(car)
library(performance)
#plotting model results
library(effects)

##-------------------------------Fit models testing interaction between biodiversity and drought intensity


##----------------------------------------------------------------Recovery


##----------------------------------------------------------------SPEI-3

#join data.frame with annual log ratio and predictors for modelling to SPEI-3 categories
TA_LogR_SPEI3 <- dplyr::left_join(TA_biom.LogR, SPEI_cat.ls$SPEI3[c("Explo", "SPEI_final_cat", "Year")], by = c("Explo", "Year"))

#check years to be sure the whole time-series is included
unique(TA_LogR_SPEI3$Year)

#get rid of 2019, which has NAs for SPEI
TA_LogR_SPEI3 <- TA_LogR_SPEI3[which(TA_LogR_SPEI3$Year != "2019"), ]

#select only years for testing recovery at SPEI3
TA_LogR_SPEI3 <- TA_LogR_SPEI3[(TA_LogR_SPEI3$Year %in% Recovery_yrs$SPEI3$Yrs[Recovery_yrs$SPEI3$Explo == "ALB"] & TA_LogR_SPEI3$Explo == "ALB") |
                                 (TA_LogR_SPEI3$Year %in% Recovery_yrs$SPEI3$Yrs[Recovery_yrs$SPEI3$Explo == "HAI"] & TA_LogR_SPEI3$Explo == "HAI") |
                                 (TA_LogR_SPEI3$Year %in% Recovery_yrs$SPEI3$Yrs[Recovery_yrs$SPEI3$Explo == "SCH"] & TA_LogR_SPEI3$Explo == "SCH"), ]

#code SPEI_final_cat (categorical variable representing SPEI categories) so that it has NormalWet as reference level
TA_LogR_SPEI3$SPEI_final_cat <- factor(TA_LogR_SPEI3$SPEI_final_cat, levels = c("Normal", "ModerateWet", "ExtremeWet"))

#check
levels(TA_LogR_SPEI3$SPEI_final_cat)

#center all numeric predictors before fitting the model

#added on 27th Nov 2023 to get scaled factors (averages) for plots of partial effects
LogR_sc_fact_SPEI3 <- attr(scale(TA_LogR_SPEI3[!colnames(TA_LogR_SPEI3) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")], center = T, scale = F), "scaled:center")

TA_LogR_SPEI3 <- data.frame(TA_LogR_SPEI3[c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                            scale(TA_LogR_SPEI3[!colnames(TA_LogR_SPEI3) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  center = T, scale = F))

#check NAs
anyNA(TA_LogR_SPEI3)


##model for functional component

#fit model with lm
Recovery_SPEI3 <- lm(LogR ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                            day_of_year + SoilMoist10,
                          data = TA_LogR_SPEI3)


#test significance of interactions using Type II Anova
Anova(Recovery_SPEI3)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Recovery_SPEI3), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Recovery_SPEI3), multiline = T)

#diagnostics
check_model(Recovery_SPEI3)
r2(Recovery_SPEI3) #.117 (R2) / .076 (R2 adj)

##model for taxonomic component

#fit model with lm
Recovery_SPEI3.tx <- lm(LogR ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                             day_of_year + SoilMoist10,
                           data = TA_LogR_SPEI3)


Anova(Recovery_SPEI3.tx)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Recovery_SPEI3.tx), multiline = T)


#diagnostics
check_model(Recovery_SPEI3.tx)
r2(Recovery_SPEI3.tx) #.179 (R2) / .150 (R2 adj)


##----------------------------------------------------------------SPEI-12

TA_LogR_SPEI12 <- dplyr::left_join(TA_biom.LogR, SPEI_cat.ls$SPEI12[c("Explo", "SPEI_final_cat", "Year")], by = c("Explo", "Year"))

unique(TA_LogR_SPEI12$Year)

#get rid of 2019, which has NAs for SPEI
TA_LogR_SPEI12 <- TA_LogR_SPEI12[which(TA_LogR_SPEI12$Year != "2019"), ]

#select only years for testing recovery at SPEI12
TA_LogR_SPEI12 <- TA_LogR_SPEI12[(TA_LogR_SPEI12$Year %in% Recovery_yrs$SPEI12$Yrs[Recovery_yrs$SPEI12$Explo == "ALB"] & TA_LogR_SPEI12$Explo == "ALB") |
                                 (TA_LogR_SPEI12$Year %in% Recovery_yrs$SPEI12$Yrs[Recovery_yrs$SPEI12$Explo == "HAI"] & TA_LogR_SPEI12$Explo == "HAI") |
                                 (TA_LogR_SPEI12$Year %in% Recovery_yrs$SPEI12$Yrs[Recovery_yrs$SPEI12$Explo == "SCH"] & TA_LogR_SPEI12$Explo == "SCH"), ]

#code so that it has NormalWet as reference level
TA_LogR_SPEI12$SPEI_final_cat <- factor(TA_LogR_SPEI12$SPEI_final_cat, levels = c("Normal", "ModerateWet", "ExtremeWet"))
#check
levels(TA_LogR_SPEI12$SPEI_final_cat)

#center all numeric predictors before fitting the model

#added on 27th Nov 2023 to get scaled factors (averages) for plots of partial effects
LogR_sc_fact_SPEI12 <- attr(scale(TA_LogR_SPEI12[!colnames(TA_LogR_SPEI12) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")], center = T, scale = F), "scaled:center")

TA_LogR_SPEI12 <- data.frame(TA_LogR_SPEI12[c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                            scale(TA_LogR_SPEI12[!colnames(TA_LogR_SPEI12) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  center = T, scale = F))

#check NAs
anyNA(TA_LogR_SPEI12)

##model for functional component

#fit model with lm
Recovery_SPEI12 <- lm(LogR ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                             day_of_year + SoilMoist10,
                           data = TA_LogR_SPEI12)


Anova(Recovery_SPEI12)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Recovery_SPEI12), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Recovery_SPEI12), multiline = T)


#diagnostics
check_model(Recovery_SPEI12)
r2(Recovery_SPEI12) #.292 (R2) / .270 (R2 adj)

##model for taxonomic component

#fit model with lm
Recovery_SPEI12.tx <- lm(LogR ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                day_of_year + SoilMoist10,
                              data = TA_LogR_SPEI12)


Anova(Recovery_SPEI12.tx)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Recovery_SPEI12.tx), multiline = T)


#diagnostics
check_model(Recovery_SPEI12.tx)
r2(Recovery_SPEI12.tx) #.290 (R2) / .273 (R2 adj)

##----------------------------------------------------------------SPEI-24

TA_LogR_SPEI24 <- dplyr::left_join(TA_biom.LogR, SPEI_cat.ls$SPEI24[c("Explo", "SPEI_final_cat", "Year")], by = c("Explo", "Year"))

unique(TA_LogR_SPEI24$Year)

#get rid of 2019, which has NAs for SPEI
TA_LogR_SPEI24 <- TA_LogR_SPEI24[which(TA_LogR_SPEI24$Year != "2019"), ]

#select only years for testing recovery at SPEI24
TA_LogR_SPEI24 <- TA_LogR_SPEI24[(TA_LogR_SPEI24$Year %in% Recovery_yrs$SPEI24$Yrs[Recovery_yrs$SPEI24$Explo == "ALB"] & TA_LogR_SPEI24$Explo == "ALB") |
                                   (TA_LogR_SPEI24$Year %in% Recovery_yrs$SPEI24$Yrs[Recovery_yrs$SPEI24$Explo == "HAI"] & TA_LogR_SPEI24$Explo == "HAI") |
                                   (TA_LogR_SPEI24$Year %in% Recovery_yrs$SPEI24$Yrs[Recovery_yrs$SPEI24$Explo == "SCH"] & TA_LogR_SPEI24$Explo == "SCH"), ]

#code so that it has NormalWet as reference level
TA_LogR_SPEI24$SPEI_final_cat <- factor(TA_LogR_SPEI24$SPEI_final_cat, levels = c("Normal")) # -> only 1 level (Normal)
#check
levels(TA_LogR_SPEI24$SPEI_final_cat)

#center all numeric predictors before fitting the model

#added on 27th Nov 2023 to get scaled factors (averages) for plots of partial effects
LogR_sc_fact_SPEI24 <- attr(scale(TA_LogR_SPEI24[!colnames(TA_LogR_SPEI24) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")], center = T, scale = F), "scaled:center")

TA_LogR_SPEI24 <- data.frame(TA_LogR_SPEI24[c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                             scale(TA_LogR_SPEI24[!colnames(TA_LogR_SPEI24) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                   center = T, scale = F))

#check NAs
anyNA(TA_LogR_SPEI24)

##models for functional component

#only 1 level for SPEI_final_cat, so no interaction is fitted
#no repetitions of plotID (1 year/region), so lm is fitted

#fit model with lm
Recovery_SPEI24 <- lm(LogR ~ Across_PC1_CWM + Yr_Comb_rao_FP + Yr_lui + Explo +
                              day_of_year + SoilMoist10,
                            data = TA_LogR_SPEI24)

#Anova F-test
Anova(Recovery_SPEI24)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Recovery_SPEI24))
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Recovery_SPEI24))


#diagnostics
check_model(Recovery_SPEI24)
r2(Recovery_SPEI24) #0.356 (adj - 0.324)


##models for taxonomic component

#fit model with lm
Recovery_SPEI24.tx <- lm(LogR ~ Yr_species_rich + Yr_lui + Explo +
                                 day_of_year + SoilMoist10,
                               data = TA_LogR_SPEI24)

#Anova with F-test
Anova(Recovery_SPEI24.tx)

#check partial effect for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Recovery_SPEI24.tx))


#diagnostics
check_model(Recovery_SPEI24.tx)
r2(Recovery_SPEI24.tx) #0.358 (R2) / 0.331 (R2 adj)

##----------------------------------------------------------------Resistance


##----------------------------------------------------------------SPEI-3

TA_LogR_trpl_SPEI3 <- dplyr::left_join(TA_biom.LogR_tr_plot, SPEI_cat.ls$SPEI3[c("Explo", "SPEI_final_cat", "Year")], by = c("Explo", "Year"))

unique(TA_LogR_trpl_SPEI3$Year)

#get rid of 2019, which has NAs for SPEI
TA_LogR_trpl_SPEI3 <- TA_LogR_trpl_SPEI3[which(TA_LogR_trpl_SPEI3$Year != "2019"), ]

#select only years for testing resistance at SPEI3
TA_LogR_trpl_SPEI3 <- TA_LogR_trpl_SPEI3[(TA_LogR_trpl_SPEI3$Year %in% Resist_yrs$SPEI3$Year[Resist_yrs$SPEI3$Explo == "ALB"] & TA_LogR_trpl_SPEI3$Explo == "ALB") |
                                 (TA_LogR_trpl_SPEI3$Year %in% Resist_yrs$SPEI3$Year[Resist_yrs$SPEI3$Explo == "HAI"] & TA_LogR_trpl_SPEI3$Explo == "HAI") |
                                 (TA_LogR_trpl_SPEI3$Year %in% Resist_yrs$SPEI3$Year[Resist_yrs$SPEI3$Explo == "SCH"] & TA_LogR_trpl_SPEI3$Explo == "SCH"), ]

#code so that it has ModerateDry as reference level
TA_LogR_trpl_SPEI3$SPEI_final_cat <- factor(TA_LogR_trpl_SPEI3$SPEI_final_cat, levels = c("ModerateDry", "ExtremeDry"))
#check
levels(TA_LogR_trpl_SPEI3$SPEI_final_cat)


#center all numeric predictors before fitting the model

#added on 27th Nov 2023 to get scaled factors (averages) for plots of partial effects
LogRtrpl_sc_fact_SPEI3 <- attr(scale(TA_LogR_trpl_SPEI3[!colnames(TA_LogR_trpl_SPEI3) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")], center = T, scale = F), "scaled:center")

TA_LogR_trpl_SPEI3 <- data.frame(TA_LogR_trpl_SPEI3[c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                            scale(TA_LogR_trpl_SPEI3[!colnames(TA_LogR_trpl_SPEI3) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  center = T, scale = F))

#check NAs
anyNA(TA_LogR_trpl_SPEI3)


##model for functional component

#fit model with lm
Resistance_SPEI3 <- lm(LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                             day_of_year + SoilMoist10,
                           data = TA_LogR_trpl_SPEI3)


Anova(Resistance_SPEI3)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Resistance_SPEI3), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Resistance_SPEI3), multiline = T)


#diagnostics
check_model(Resistance_SPEI3)
r2(Resistance_SPEI3) #.134 (R2) / .123 (R2 adj)

##model for taxonomic component

#fit model with lm
Resistance_SPEI3.tx <- lm(LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                day_of_year + SoilMoist10,
                              data = TA_LogR_trpl_SPEI3)


Anova(Resistance_SPEI3.tx)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Resistance_SPEI3.tx), multiline = T)


#diagnostics
check_model(Resistance_SPEI3.tx)
r2(Resistance_SPEI3.tx) #.162 (R2) / .153 (R2 adj)

##----------------------------------------------------------------SPEI-12

TA_LogR_trpl_SPEI12 <- dplyr::left_join(TA_biom.LogR_tr_plot, SPEI_cat.ls$SPEI12[c("Explo", "SPEI_final_cat", "Year")], by = c("Explo", "Year"))

unique(TA_LogR_trpl_SPEI12$Year)

#get rid of 2019, which has NAs for SPEI
TA_LogR_trpl_SPEI12 <- TA_LogR_trpl_SPEI12[which(TA_LogR_trpl_SPEI12$Year != "2019"), ]

#select only years for testing resistance at SPEI12
TA_LogR_trpl_SPEI12 <- TA_LogR_trpl_SPEI12[(TA_LogR_trpl_SPEI12$Year %in% Resist_yrs$SPEI12$Year[Resist_yrs$SPEI12$Explo == "ALB"] & TA_LogR_trpl_SPEI12$Explo == "ALB") |
                                           (TA_LogR_trpl_SPEI12$Year %in% Resist_yrs$SPEI12$Year[Resist_yrs$SPEI12$Explo == "HAI"] & TA_LogR_trpl_SPEI12$Explo == "HAI") |
                                           (TA_LogR_trpl_SPEI12$Year %in% Resist_yrs$SPEI12$Year[Resist_yrs$SPEI12$Explo == "SCH"] & TA_LogR_trpl_SPEI12$Explo == "SCH"), ]

#code so that it has ModerateDry as reference level
TA_LogR_trpl_SPEI12$SPEI_final_cat <- factor(TA_LogR_trpl_SPEI12$SPEI_final_cat, levels = c("ModerateDry", "ExtremeDry"))
#check
levels(TA_LogR_trpl_SPEI12$SPEI_final_cat)


#center all numeric predictors before fitting the model

#added on 27th Nov 2023 to get scaled factors (averages) for plots of partial effects
LogRtrpl_sc_fact_SPEI12 <- attr(scale(TA_LogR_trpl_SPEI12[!colnames(TA_LogR_trpl_SPEI12) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")], center = T, scale = F), "scaled:center")

TA_LogR_trpl_SPEI12 <- data.frame(TA_LogR_trpl_SPEI12[c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                 scale(TA_LogR_trpl_SPEI12[!colnames(TA_LogR_trpl_SPEI12) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                       center = T, scale = F))

#check NAs
anyNA(TA_LogR_trpl_SPEI12)


##model for functional component

#fit model with lm
Resistance_SPEI12 <- lm(LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                               day_of_year + SoilMoist10,
                             data = TA_LogR_trpl_SPEI12)

Anova(Resistance_SPEI12)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Resistance_SPEI12), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Resistance_SPEI12), multiline = T)


#diagnostics
check_model(Resistance_SPEI12)
r2(Resistance_SPEI12) #.241 (R2) / .232 (R2 adj)

##model for taxonomic component

#fit model with lm
Resistance_SPEI12.tx <- lm(LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                  day_of_year + SoilMoist10,
                                data = TA_LogR_trpl_SPEI12)


Anova(Resistance_SPEI12.tx)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Resistance_SPEI12.tx), multiline = T)


#diagnostics
check_model(Resistance_SPEI12.tx)
r2(Resistance_SPEI12.tx) #.263 (R2) / .256 (R2 adj)


##----------------------------------------------------------------SPEI-24

TA_LogR_trpl_SPEI24 <- dplyr::left_join(TA_biom.LogR_tr_plot, SPEI_cat.ls$SPEI24[c("Explo", "SPEI_final_cat", "Year")], by = c("Explo", "Year"))

unique(TA_LogR_trpl_SPEI24$Year)

#get rid of 2019, which has NAs for SPEI
TA_LogR_trpl_SPEI24 <- TA_LogR_trpl_SPEI24[which(TA_LogR_trpl_SPEI24$Year != "2019"), ]

#select only years for testing resistance at SPEI24
TA_LogR_trpl_SPEI24 <- TA_LogR_trpl_SPEI24[(TA_LogR_trpl_SPEI24$Year %in% Resist_yrs$SPEI24$Year[Resist_yrs$SPEI24$Explo == "ALB"] & TA_LogR_trpl_SPEI24$Explo == "ALB") |
                                             (TA_LogR_trpl_SPEI24$Year %in% Resist_yrs$SPEI24$Year[Resist_yrs$SPEI24$Explo == "HAI"] & TA_LogR_trpl_SPEI24$Explo == "HAI") |
                                             (TA_LogR_trpl_SPEI24$Year %in% Resist_yrs$SPEI24$Year[Resist_yrs$SPEI24$Explo == "SCH"] & TA_LogR_trpl_SPEI24$Explo == "SCH"), ]

#code so that it has ModerateDry as reference level
TA_LogR_trpl_SPEI24$SPEI_final_cat <- factor(TA_LogR_trpl_SPEI24$SPEI_final_cat, levels = c("ModerateDry", "ExtremeDry"))
#check
levels(TA_LogR_trpl_SPEI24$SPEI_final_cat)


#center all numeric predictors before fitting the model

#added on 27th Nov 2023 to get scaled factors (averages) for plots of partial effects
LogRtrpl_sc_fact_SPEI24 <- attr(scale(TA_LogR_trpl_SPEI24[!colnames(TA_LogR_trpl_SPEI24) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")], center = T, scale = F), "scaled:center")

TA_LogR_trpl_SPEI24 <- data.frame(TA_LogR_trpl_SPEI24[c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  scale(TA_LogR_trpl_SPEI24[!colnames(TA_LogR_trpl_SPEI24) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                        center = T, scale = F))

#check NAs
anyNA(TA_LogR_trpl_SPEI24)


##model for functional component

#fit model with lm
Resistance_SPEI24 <- lm(LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                                day_of_year + SoilMoist10,
                              data = TA_LogR_trpl_SPEI24)


Anova(Resistance_SPEI24)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Resistance_SPEI24), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Resistance_SPEI24), multiline = T)


#diagnostics
check_model(Resistance_SPEI24)
r2(Resistance_SPEI24) #.254 (R2) / .245 (R2 adj) 

##model for taxonomic component

#fit model with lm
Resistance_SPEI24.tx <- lm(LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                   day_of_year + SoilMoist10,
                                 data = TA_LogR_trpl_SPEI24)

Anova(Resistance_SPEI24.tx)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Resistance_SPEI24.tx), multiline = T)


#diagnostics
check_model(Resistance_SPEI24.tx)
r2(Resistance_SPEI24.tx) #.257 (R2) / .250 (R2 adj)


##---------------------------Plots of confidence intervals for regression parameters

#create list including models for functional components

Functional_models <- list(RC_3 = Recovery_SPEI3, RC_12 = Recovery_SPEI12, RC_24 = Recovery_SPEI24,
                          RS_3 = Resistance_SPEI3, RS_12 = Resistance_SPEI12, RS_24 = Resistance_SPEI24)


#create table with confidence intervals for regression parameters of models including functional components

Fun_mod_CI <- do.call(rbind, lapply(names(Functional_models), function(nm) {
  fun_mod <- Functional_models[[nm]]
  fun_coef <- coef(fun_mod)
  fun_ci <- confint(fun_mod)
  #check estimates are associated with correct CI
  if(!isTRUE(identical(names(fun_coef), rownames(fun_ci)))) stop("Estimates' names not matching CI's names")
  ci_tab_sf <- data.frame(estim = unname(fun_coef[grepl(pattern = "Across_PC1_CWM", x = names(fun_coef))]),
                          fun_ci[grepl(pattern = "Across_PC1_CWM", x = rownames(fun_ci)), , drop = F],
                          predictor = "Across_PC1_CWM")
  ci_tab_fd <- data.frame(estim = unname(fun_coef[grepl(pattern = "Yr_Comb_rao_FP", x = names(fun_coef))]),
                          fun_ci[grepl(pattern = "Yr_Comb_rao_FP", x = rownames(fun_ci)), , drop = F],
                          predictor = "Yr_Comb_rao_FP")
  ci_tab <- rbind(ci_tab_sf, ci_tab_fd)
  colnames(ci_tab)[c(2, 3)] <- c("lower", "upper")
  ci_tab$SPEI_resp <- nm
  return(ci_tab)
}))

#create a new column reporting the SPEI category associated with each parameter's estimate and confidence interval
#this table will be used to generate plots of confidence intervals

Fun_mod_CI$SPEI_cat <- NA

Fun_mod_CI[grep("ModerateWet", row.names(Fun_mod_CI)), "SPEI_cat"] <- "ModerateWet"
Fun_mod_CI[grep("ExtremeWet", row.names(Fun_mod_CI)), "SPEI_cat"] <- "ExtremeWet"
Fun_mod_CI[grep("ExtremeDry", row.names(Fun_mod_CI)), "SPEI_cat"] <- "ExtremeDry"

Fun_mod_CI[is.na(Fun_mod_CI$SPEI_cat) & grepl("RC_", Fun_mod_CI$SPEI_resp), "SPEI_cat"] <- "Normal"
Fun_mod_CI[is.na(Fun_mod_CI$SPEI_cat), "SPEI_cat"] <- "ModerateDry"


#the same procedure is carried out for the models including species richness (the taxonomic component)

Taxonomic_models <- list(RC_3 = Recovery_SPEI3.tx, RC_12 = Recovery_SPEI12.tx, RC_24 = Recovery_SPEI24.tx,
                          RS_3 = Resistance_SPEI3.tx, RS_12 = Resistance_SPEI12.tx, RS_24 = Resistance_SPEI24.tx)


#create table with confidence intervals for regression parameters of models including the taxonomic component

Tax_mod_CI <- do.call(rbind, lapply(names(Taxonomic_models), function(nm) {
  tax_mod <- Taxonomic_models[[nm]]
  tax_coef <- coef(tax_mod)
  tax_ci <- confint(tax_mod)
  #check estimates are associated with correct CI
  if(!isTRUE(identical(names(tax_coef), rownames(tax_ci)))) stop("Estimates' names not matching CI's names")
  ci_tab_sr <- data.frame(estim = unname(tax_coef[grepl(pattern = "Yr_species_rich", x = names(tax_coef))]),
                          tax_ci[grepl(pattern = "Yr_species_rich", x = rownames(tax_ci)), , drop = F],
                          predictor = "Yr_species_rich")
  colnames(ci_tab_sr)[c(2, 3)] <- c("lower", "upper")
  ci_tab_sr$SPEI_resp <- nm
  return(ci_tab_sr)
}))


#create a new column reporting the SPEI category associated with each parameter's estimate and confidence interval
#this table will be used to generate plots of confidence intervals

Tax_mod_CI$SPEI_cat <- NA

Tax_mod_CI[grep("ModerateWet", row.names(Tax_mod_CI)), "SPEI_cat"] <- "ModerateWet"
Tax_mod_CI[grep("ExtremeWet", row.names(Tax_mod_CI)), "SPEI_cat"] <- "ExtremeWet"
Tax_mod_CI[grep("ExtremeDry", row.names(Tax_mod_CI)), "SPEI_cat"] <- "ExtremeDry"

Tax_mod_CI[is.na(Tax_mod_CI$SPEI_cat) & grepl("RC_", Tax_mod_CI$SPEI_resp), "SPEI_cat"] <- "Normal"
Tax_mod_CI[is.na(Tax_mod_CI$SPEI_cat), "SPEI_cat"] <- "ModerateDry"


#put together CI(s) for recovery and resistance

Recovery_CIs <- rbind(Fun_mod_CI[grep("RC_", x = Fun_mod_CI$SPEI_resp), ],
                      Tax_mod_CI[grep("RC_", x = Tax_mod_CI$SPEI_resp), ])


Resistance_CIs <- rbind(Fun_mod_CI[grep("RS_", x = Fun_mod_CI$SPEI_resp), ],
                        Tax_mod_CI[grep("RS_", x = Tax_mod_CI$SPEI_resp), ])


#create plots for models of recovery (last update for review: Nov. 21st 2023)

#order SPEI_resp for increasing SPEI time-scale
Recovery_CIs$SPEI_resp <- factor(Recovery_CIs$SPEI_resp, levels = unique(Recovery_CIs$SPEI_resp))

#order SPEI_cat for increasing water availability
Recovery_CIs$SPEI_cat <- factor(Recovery_CIs$SPEI_cat, levels = unique(Recovery_CIs$SPEI_cat))

Recovery_plot <- ggplot(Recovery_CIs, aes(x = SPEI_cat, y = estim, col = SPEI_cat)) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  geom_point(cex = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  scale_color_manual(values = Drought_pal[c(3, 4, 5)], name = "SPEI category") +
  facet_grid(predictor ~ SPEI_resp, scales = "free",
             labeller = labeller(SPEI_resp = c("RC_3" = "SPEI-3",
                                               "RC_12" = "SPEI-12",
                                               "RC_24" = "SPEI-24"), 
                                 predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness"))) +
  ylab("LogR - Recovery, model coefficient (95% CI)") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "bottom", axis.text.x.bottom = element_blank(),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = "black"))

ggsave(plot = Recovery_plot, "~/Documents/ClimateExtremes/CleanCode/New_figures/Recovery_CI.jpeg", device = "jpeg", dpi = 300, width = 22, height = 20, units = "cm")

#create plots for models of resistance (last update for review: Nov. 21st 2023)

#order SPEI_Resp for increasing SPEI time-scale
Resistance_CIs$SPEI_resp <- factor(Resistance_CIs$SPEI_resp, levels = unique(Resistance_CIs$SPEI_resp))

#order SPEI_cat for increasing water availability
Resistance_CIs$SPEI_cat <- factor(Resistance_CIs$SPEI_cat, levels = rev(unique(Resistance_CIs$SPEI_cat)))

Resistance_plot <- ggplot(Resistance_CIs, aes(x = SPEI_cat, y = estim, col = SPEI_cat)) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  geom_point(cex = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  scale_color_manual(values = Drought_pal[c(1, 2)], name = "SPEI category") +
  facet_grid(predictor ~ SPEI_resp, scales = "free",
             labeller = labeller(SPEI_resp = c("RS_3" = "SPEI-3",
                                               "RS_12" = "SPEI-12",
                                               "RS_24" = "SPEI-24"), 
                                 predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness"))) +
  ylab("LogR ref-plot - Resistance, model coefficient (95% CI)") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "bottom", axis.text.x.bottom = element_blank(),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 12),
        panel.background = element_rect(fill = NA, color = "black"))

ggsave(plot = Resistance_plot, "~/Documents/ClimateExtremes/CleanCode/New_figures/Resistance_CI.jpeg", device = "jpeg", dpi = 300, width = 22, height = 20, units = "cm")


##---------------Prediction plots showing the (partial) effects of functional and taxonomic components

#get predictions for models including functional components
#predictions for the different SPEI time-scales are then put together in a single data.frame

Fun_mod_fitted <- do.call(rbind, lapply(names(Functional_models), function(nm) {
  Fun_mod <- Functional_models[[nm]]
  Slow_fast.fit <- as.data.frame(predictorEffect(predictor = "Across_PC1_CWM", mod = Fun_mod))
  Slow_fast.fit <- data.frame(Slow_fast.fit, predictor = "Across_PC1_CWM", SPEI_resp = nm)
  Fd.fit <- as.data.frame(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Fun_mod))
  Fd.fit <- data.frame(Fd.fit, predictor = "Yr_Comb_rao_FP", SPEI_resp = nm)
  colnames(Slow_fast.fit)[1] <- colnames(Fd.fit)[1] <- "X"
  Df <- rbind(Slow_fast.fit, Fd.fit)
  if(nm %in% "RC_24") {
    Df$SPEI_final_cat <- "Normal"
    Df <- Df[c(1, 8, 2, 3, 4, 5, 6, 7)]
  }
  return(Df)
}))


#same as above but for models including the taxonomic component (species richness)

Tax_mod_fitted <- do.call(rbind, lapply(names(Taxonomic_models), function(nm) {
  Tax_mod <- Taxonomic_models[[nm]]
  Sr.fit <- as.data.frame(predictorEffect(predictor = "Yr_species_rich", mod = Tax_mod))
  Sr.fit <- data.frame(Sr.fit, predictor = "Yr_species_rich", SPEI_resp = nm)
  colnames(Sr.fit)[1] <- "X"
  if(nm %in% "RC_24") {
    Sr.fit$SPEI_final_cat <- "Normal"
    Sr.fit <- Sr.fit[c(1, 8, 2, 3, 4, 5, 6, 7)]
  }
  return(Sr.fit)
}))


#put together data.frames of predictions for the models of recovery and resistance

Recovery_FIT <- rbind(Fun_mod_fitted[grep("RC_", x = Fun_mod_fitted$SPEI_resp), ],
                      Tax_mod_fitted[grep("RC_", x = Tax_mod_fitted$SPEI_resp), ])


Resistance_FIT <- rbind(Fun_mod_fitted[grep("RS_", x = Fun_mod_fitted$SPEI_resp), ],
                        Tax_mod_fitted[grep("RS_", x = Tax_mod_fitted$SPEI_resp), ])


#prediction plots for models of recovery (last update for review: Nov. 20th 2023)

#for changing position of facets' strip text
#see: https://stackoverflow.com/questions/69199276/ggplot2-facet-wrap-graph-with-custom-x-axis-labels

#order SPEI_final_cat for increasing water availability
Recovery_FIT$SPEI_final_cat <- factor(Recovery_FIT$SPEI_final_cat, levels = as.character(unique(Recovery_FIT$SPEI_final_cat)))

#SPEI-3 - put apart predictions at SPEI-3 as this will serve to create a separate figure
Recovery_FIT_SP3 <- Recovery_FIT[Recovery_FIT$SPEI_resp == "RC_3", ]

#SPEI-12 & -24
Recovery_FIT <- Recovery_FIT[Recovery_FIT$SPEI_resp != "RC_3", ]

#order factor for increasing SPEI time-scale
Recovery_FIT$SPEI_resp <- factor(Recovery_FIT$SPEI_resp, levels = unique(Recovery_FIT$SPEI_resp))

#transform predictors back to the original scale (sum average used as scaling factor when centering)

#add an empty column, which will be filled with original values for the predictors
Recovery_FIT$X_original <- NA

#SPEI-12
#all(as.character(Recovery_FIT$Predictor[Recovery_FIT$SPEI_resp == "RC_12"]) %in% names(LogR_sc_fact_SPEI12)) #T

Recovery_FIT$X_original[Recovery_FIT$SPEI_resp == "RC_12"] <- (Recovery_FIT$X[Recovery_FIT$SPEI_resp == "RC_12"] + unname(LogR_sc_fact_SPEI12[as.character(Recovery_FIT$predictor[Recovery_FIT$SPEI_resp == "RC_12"])]))

#SPEI-24
#all(as.character(Recovery_FIT$Predictor[Recovery_FIT$SPEI_resp == "RC_24"]) %in% names(LogR_sc_fact_SPEI24)) #T

Recovery_FIT$X_original[Recovery_FIT$SPEI_resp == "RC_24"] <- (Recovery_FIT$X[Recovery_FIT$SPEI_resp == "RC_24"] + unname(LogR_sc_fact_SPEI24[as.character(Recovery_FIT$predictor[Recovery_FIT$SPEI_resp == "RC_24"])]))


Recovery_partial_plot <- ggplot(Recovery_FIT, aes(x = X_original, y = fit, fill = SPEI_final_cat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  geom_line(aes(col = SPEI_final_cat), lwd = 1) +
  scale_color_manual(values = Drought_pal[c(3, 4, 5)], name = "SPEI category") +
  scale_fill_manual(values = Drought_pal[c(3, 4, 5)], name = "SPEI category") +
  facet_grid(SPEI_resp ~ predictor, scales = "free",
             labeller = labeller(SPEI_resp = c("RC_12" = "SPEI-12",
                                               "RC_24" = "SPEI-24"), 
                                 predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness")), switch = "x") +
  ylab("LogR - Recovery") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "top", axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16), strip.placement = "outside", strip.background.x = element_blank(),
        legend.text = element_text(size = 14), legend.title = element_text(size = 12))


ggsave(plot = Recovery_partial_plot, "~/Documents/ClimateExtremes/CleanCode/New_figures/Recovery_partial.jpeg", device = "jpeg", dpi = 300, width = 24, height = 22, units = "cm")


#prediction plots for models of resistance (last update for review: Nov. 20th 2023)

#order SPEI_final_cat for increasing water availability
Resistance_FIT$SPEI_final_cat <- factor(Resistance_FIT$SPEI_final_cat, levels = rev(as.character(unique(Resistance_FIT$SPEI_final_cat))))

#SPEI-3
Resistance_FIT_SP3 <- Resistance_FIT[Resistance_FIT$SPEI_resp == "RS_3", ]

#SPEI-12 & -24
Resistance_FIT <- Resistance_FIT[Resistance_FIT$SPEI_resp != "RS_3", ]

#order SPEI_resp for increasing SPEI time-scale
Resistance_FIT$SPEI_resp <- factor(Resistance_FIT$SPEI_resp, levels = unique(Resistance_FIT$SPEI_resp))

#transform predictors back to original scale (sum average used as scaling factor when centering)

Resistance_FIT$X_original <- NA

#SPEI-12
Resistance_FIT$X_original[Resistance_FIT$SPEI_resp == "RS_12"] <- (Resistance_FIT$X[Resistance_FIT$SPEI_resp == "RS_12"] +
                                                                     unname(LogRtrpl_sc_fact_SPEI12[Resistance_FIT$predictor[Resistance_FIT$SPEI_resp == "RS_12"]]))

#SPEI-24
Resistance_FIT$X_original[Resistance_FIT$SPEI_resp == "RS_24"] <- (Resistance_FIT$X[Resistance_FIT$SPEI_resp == "RS_24"] +
                                                                     unname(LogRtrpl_sc_fact_SPEI24[Resistance_FIT$predictor[Resistance_FIT$SPEI_resp == "RS_24"]]))

Resistance_partial_plot <- ggplot(Resistance_FIT, aes(x = X_original, y = fit, fill = SPEI_final_cat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  geom_line(aes(col = SPEI_final_cat), lwd = 1) +
  scale_color_manual(values = Drought_pal[c(1, 2)], name = "SPEI category") +
  scale_fill_manual(values = Drought_pal[c(1, 2)], name = "SPEI category") +
  facet_grid(SPEI_resp ~ predictor, scales = "free",
             labeller = labeller(SPEI_resp = c("RS_12" = "SPEI-12",
                                               "RS_24" = "SPEI-24"), 
                                 predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness")), switch = "x") +
  ylab("LogR ref-plot - Resistance") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "top", axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16), strip.placement = "outside", strip.background.x = element_blank(),
        legend.text = element_text(size = 14), legend.title = element_text(size = 12))

ggsave(plot = Resistance_partial_plot, "~/Documents/ClimateExtremes/CleanCode/New_figures/Resistance_partial.jpeg", device = "jpeg", dpi = 300, width = 24, height = 22, units = "cm")


#prediction plots for models of recovery and resistance for SPEI-3

#put together predictions for recovery and resistance
Reco_resi_SP3 <- rbind(Recovery_FIT_SP3, Resistance_FIT_SP3)

#re-code SPEI_final_cat for increasing water availability
Reco_resi_SP3$SPEI_final_cat <- factor(Reco_resi_SP3$SPEI_final_cat, levels = names(Drought_pal))

#transform predictors back to original scale (sum average used as scaling factor when centering)

Reco_resi_SP3$X_original <- NA


Reco_resi_SP3$X_original[Reco_resi_SP3$SPEI_resp == "RC_3"] <- (Reco_resi_SP3$X[Reco_resi_SP3$SPEI_resp == "RC_3"] +
                                                                  unname(LogR_sc_fact_SPEI3[Reco_resi_SP3$predictor[Reco_resi_SP3$SPEI_resp == "RC_3"]])) 

Reco_resi_SP3$X_original[Reco_resi_SP3$SPEI_resp == "RS_3"] <- (Reco_resi_SP3$X[Reco_resi_SP3$SPEI_resp == "RS_3"] +
                                                                  unname(LogRtrpl_sc_fact_SPEI3[Reco_resi_SP3$predictor[Reco_resi_SP3$SPEI_resp == "RS_3"]]))


Reco_resi_plot <- ggplot(Reco_resi_SP3, aes(x = X_original, y = fit, fill = SPEI_final_cat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .3) +
  geom_line(aes(col = SPEI_final_cat), lwd = 1) +
  scale_color_manual(values = Drought_pal, name = "SPEI category") +
  scale_fill_manual(values = Drought_pal, name = "SPEI category") +
  facet_grid(SPEI_resp ~ predictor, scales = "free",
             labeller = labeller(SPEI_resp = c("RC_3" = "LogR - Recovery",
                                               "RS_3" = "LogR ref-plot - Resistance"), 
                                 predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness")), switch = "both") +
  ylab(NULL) + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "top", axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16), strip.placement = "outside", strip.background = element_blank(),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))


Reco_resi_plot <- ggdraw(Reco_resi_plot) +
  draw_image(image = Icon_eupesu, x = 0.29, y = 0.075, hjust = 0, valign = 0, width = 0.12, height = 0.12) +
  draw_image(image = Icon_fesspp, x = 0.04, y = 0.055, hjust = 0, valign = 0, width = 0.12, height = 0.11)


ggsave(plot = Reco_resi_plot, "~/Documents/ClimateExtremes/CleanCode/New_figures/Reco_resi_plot.jpeg", device = "jpeg", dpi = 500, width = 29, height = 25, units = "cm")



