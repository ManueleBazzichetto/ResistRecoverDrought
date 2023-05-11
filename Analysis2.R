##Code to replicate results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought" by Bazzichetto et al.
##Statistical analyses - models testing the statistical interaction between functional composition and diversity, taxonomic diversity and drought intensity

library(dplyr)
#for plotting
library(ggplot2)
library(ggpubr)
library(patchwork)
#library(GGally)
#for modelling
library(car)
#library(MuMIn)
#library(MASS)
library(performance)
library(nlme) 
library(lme4)
library(lmerTest)
#plotting model results
library(effects)

##----------------------------------------------------------------Fit models testing interaction between biodiversity and drought intensity

#compute minmax day of biomass cut
Min_day_harv <- with(biomass, tapply(X = day_of_year, INDEX = list(Explo, Year), min, na.rm = T))
round(rowMeans(Min_day_harv), digits = 0)

Max_day_harv <- with(biomass, tapply(X = day_of_year, INDEX = list(Explo, Year), max, na.rm = T))
round(rowMeans(Max_day_harv), digits = 0)

#create data.frame with predictors
TA_biom.LogRatios <- biomass[!(biomass$Year %in% c("2020")),
                                 c("LogR", "LogR_tr_plot", "Explo",
                                   "Yr_lui", "Across_PC1_CWM",
                                   "Yr_Comb_rao_FP", "Yr_Myco_intensity",
                                   "Yr_species_rich",
                                   "Year", "Useful_EP_PlotID", "day_of_year")]

#NAs for LogR
sum(is.na(TA_biom.LogRatios$LogR)); which(is.na(TA_biom.LogRatios$LogR) & !TA_biom.LogRatios$Year %in% c("2009"))
#NAs for LogR_tr_plot
sum(is.na(TA_biom.LogRatios$LogR_tr_plot)); which(is.na(TA_biom.LogRatios$LogR_tr_plot))
#NAs for predictors
lapply(TA_biom.LogRatios[c(4, 5, 6, 7, 8, 11)], function(i) which(is.na(i)))

#exclude NAs for Yr_lui (this will exclude NAs for all other predictors, but day_of_year - NAs are all at the same position)
TA_biom.LogRatios <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$Yr_lui), ]

#remove NAs for day_of_year -> these are only 3 for day_of_year
TA_biom.LogRatios <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$day_of_year), ]

#generate LogR specific datasets

#dataset for LogR (2009 is excluded)
TA_biom.LogR <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$LogR), !colnames(TA_biom.LogRatios) %in% c("LogR_tr_plot")]

#re-order df: alphab. for plots, chronolog. for year
TA_biom.LogR <- do.call(rbind, lapply(sort(unique(TA_biom.LogR$Useful_EP_PlotID)), function(nm) {
  df.plot <- TA_biom.LogR[TA_biom.LogR$Useful_EP_PlotID == nm, ]
  df.plot <- df.plot[order(df.plot$Year), ]
  return(df.plot)
})
)

#check
anyNA(TA_biom.LogR)
unique(TA_biom.LogR$Year)
unique(TA_biom.LogR$Useful_EP_PlotID)

#dataset for LogR_tr_plot
TA_biom.LogR_tr_plot <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$LogR_tr_plot), !colnames(TA_biom.LogRatios) %in% c("LogR")]

#re-order df: alphab. for plots, chronolog. for year
TA_biom.LogR_tr_plot <- do.call(rbind, lapply(sort(unique(TA_biom.LogR_tr_plot$Useful_EP_PlotID)), function(nm) {
  Df <- TA_biom.LogR_tr_plot[TA_biom.LogR_tr_plot$Useful_EP_PlotID == nm, ]
  Df <- Df[order(Df$Year), ]
  return(Df)
})
)

#check
anyNA(TA_biom.LogR_tr_plot)
unique(TA_biom.LogR_tr_plot$Year)
unique(TA_biom.LogR_tr_plot$Useful_EP_PlotID)


##----------------------------------------------------------------Recovery


##----------------------------------------------------------------SPEI-3

TA_LogR_SPEI3 <- dplyr::left_join(TA_biom.LogR, SPEI_cat.ls$SPEI3[c("Explo", "SPEI_final_cat", "Year")], by = c("Explo", "Year"))

unique(TA_LogR_SPEI3$Year)

#get rid of 2019, which has NAs for SPEI
TA_LogR_SPEI3 <- TA_LogR_SPEI3[which(TA_LogR_SPEI3$Year != "2019"), ]

#select only years for testing recovery at SPEI3
TA_LogR_SPEI3 <- TA_LogR_SPEI3[(TA_LogR_SPEI3$Year %in% Recovery_yrs$SPEI3$Yrs[Recovery_yrs$SPEI3$Explo == "ALB"] & TA_LogR_SPEI3$Explo == "ALB") |
                                 (TA_LogR_SPEI3$Year %in% Recovery_yrs$SPEI3$Yrs[Recovery_yrs$SPEI3$Explo == "HAI"] & TA_LogR_SPEI3$Explo == "HAI") |
                                 (TA_LogR_SPEI3$Year %in% Recovery_yrs$SPEI3$Yrs[Recovery_yrs$SPEI3$Explo == "SCH"] & TA_LogR_SPEI3$Explo == "SCH"), ]

#code so that it has NormalWet as reference level
TA_LogR_SPEI3$SPEI_final_cat <- factor(TA_LogR_SPEI3$SPEI_final_cat, levels = c("Normal", "ModerateWet", "ExtremeWet"))
#check
levels(TA_LogR_SPEI3$SPEI_final_cat)


#center all numeric predictors before fitting the model
TA_LogR_SPEI3 <- data.frame(TA_LogR_SPEI3[c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                            scale(TA_LogR_SPEI3[!colnames(TA_LogR_SPEI3) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  center = T, scale = F))

#check NAs
anyNA(TA_LogR_SPEI3)


##models for functional component

#fit model with nlme
Recovery_SPEI3.nlme <- lme(fixed = LogR ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                            day_of_year + Yr_Myco_intensity,
                          data = TA_LogR_SPEI3, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Recovery_SPEI3.lmer <- lmer(LogR ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                                  day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                                data = TA_LogR_SPEI3)


#Anova with corrected p-val Kenward-Roger
Anova(Recovery_SPEI3.lmer, test = "F")
#Anova with Wald's
Anova(Recovery_SPEI3.lmer)
#same but nlme
Anova(Recovery_SPEI3.nlme)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Recovery_SPEI3.nlme), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Recovery_SPEI3.nlme), multiline = T)

#diagnostics
plot(Recovery_SPEI3.nlme)
qqPlot(residuals(Recovery_SPEI3.nlme, type = "normalized"))
vif(Recovery_SPEI3.nlme)
performance::r2(Recovery_SPEI3.nlme) #.126 (marginal) / .126 (conditional)


##models for taxonomic component

#fit model with nlme
Recovery_SPEI3.tx.nlme <- lme(fixed = LogR ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                             day_of_year + Yr_Myco_intensity,
                           data = TA_LogR_SPEI3, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Recovery_SPEI3.tx.lmer <- lmer(LogR ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                              day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                            data = TA_LogR_SPEI3)


#Anova with corrected p-val Kenward-Roger
Anova(Recovery_SPEI3.tx.lmer, test = "F")
#Anova with Wald's
Anova(Recovery_SPEI3.tx.lmer)
#same but nlme
Anova(Recovery_SPEI3.tx.nlme)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Recovery_SPEI3.tx.nlme), multiline = T)

#diagnostics
plot(Recovery_SPEI3.tx.nlme)
qqPlot(residuals(Recovery_SPEI3.tx.nlme, type = "normalized"))
vif(Recovery_SPEI3.tx.nlme)
performance::r2(Recovery_SPEI3.tx.nlme) #.185 (marginal) / .185 (conditional)


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
TA_LogR_SPEI12 <- data.frame(TA_LogR_SPEI12[c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                            scale(TA_LogR_SPEI12[!colnames(TA_LogR_SPEI12) %in% c("LogR", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  center = T, scale = F))

#check NAs
anyNA(TA_LogR_SPEI12)

##models for functional component

#fit model with nlme
Recovery_SPEI12.nlme <- lme(fixed = LogR ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                             day_of_year + Yr_Myco_intensity,
                           data = TA_LogR_SPEI12, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Recovery_SPEI12.lmer <- lmer(LogR ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                              day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                            data = TA_LogR_SPEI12)


#Anova with corrected p-val Kenward-Roger
Anova(Recovery_SPEI12.lmer, test = "F")
#Anova with Wald's
Anova(Recovery_SPEI12.lmer)
#same but nlme
Anova(Recovery_SPEI12.nlme)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Recovery_SPEI12.nlme), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Recovery_SPEI12.nlme), multiline = T)

#diagnostics
plot(Recovery_SPEI12.nlme)
qqPlot(residuals(Recovery_SPEI12.nlme, type = "normalized"))
vif(Recovery_SPEI12.nlme)
performance::r2(Recovery_SPEI12.nlme) #.286 (marginal) / .286 (conditional)


##models for taxonomic component

#fit model with nlme
Recovery_SPEI12.tx.nlme <- lme(fixed = LogR ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                day_of_year + Yr_Myco_intensity,
                              data = TA_LogR_SPEI12, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Recovery_SPEI12.tx.lmer <- lmer(LogR ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                 day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                               data = TA_LogR_SPEI12)


#Anova with corrected p-val Kenward-Roger
Anova(Recovery_SPEI12.tx.lmer, test = "F")
#Anova with Wald's
Anova(Recovery_SPEI12.tx.lmer)
#same but nlme
Anova(Recovery_SPEI12.tx.nlme)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Recovery_SPEI12.tx.nlme), multiline = T)

#diagnostics
plot(Recovery_SPEI12.tx.nlme)
qqPlot(residuals(Recovery_SPEI12.tx.nlme, type = "normalized"))
vif(Recovery_SPEI12.tx.nlme)
performance::r2(Recovery_SPEI12.tx.nlme) #.287 (marginal) / .287 (conditional)


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
                              day_of_year + Yr_Myco_intensity,
                            data = TA_LogR_SPEI24)

#Anova F-test
Anova(Recovery_SPEI24)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Recovery_SPEI24))
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Recovery_SPEI24))

#diagnostics
plot(Recovery_SPEI24)
qqPlot(Recovery_SPEI24)
vif(Recovery_SPEI24)
performance::r2(Recovery_SPEI24) #0.345 (adj - 0.313)


##models for taxonomic component

#fit model with lm
Recovery_SPEI24.tx <- lm(LogR ~ Yr_species_rich + Yr_lui + Explo +
                                 day_of_year + Yr_Myco_intensity,
                               data = TA_LogR_SPEI24)

#Anova with F-test
Anova(Recovery_SPEI24.tx)

#check partial effect for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Recovery_SPEI24.tx))

#diagnostics
plot(Recovery_SPEI24.tx)
qqPlot(Recovery_SPEI24.tx)
vif(Recovery_SPEI24.tx)
performance::r2(Recovery_SPEI24.tx) #0.350 (adj 0.323)


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
TA_LogR_trpl_SPEI3 <- data.frame(TA_LogR_trpl_SPEI3[c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                            scale(TA_LogR_trpl_SPEI3[!colnames(TA_LogR_trpl_SPEI3) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  center = T, scale = F))

#check NAs
anyNA(TA_LogR_trpl_SPEI3)


##models for functional component

#fit model with nlme
Resistance_SPEI3.nlme <- lme(fixed = LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                             day_of_year + Yr_Myco_intensity,
                           data = TA_LogR_trpl_SPEI3, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Resistance_SPEI3.lmer <- lmer(LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                              day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                            data = TA_LogR_trpl_SPEI3)


#Anova with corrected p-val Kenward-Roger
Anova(Resistance_SPEI3.lmer, test = "F")
#Anova with Wald's
Anova(Resistance_SPEI3.lmer)
#same but nlme
Anova(Resistance_SPEI3.nlme)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Resistance_SPEI3.nlme), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Resistance_SPEI3.nlme), multiline = T)

#diagnostics
plot(Resistance_SPEI3.nlme)
qqPlot(residuals(Resistance_SPEI3.nlme, type = "normalized"))
vif(Resistance_SPEI3.nlme)
performance::r2(Resistance_SPEI3.nlme) #.130 (marginal) / .130 (conditional)


##models for taxonomic component

#fit model with nlme
Resistance_SPEI3.tx.nlme <- lme(fixed = LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                day_of_year + Yr_Myco_intensity,
                              data = TA_LogR_trpl_SPEI3, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Resistance_SPEI3.tx.lmer <- lmer(LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                 day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                               data = TA_LogR_trpl_SPEI3)


#Anova with corrected p-val Kenward-Roger
Anova(Resistance_SPEI3.tx.lmer, test = "F")
#Anova with Wald's
Anova(Resistance_SPEI3.tx.lmer)
#same but nlme
Anova(Resistance_SPEI3.tx.nlme)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Resistance_SPEI3.tx.nlme), multiline = T)

#diagnostics
plot(Resistance_SPEI3.tx.nlme)
qqPlot(residuals(Resistance_SPEI3.tx.nlme, type = "normalized"))
vif(Resistance_SPEI3.tx.nlme)
performance::r2(Resistance_SPEI3.tx.nlme) #.156 (marginal) / .156 (conditional)


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
TA_LogR_trpl_SPEI12 <- data.frame(TA_LogR_trpl_SPEI12[c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                 scale(TA_LogR_trpl_SPEI12[!colnames(TA_LogR_trpl_SPEI12) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                       center = T, scale = F))

#check NAs
anyNA(TA_LogR_trpl_SPEI12)


##models for functional component

#fit model with nlme
Resistance_SPEI12.nlme <- lme(fixed = LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                               day_of_year + Yr_Myco_intensity,
                             data = TA_LogR_trpl_SPEI12, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Resistance_SPEI12.lmer <- lmer(LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                                day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                              data = TA_LogR_trpl_SPEI12)


#Anova with corrected p-val Kenward-Roger
Anova(Resistance_SPEI12.lmer, test = "F")
#Anova with Wald's
Anova(Resistance_SPEI12.lmer)
#same but nlme
Anova(Resistance_SPEI12.nlme)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Resistance_SPEI12.nlme), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Resistance_SPEI12.nlme), multiline = T)

#diagnostics
plot(Resistance_SPEI12.nlme)
qqPlot(residuals(Resistance_SPEI12.nlme, type = "normalized"))
vif(Resistance_SPEI12.nlme)
performance::r2(Resistance_SPEI12.nlme) #.240 (marginal) / .240 (conditional)


##models for taxonomic component

#fit model with nlme
Resistance_SPEI12.tx.nlme <- lme(fixed = LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                  day_of_year + Yr_Myco_intensity,
                                data = TA_LogR_trpl_SPEI12, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Resistance_SPEI12.tx.lmer <- lmer(LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                   day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                                 data = TA_LogR_trpl_SPEI12)


#Anova with corrected p-val Kenward-Roger
Anova(Resistance_SPEI12.tx.lmer, test = "F")
#Anova with Wald's
Anova(Resistance_SPEI12.tx.lmer)
#same but nlme
Anova(Resistance_SPEI12.tx.nlme)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Resistance_SPEI12.tx.nlme), multiline = T)

#diagnostics
plot(Resistance_SPEI12.tx.nlme)
qqPlot(residuals(Resistance_SPEI12.tx.nlme, type = "normalized"))
vif(Resistance_SPEI12.tx.nlme)
performance::r2(Resistance_SPEI12.tx.nlme) #0.260 (marginal) / 0.260 (conditional)


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
TA_LogR_trpl_SPEI24 <- data.frame(TA_LogR_trpl_SPEI24[c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                  scale(TA_LogR_trpl_SPEI24[!colnames(TA_LogR_trpl_SPEI24) %in% c("LogR_tr_plot", "Explo", "Useful_EP_PlotID", "Year", "SPEI_final_cat")],
                                        center = T, scale = F))

#check NAs
anyNA(TA_LogR_trpl_SPEI24)


##models for functional component

#fit model with nlme
Resistance_SPEI24.nlme <- lme(fixed = LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                                day_of_year + Yr_Myco_intensity,
                              data = TA_LogR_trpl_SPEI24, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Resistance_SPEI24.lmer <- lmer(LogR_tr_plot ~ SPEI_final_cat*(Across_PC1_CWM + Yr_Comb_rao_FP) + Yr_lui + Explo +
                                 day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                               data = TA_LogR_trpl_SPEI24)


#Anova with corrected p-val Kenward-Roger
Anova(Resistance_SPEI24.lmer, test = "F")
#Anova with Wald's
Anova(Resistance_SPEI24.lmer)
#same but nlme
Anova(Resistance_SPEI24.nlme)

#check interaction for slow-fast
plot(predictorEffect(predictor = "Across_PC1_CWM", mod = Resistance_SPEI24.nlme), multiline = T)
#check interaction for functional diversity
plot(predictorEffect(predictor = "Yr_Comb_rao_FP", mod = Resistance_SPEI24.nlme), multiline = T)

#diagnostics
plot(Resistance_SPEI24.nlme)
qqPlot(residuals(Resistance_SPEI24.nlme, type = "normalized"))
vif(Resistance_SPEI24.nlme)
performance::r2(Resistance_SPEI24.nlme) #.252 (marginal) / .252 (conditional)


##models for taxonomic component

#fit model with nlme
Resistance_SPEI24.tx.nlme <- lme(fixed = LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                   day_of_year + Yr_Myco_intensity,
                                 data = TA_LogR_trpl_SPEI24, random = ~ 1 | Useful_EP_PlotID)

#fit model with lmer
Resistance_SPEI24.tx.lmer <- lmer(LogR_tr_plot ~ SPEI_final_cat*Yr_species_rich + Yr_lui + Explo +
                                    day_of_year + Yr_Myco_intensity + (1 | Useful_EP_PlotID),
                                  data = TA_LogR_trpl_SPEI24)


#Anova with corrected p-val Kenward-Roger
Anova(Resistance_SPEI24.tx.lmer, test = "F")
#Anova with Wald's
Anova(Resistance_SPEI24.tx.lmer)
#same but nlme
Anova(Resistance_SPEI24.tx.nlme)

#check interaction for species richness
plot(predictorEffect(predictor = "Yr_species_rich", mod = Resistance_SPEI24.tx.nlme), multiline = T)

#diagnostics
plot(Resistance_SPEI24.tx.nlme)
qqPlot(residuals(Resistance_SPEI24.tx.nlme, type = "normalized"))
vif(Resistance_SPEI24.tx.nlme)
performance::r2(Resistance_SPEI24.tx.nlme) #.257 (marginal) / .257 (conditional)


##----------------------------------------------------------------Plot confidence intervals

##confidence intervals of models considering functional components

Functional_models <- list(RC_3 = Recovery_SPEI3.nlme, RC_12 = Recovery_SPEI12.nlme, 
                          RS_3 = Resistance_SPEI3.nlme, RS_12 = Resistance_SPEI12.nlme, RS_24 = Resistance_SPEI24.nlme)

Fun_mod_CI <- do.call(rbind, lapply(names(Functional_models), function(nm) {
  Mod <- Functional_models[[nm]]
  CI_tab <- as.data.frame(intervals(Mod, which = "fixed")[[1]])
  CI_sf <- CI_tab[grep("Across_PC1_CWM", x = row.names(CI_tab)), ]
  CI_fd <- CI_tab[grep("Yr_Comb_rao_FP", x = row.names(CI_tab)), ]
  CI_res <- data.frame(rbind(CI_sf, CI_fd), Predictor = c(rep("Across_PC1_CWM", nrow(CI_sf)),
                                                          rep("Yr_Comb_rao_FP", nrow(CI_fd))),
                       SPEI_Resp = nm)
  return(CI_res)
  }))

Fun_mod_CI$SPEI_cat <- NA

Fun_mod_CI[grep("ModerateWet", row.names(Fun_mod_CI)), "SPEI_cat"] <- "ModerateWet"
Fun_mod_CI[grep("ExtremeWet", row.names(Fun_mod_CI)), "SPEI_cat"] <- "ExtremeWet"
Fun_mod_CI[grep("ExtremeDry", row.names(Fun_mod_CI)), "SPEI_cat"] <- "ExtremeDry"

Fun_mod_CI[is.na(Fun_mod_CI$SPEI_cat) & grepl("RC_", Fun_mod_CI$SPEI_Resp), "SPEI_cat"] <- "Normal"
Fun_mod_CI[is.na(Fun_mod_CI$SPEI_cat), "SPEI_cat"] <- "ModerateDry"

#attach confidence intervals for recovery model SPEI24, which was fitted using lm
Recovery_SPEI24.CI <- data.frame(confint(Recovery_SPEI24),
                                 est. = unname(coef(Recovery_SPEI24)),
                                 Predictor = names(coef(Recovery_SPEI24)), 
                                 SPEI_Resp = "RC_24",
                                 SPEI_cat = "Normal")

Recovery_SPEI24.CI <- Recovery_SPEI24.CI[Recovery_SPEI24.CI$Predictor %in% c("Across_PC1_CWM", "Yr_Comb_rao_FP"),
                                         c(1, 3, 2, 4, 5, 6)]

colnames(Recovery_SPEI24.CI) <- colnames(Fun_mod_CI)


Fun_mod_CI <- rbind(Fun_mod_CI, Recovery_SPEI24.CI)

##confidence intervals of models considering species richness

Taxonomic_models <- list(RC_3 = Recovery_SPEI3.tx.nlme, RC_12 = Recovery_SPEI12.tx.nlme, 
                          RS_3 = Resistance_SPEI3.tx.nlme, RS_12 = Resistance_SPEI12.tx.nlme, RS_24 = Resistance_SPEI24.tx.nlme)


Tax_mod_CI <- do.call(rbind, lapply(names(Taxonomic_models), function(nm) {
  Mod <- Taxonomic_models[[nm]]
  CI_tab <- as.data.frame(intervals(Mod, which = "fixed")[[1]])
  CI_sr <- CI_tab[grep("Yr_species_rich", x = row.names(CI_tab)), ]
  CI_res <- data.frame(CI_sr, Predictor = "Yr_species_rich",
                       SPEI_Resp = nm)
  return(CI_res)
  }))

Tax_mod_CI$SPEI_cat <- NA

Tax_mod_CI[grep("ModerateWet", row.names(Tax_mod_CI)), "SPEI_cat"] <- "ModerateWet"
Tax_mod_CI[grep("ExtremeWet", row.names(Tax_mod_CI)), "SPEI_cat"] <- "ExtremeWet"
Tax_mod_CI[grep("ExtremeDry", row.names(Tax_mod_CI)), "SPEI_cat"] <- "ExtremeDry"

Tax_mod_CI[is.na(Tax_mod_CI$SPEI_cat) & grepl("RC_", Tax_mod_CI$SPEI_Resp), "SPEI_cat"] <- "Normal"
Tax_mod_CI[is.na(Tax_mod_CI$SPEI_cat), "SPEI_cat"] <- "ModerateDry"

#attach confidence intervals for recovery model SPEI24, which was fitted using lm
Recovery_SPEI24.tx.CI <- data.frame(confint(Recovery_SPEI24.tx),
                                 est. = unname(coef(Recovery_SPEI24.tx)),
                                 Predictor = names(coef(Recovery_SPEI24.tx)), 
                                 SPEI_Resp = "RC_24",
                                 SPEI_cat = "Normal")

Recovery_SPEI24.tx.CI <- Recovery_SPEI24.tx.CI[Recovery_SPEI24.tx.CI$Predictor %in% c("Yr_species_rich"),
                                         c(1, 3, 2, 4, 5, 6)]

colnames(Recovery_SPEI24.tx.CI) <- colnames(Tax_mod_CI)

Tax_mod_CI <- rbind(Tax_mod_CI, Recovery_SPEI24.tx.CI)

#put together CI(s) for recovery and resistance

Recovery_CIs <- rbind(Fun_mod_CI[grep("RC_", x = Fun_mod_CI$SPEI_Resp), ],
                      Tax_mod_CI[grep("RC_", x = Tax_mod_CI$SPEI_Resp), ])


Resistance_CIs <- rbind(Fun_mod_CI[grep("RS_", x = Fun_mod_CI$SPEI_Resp), ],
                      Tax_mod_CI[grep("RS_", x = Tax_mod_CI$SPEI_Resp), ])


#plot recovery models

#order factor for increasing SPEI time-scale
Recovery_CIs$SPEI_Resp <- factor(Recovery_CIs$SPEI_Resp, levels = unique(Recovery_CIs$SPEI_Resp))

#order factor for increasing water availability
Recovery_CIs$SPEI_cat <- factor(Recovery_CIs$SPEI_cat, levels = unique(Recovery_CIs$SPEI_cat))

Recovery_plot <- ggplot(Recovery_CIs, aes(x = SPEI_cat, y = est., col = SPEI_cat)) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  geom_point(cex = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  scale_color_manual(values = Drought_pal[c(3, 4, 5)], name = "Standardised Precip. Evapotransp. Index category") +
  facet_grid(Predictor ~ SPEI_Resp, scales = "free",
             labeller = labeller(SPEI_Resp = c("RC_3" = "SPEI-3",
                                               "RC_12" = "SPEI-12",
                                               "RC_24" = "SPEI-24"), 
                                 Predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                          "Yr_Comb_rao_FP" = "Functional diversity",
                                          "Yr_species_rich" = "Species richness"))) +
  ggtitle("Recovery") + ylab("LogR, model coefficient (95% CI)") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "bottom", axis.text.x.bottom = element_blank(),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        title = element_text(size = 16),
        panel.background = element_rect(fill = NA, color = "black"))

ggsave(plot = Recovery_plot, "~/Documents/ClimateExtremes/CleanCode/Recovery_CI.jpeg", device = "jpeg", dpi = 300, width = 22, height = 20, units = "cm")

#plot resistance models

#order factor for increasing SPEI time-scale
Resistance_CIs$SPEI_Resp <- factor(Resistance_CIs$SPEI_Resp, levels = unique(Resistance_CIs$SPEI_Resp))

#order factor for increasing water availability
Resistance_CIs$SPEI_cat <- factor(Resistance_CIs$SPEI_cat, levels = rev(unique(Resistance_CIs$SPEI_cat)))

Resistance_plot <- ggplot(Resistance_CIs, aes(x = SPEI_cat, y = est., col = SPEI_cat)) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  geom_point(cex = 2.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  scale_color_manual(values = Drought_pal[c(1, 2)], name = "Standardised Precip. Evapotransp. Index category") +
  facet_grid(Predictor ~ SPEI_Resp, scales = "free",
             labeller = labeller(SPEI_Resp = c("RS_3" = "SPEI-3",
                                               "RS_12" = "SPEI-12",
                                               "RS_24" = "SPEI-24"), 
                                 Predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness"))) +
  ggtitle("Resistance") + ylab("LogR ref-plot, model coefficient (95% CI)") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "bottom", axis.text.x.bottom = element_blank(),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        title = element_text(size = 16),
        panel.background = element_rect(fill = NA, color = "black"))

ggsave(plot = Resistance_plot, "~/Documents/ClimateExtremes/CleanCode/Resistance_CI.jpeg", device = "jpeg", dpi = 300, width = 22, height = 20, units = "cm")

##----------------------------------------------------------------Plot partial effects

#extract fitted values for functional models
Fun_mod_fitted <- do.call(rbind, lapply(names(Functional_models), function(nm) {
  Slow_fast.fit <- as.data.frame(predictorEffect(predictor = "Across_PC1_CWM",
                                                 mod = Functional_models[[nm]]))
  Slow_fast.fit <- data.frame(Slow_fast.fit, Predictor = "Across_PC1_CWM", SPEI_resp = nm)
  FD.fit <- as.data.frame(predictorEffect(predictor = "Yr_Comb_rao_FP",
                                          mod = Functional_models[[nm]]))
  FD.fit <- data.frame(FD.fit, Predictor = "Yr_Comb_rao_FP", SPEI_resp = nm)
  colnames(Slow_fast.fit)[1] <- colnames(FD.fit)[1] <- "X"
  Df <- rbind(Slow_fast.fit, FD.fit)
  return(Df)
  }))

#attach fitted values for recovery model SPEI-24
Recovery_SPEI24.fit <- as.data.frame(predictorEffects(mod = Recovery_SPEI24))[c("Across_PC1_CWM", "Yr_Comb_rao_FP")]

Recovery_SPEI24.fit <- do.call(rbind, lapply(names(Recovery_SPEI24.fit), function(nm) {
  Df <- Recovery_SPEI24.fit[[nm]]
  Df <- data.frame(Df, Predictor = nm, SPEI_resp = "RC_24")
  colnames(Df)[1] <- "X"
  Df$SPEI_final_cat <- "Normal"
  Df <- Df[c(1, 8, 2, 3, 4, 5, 6, 7)]
  return(Df)
  }))

Fun_mod_FIT <- rbind(Fun_mod_fitted, Recovery_SPEI24.fit)

#extract fitted values for taxonomic models
Tax_mod_fitted <- do.call(rbind, lapply(names(Taxonomic_models), function(nm) {
  Sr.fit <- as.data.frame(predictorEffect(predictor = "Yr_species_rich", mod = Taxonomic_models[[nm]]))
  Sr.fit <- data.frame(Sr.fit, Predictor = "Yr_species_rich", SPEI_resp = nm)
  colnames(Sr.fit)[1] <- "X"
  return(Sr.fit)
}))

#attach fitted values for recovery model SPEI-24
Recovery_SPEI24.tx.fit <- as.data.frame(predictorEffect(predictor = "Yr_species_rich", mod = Recovery_SPEI24.tx))

Recovery_SPEI24.tx.fit <- data.frame(Recovery_SPEI24.tx.fit, Predictor = "Yr_species_rich",
                                     SPEI_resp = "RC_24", SPEI_final_cat = "Normal")

colnames(Recovery_SPEI24.tx.fit)[1] <- "X"

Recovery_SPEI24.tx.fit <- Recovery_SPEI24.tx.fit[c(1, 8, 2:7)]

Tax_mod_FIT <- rbind(Tax_mod_fitted, Recovery_SPEI24.tx.fit)

#put together data.frame(s) of fitted values for recovery and resistance

Recovery_FIT <- rbind(Fun_mod_FIT[grep("RC_", x = Fun_mod_FIT$SPEI_resp), ],
                      Tax_mod_FIT[grep("RC_", x = Tax_mod_FIT$SPEI_resp), ])


Resistance_FIT <- rbind(Fun_mod_FIT[grep("RS_", x = Fun_mod_FIT$SPEI_resp), ],
                        Tax_mod_FIT[grep("RS_", x = Tax_mod_FIT$SPEI_resp), ])



#plot recovery models

#order factor for increasing SPEI time-scale
Recovery_FIT$SPEI_resp <- factor(Recovery_FIT$SPEI_resp, levels = unique(Recovery_FIT$SPEI_resp))

#order factor for increasing water availability
Recovery_FIT$SPEI_final_cat <- factor(Recovery_FIT$SPEI_final_cat, levels = as.character(unique(Recovery_FIT$SPEI_final_cat)))

Recovery_partial_plot <- ggplot(Recovery_FIT, aes(x = X, y = fit, fill = SPEI_final_cat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  geom_line(aes(col = SPEI_final_cat), lwd = 1) +
  scale_color_manual(values = Drought_pal[c(3, 4, 5)], name = "SPEI category") +
  scale_fill_manual(values = Drought_pal[c(3, 4, 5)], name = "SPEI category") +
  facet_grid(SPEI_resp ~ Predictor, scales = "free",
             labeller = labeller(SPEI_resp = c("RC_3" = "SPEI-3",
                                               "RC_12" = "SPEI-12",
                                               "RC_24" = "SPEI-24"), 
                                 Predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness"))) +
  ggtitle("Recovery") + ylab("LogR") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "bottom", axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        title = element_text(size = 16))


ggsave(plot = Recovery_partial_plot, "~/Documents/ClimateExtremes/CleanCode/Recovery_partial.jpeg", device = "jpeg", dpi = 300, width = 24, height = 22, units = "cm")

#plot resistance models

#order factor for increasing SPEI time-scale
Resistance_FIT$SPEI_resp <- factor(Resistance_FIT$SPEI_resp, levels = unique(Resistance_FIT$SPEI_resp))

#order factor for increasing water availability
Resistance_FIT$SPEI_final_cat <- factor(Resistance_FIT$SPEI_final_cat, levels = as.character(unique(Resistance_FIT$SPEI_final_cat)))

Resistance_partial_plot <- ggplot(Resistance_FIT, aes(x = X, y = fit, fill = SPEI_final_cat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .2) +
  geom_line(aes(col = SPEI_final_cat), lwd = 1) +
  scale_color_manual(values = Drought_pal[c(1, 2)], name = "SPEI category") +
  scale_fill_manual(values = Drought_pal[c(1, 2)], name = "SPEI category") +
  facet_grid(SPEI_resp ~ Predictor, scales = "free",
             labeller = labeller(SPEI_resp = c("RS_3" = "SPEI-3",
                                               "RS_12" = "SPEI-12",
                                               "RS_24" = "SPEI-24"), 
                                 Predictor = c("Across_PC1_CWM" = "Slow-fast continuum",
                                               "Yr_Comb_rao_FP" = "Functional diversity",
                                               "Yr_species_rich" = "Species richness"))) +
  ggtitle("Resistance") + ylab("LogR ref-plot") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "bottom", axis.text.x.bottom = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        title = element_text(size = 16))

ggsave(plot = Resistance_partial_plot, "~/Documents/ClimateExtremes/CleanCode/Resistance_partial.jpeg", device = "jpeg", dpi = 300, width = 24, height = 22, units = "cm")
