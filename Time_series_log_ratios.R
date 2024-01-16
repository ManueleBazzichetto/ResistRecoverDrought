##Code to replicate results presented in the manuscript "Biodiversity promotes resistance but dominant species shape recovery of grasslands under extreme drought" by Bazzichetto et al.

##Time-series analysis of annual (LogR) and plot reference (LogRref-plot) log ratio

library(ggplot2)
library(ggpubr)
#library(lme4)
library(nlme)
library(car)
library(effects)

##----------------------------------Create datasets to model LogR and LogRref-plot


#info on each variable's meaning (variables selected below):
#LogR: annual log ratio
#LogR_tr_plot: plot reference log ratio
#Explo: 3-levels categorical variable, each level is a region (represented by the acronym of the region name)
#Yr_lui: (yearly values of) land-use intensity
#Across_PC1_CWM: slow-fast continuum. computed performing a PCA on (yearly values of) CWMs
#Yr_Comb_rao_FP: (yearly values of) functional diversity
#Yr_species_rich: (yearly values of) species richness
#Year: sampling year
#Useful_EP_PlotID: ID of the vegetaton plots
#day_of_year: day of biomass harvest
#SoilMoist10: soil moisture at first 10 cm depth

#subset the 'biomass' data.frame to only keep predictors used in the analyses
TA_biom.LogRatios <- biomass[!(biomass$Year %in% c("2020")),
                             c("LogR", "LogR_tr_plot", "Explo",
                               "Yr_lui", "Across_PC1_CWM",
                               "Yr_Comb_rao_FP",
                               "Yr_species_rich",
                               "Year", "Useful_EP_PlotID", "day_of_year", "SoilMoist10")]

#check NAs and remove them in a way that maximizes sample size to model LogR and LogRref-plot 

#check NAs for LogR
sum(is.na(TA_biom.LogRatios$LogR)); which(is.na(TA_biom.LogRatios$LogR) & !TA_biom.LogRatios$Year %in% c("2009"))

#check NAs for LogR_tr_plot
sum(is.na(TA_biom.LogRatios$LogR_tr_plot)); which(is.na(TA_biom.LogRatios$LogR_tr_plot))

#check NAs for predictors
lapply(TA_biom.LogRatios[c(4, 5, 6, 7, 10, 11)], function(i) which(is.na(i)))

#land-use intensity, slow-fast continuum, functional diversity, species richness have all same NAs
#day_of_year has 3 NAs in unique positions

#exclude NAs for land-use intensity (Yr_lui). Note that this will exclude NAs for all other predictors, but day_of_year - NAs are all at the same position
TA_biom.LogRatios <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$Yr_lui), ]

#exclude NAs for day_of_year -> these are only 3
TA_biom.LogRatios <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$day_of_year), ]

#check correlation species richness vs. slow-fast continuum, functional diversity and land-use intensity
Cor_sr_vs_all <- do.call(rbind, lapply(sort(unique(TA_biom.LogRatios$Year)), function(yr) {
  df_for_cor <- TA_biom.LogRatios[TA_biom.LogRatios$Year == yr,
                                  c("Yr_lui", "Across_PC1_CWM", "Yr_Comb_rao_FP", "Yr_species_rich")]
  sr_vs_all <- list(cor_sr_lui = cor(df_for_cor$Yr_species_rich, df_for_cor$Yr_lui),
                 cor_sr_sf = cor(df_for_cor$Yr_species_rich, df_for_cor$Across_PC1_CWM),
                 cor_sr_fd = cor(df_for_cor$Yr_species_rich, df_for_cor$Yr_Comb_rao_FP))
  df_res <- as.data.frame(c(sr_vs_all, year = yr))
  return(df_res)
}))

round(colMeans(Cor_sr_vs_all[c(1, 2, 3)]), 2)

#generate specific datasets for each of the two log ratios (LogR and LogRref-plot)

#dataset for LogR (2009 is excluded as there is no biomass data available for 2008, so all values of LogR for 2009 are NA)
TA_biom.LogR <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$LogR), !colnames(TA_biom.LogRatios) %in% c("LogR_tr_plot")]

#re-order dataset so that vegetation plots are alphabetically ordered (by plot name), and also chronologically ordered (by year)
TA_biom.LogR <- do.call(rbind, lapply(sort(unique(TA_biom.LogR$Useful_EP_PlotID)), function(nm) {
  df.plot <- TA_biom.LogR[TA_biom.LogR$Useful_EP_PlotID == nm, ]
  df.plot <- df.plot[order(df.plot$Year), ]
  return(df.plot)
})
)

#double-check NAs
anyNA(TA_biom.LogR)
unique(TA_biom.LogR$Year)
unique(TA_biom.LogR$Useful_EP_PlotID)

#the same procedure is done for the dataset of the LogR_tr_plot
TA_biom.LogR_tr_plot <- TA_biom.LogRatios[!is.na(TA_biom.LogRatios$LogR_tr_plot), !colnames(TA_biom.LogRatios) %in% c("LogR")]

#re-order dataset so that vegetation plots are alphabetically ordered (by plot name), and also chronologically ordered (by year)
TA_biom.LogR_tr_plot <- do.call(rbind, lapply(sort(unique(TA_biom.LogR_tr_plot$Useful_EP_PlotID)), function(nm) {
  Df <- TA_biom.LogR_tr_plot[TA_biom.LogR_tr_plot$Useful_EP_PlotID == nm, ]
  Df <- Df[order(Df$Year), ]
  return(Df)
})
)

#double-check NAs
anyNA(TA_biom.LogR_tr_plot)
unique(TA_biom.LogR_tr_plot$Year)
unique(TA_biom.LogR_tr_plot$Useful_EP_PlotID)

#create labels for predictors used in time-series models for LogR and LogRref-plot
Lab_vars.f <- c('Across_PC1_CWM' = "Slow-fast continuum",
                'Yr_Comb_rao_FP' = "Functional diversity",
                'Yr_lui' = "Land-use intensity",
                'day_of_year' = "Day of year", 
                "SoilMoist10" = "Soil moisture")

Lab_vars.sr <- c('Yr_species_rich' = "Species richness",
                 'Yr_lui' = "Land-use intensity",
                 'day_of_year' = "Day of year", 
                 "SoilMoist10" = "Soil moisture")

##------------------------------------------------------Time-series analysis

#The aim of this time-series analysis is to model the LogR and LogRref-plot as a function of
#functional and taxonomic diversity, slow-fast continuum, land-use intensity, day of biomass harvest (i.e., day of the year),
#and soil moisture. Note that the effect of taxonomic diversity is modelled separataley from that
#of functional diversity and slow-fast continuum
#To account for the non-independence among observations sampled from the same vegetation plot over time, 
#I will fit different temporal autocorrelation models and select the best fitting one
#(in terms of parsimony and actual reduction in temporal autocorrelation estimated from normalised residuals)
#Also, although the repeated measurements design of the dataset calls for a plot-specific random intercept
#both the LogR and LogRref-plot are measures 'standardized' per plot, and therefore would inevitably lead
#to zero variance of an hypothetical random intercept included in the model

#-------------------------------Analysis LogR 

#check NA
anyNA(TA_biom.LogR)
#check str of TA_biom.LogR
str(TA_biom.LogR)

#save predictors' averages in an object as these will be needed later to transform model predictions
#back to the original scale (models are fitted on centered predictors)
TA_LogR_sc_fact <- attr(scale(TA_biom.LogR[!colnames(TA_biom.LogR) %in% c("LogR", "Explo", "Year", "Useful_EP_PlotID")], center = T, scale = F), "scaled:center")

#center all numerical predictors
TA_biom.LogR.cent <- data.frame(TA_biom.LogR[c("LogR", "Explo", "Year", "Useful_EP_PlotID")],
                                scale(TA_biom.LogR[!colnames(TA_biom.LogR) %in% c("LogR", "Explo", "Year", "Useful_EP_PlotID")], center = T, scale = F))

#check averages are all 0
summary(TA_biom.LogR.cent)

#check years in the data.frame
unique(TA_biom.LogR.cent$Year)

#add a column 'Year_int', which is the integer version of Year. 
#an integer time variable is needed by nlme for the temporal autocorrelation structures
TA_biom.LogR.cent$Year_integer <- as.integer(TA_biom.LogR.cent$Year)


#-------------------------------Functional components 

#fit models for the functional components (functional diversity and slow-fast continuum)

#start with a linear model
#note that dropping the intercept and excluding main effect of predictors in interaction with year allows fitting year-specific
#relationships between the predictors and the response
#also, this avoids hypothesis testing being based on the reference level for year
TS_modLM_logr.fc <- gls(LogR ~ - 1 + Year + Explo + Year:(Across_PC1_CWM + Yr_Comb_rao_FP + Yr_lui + day_of_year + SoilMoist10), 
                         data = TA_biom.LogR.cent)


#check temporal autocorrelation in residuals - if significant temporal autocorrelation is detected,
#then this has to be accounted for by including a temporal autocorrelation structure
#there is indeed evidence of temporal autocorrelation
plot(ACF(TS_modLM_logr.fc, form = ~ Year_integer | Useful_EP_PlotID), alpha = 0.05)

Correct_acf_gaps(dtf = data.frame(Resi = resid(TS_modLM_logr.fc, type = "r"),
                                  Year = TA_biom.LogR.cent$Year_integer,
                                  PlotID = TA_biom.LogR.cent$Useful_EP_PlotID), lag = 0:5)


#fit ar(1) - first-order autoregressive model
TS_modAR1_logr.fc <- update(TS_modLM_logr.fc, . ~ .,
                            correlation = corARMA(p = 1, q = 0, form = ~ Year_integer | Useful_EP_PlotID))


#fit ar(2) - second-order autoregressive model
TS_modAR2_logr.fc <- update(TS_modLM_logr.fc, . ~ .,
                            correlation = corARMA(p = 2, q = 0, form = ~ Year_integer | Useful_EP_PlotID)) 


#fit arma(1, 1) - 1 parameter autoregressive model + 1 parameter moving average model
TS_modARMA11_logr.fc <- update(TS_modLM_logr.fc, . ~ .,
                               correlation = corARMA(p = 1, q = 1, form = ~ Year_integer | Useful_EP_PlotID)) 


#check how the different autocorrelation function modelled temporal autocorrelation
plot(ACF(TS_modAR1_logr.fc, resType = "n"), alpha = 0.05)
plot(ACF(TS_modAR2_logr.fc, resType = "n"), alpha = 0.05)
plot(ACF(TS_modARMA11_logr.fc, resType = "n"), alpha = 0.05)

#compare models including different autocorrelation structures using AIC (note that arma(1,1) is not nested in ar and vice-versa)
anova(TS_modAR1_logr.fc, TS_modAR2_logr.fc, TS_modARMA11_logr.fc) #the arma(1,1) seems to provide the best fit

#diagnostics for normality
qqnorm(residuals(TS_modARMA11_logr.fc, type = "n"))

#r2
round(rr2::R2.lik(TS_modARMA11_logr.fc), 2) #0.44

#compute model predictions
TS_fc.fitted <- as.data.frame(predictorEffects(mod = TS_modARMA11_logr.fc)) 

#check names of resulting list of data.frames including predictions for each predictor
names(TS_fc.fitted)

#subset the TS_fc.fitted list to only keep predictions for variables of interest
TS_fc.fitted <- TS_fc.fitted[c(3, 4, 5, 6, 7)]

#check all predictors of interest have been included
all(names(TS_fc.fitted) %in% names(TA_LogR_sc_fact))

#extract p-values for year-specific relationships to highlight years with significant relationships in the prediction plots

#extract tTable from the nlme model object
LogR_fc_tTab <- summary(TS_modARMA11_logr.fc)$tTable 

#this function extracts years with significant relationships
Extract_sign_yrs <- function(table_p, nm, start = 5, stop = 8) {
  sub_tab <- table_p[grepl(pattern = nm, x = rownames(table_p)), ]
  sign_yrs <- rownames(sub_tab[unname(which(sub_tab[, "p-value"] < 0.05)), , drop = F])
  sign_yrs <- substr(sign_yrs, start, stop)
  return(sign_yrs)
}

#apply function to all predictors of interest
LogR_fc_pval <- lapply(names(TS_fc.fitted), Extract_sign_yrs, table_p = LogR_fc_tTab)

#assign names to resulting list - names are always predictors' labels
names(LogR_fc_pval) <- names(TS_fc.fitted)

#walk through list of predictions and attach new columns:
#one will provide the values of the predictors back to the original scale (summing up the average)
#the other will simply assign a 'YES' to years with significant relationships (or a 'NO' otherwise)
TS_fc.fitted <- lapply(TS_fc.fitted, function(x) {
  #transform predictors back to original scale
  x$X_original <- (x[[1]] + unname(TA_LogR_sc_fact[[colnames(x)[1]]]))
  #get years with sign relationships
  sign_yrs <- LogR_fc_pval[[(colnames(x)[1])]]
  if(isTRUE(length(sign_yrs) == 0)) {
    x$sign <- "NO"
  } else {
    x$sign <- "NO"
    x[x[["Year"]] %in% sign_yrs, "sign"] <- "YES"
  }
  return(x)
})

#finally, create a list of ggplot objects (including prediction plots)
TS_fc.fitted.plots <- lapply(names(TS_fc.fitted), function(nm) {
  ggplot(TS_fc.fitted[[nm]], aes(x = X_original, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = sign), alpha = 0.4) +
    geom_line(aes(col = sign)) +
    scale_color_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    scale_fill_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    facet_wrap(~ Year, nrow = 2, ncol = 5) +
    ylab("LogR") + xlab(Lab_vars.f[[nm]]) +
    theme_pubr() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12), strip.text = element_text(size = 12))
}) 

names(TS_fc.fitted.plots) <- names(TS_fc.fitted)

#have a look at how predictions look like
ggarrange(TS_fc.fitted.plots[[1]], TS_fc.fitted.plots[[2]],
          TS_fc.fitted.plots[[3]], TS_fc.fitted.plots[[4]],
          TS_fc.fitted.plots[[5]])


#-------------------------------Taxonomic diversity 

#fit models for taxonomic diversity (species richness)

#note that the code below is not thoroughly commented as the procedure is exactly the same as the one implemented
#for the functional components. So, instructions on each step presented below can be found in the section above

TS_modLM_logr.tc <- gls(LogR ~ -1 + Year + Explo + Year:(Yr_species_rich + Yr_lui + day_of_year + SoilMoist10),
                         data = TA_biom.LogR.cent)

#check temporal autocorrelation in residuals
plot(ACF(TS_modLM_logr.tc, form = ~ Year_integer | Useful_EP_PlotID), alpha = 0.05)

Correct_acf_gaps(dtf = data.frame(Resi = resid(TS_modLM_logr.tc, type = "r"),
                                  Year = TA_biom.LogR.cent$Year_integer,
                                  PlotID = TA_biom.LogR.cent$Useful_EP_PlotID), lag = 0:5)

#fit ar(1)
TS_modAR1_logr.tc <- update(TS_modLM_logr.tc, . ~ .,
                            correlation = corARMA(p = 1, q = 0, form = ~ Year_integer | Useful_EP_PlotID)) 

#fit ar(2)
TS_modAR2_logr.tc <- update(TS_modLM_logr.tc, . ~ .,
                            correlation = corARMA(p = 2, q = 0, form = ~ Year_integer | Useful_EP_PlotID)) 

#fit arma(1, 1)
TS_modARMA11_logr.tc <- update(TS_modLM_logr.tc, . ~ .,
                               correlation = corARMA(p = 1, q = 1, form = ~ Year_integer | Useful_EP_PlotID)) 

#check temporal autocorrelation
plot(ACF(TS_modAR1_logr.tc, resType = "n"), alpha = 0.05)
plot(ACF(TS_modAR2_logr.tc, resType = "n"), alpha = 0.05)
plot(ACF(TS_modARMA11_logr.tc, resType = "n"), alpha = 0.05)

#compare models using AIC
anova(TS_modAR1_logr.tc, TS_modAR2_logr.tc, TS_modARMA11_logr.tc) #ARMA seems to provide the best fit

#diagnostics
qqnorm(residuals(TS_modARMA11_logr.tc, type = "n"))

#r2
round(rr2::R2.lik(TS_modARMA11_logr.tc), 2) #0.41

#predictions 
TS_tc.fitted <- as.data.frame(predictorEffects(mod = TS_modARMA11_logr.tc)) 

names(TS_tc.fitted)

#subset TS_tc.fitted to only keep predictions of interest
TS_tc.fitted <- TS_tc.fitted[c(3, 4, 5, 6)]

all(names(TS_tc.fitted) %in% names(TA_LogR_sc_fact)) #T

#extract years with significant relationships

LogR_tc_tTab <- summary(TS_modARMA11_logr.tc)$tTable

LogR_tc_pval <- lapply(names(TS_tc.fitted), Extract_sign_yrs, table_p = LogR_tc_tTab)

names(LogR_tc_pval) <- names(TS_tc.fitted)

TS_tc.fitted <- lapply(TS_tc.fitted, function(x) {
  x$X_original <- (x[[1]] + TA_LogR_sc_fact[[colnames(x)[1]]])
  sign_yrs <- LogR_tc_pval[[(colnames(x)[1])]]
  if(isTRUE(length(sign_yrs) == 0)) {
    x$sign <- "NO"
  } else {
    x$sign <- "NO"
    x[x[["Year"]] %in% sign_yrs, "sign"] <- "YES"
  }
  return(x)
})

TS_tc.fitted.plots <- lapply(names(TS_tc.fitted), function(nm) {
  ggplot(TS_tc.fitted[[nm]], aes(x = X_original, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = sign), alpha = 0.4) +
    geom_line(aes(col = sign)) +
    scale_color_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    scale_fill_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    facet_wrap(~ Year, nrow = 2, ncol = 5) +
    ylab("LogR") + xlab(Lab_vars.sr[[nm]]) +
    theme_pubr() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12), strip.text = element_text(size = 12))
}) 

names(TS_tc.fitted.plots) <- names(TS_tc.fitted)

ggarrange(TS_tc.fitted.plots[[1]], TS_tc.fitted.plots[[2]],
          TS_tc.fitted.plots[[3]], TS_tc.fitted.plots[[4]])


#-----------Put together prediction plots for functional and taxonomic components, and for land-use intensity, day of the year and soil moisture

#prediction plots for models including functional components (fc) and the taxonomic component (tc)
#add + scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) to prediction plot for species richness in case of nrow = 1 in ggarrange
TS_LogR_partial <- ggarrange(TS_fc.fitted.plots$Across_PC1_CWM, 
                             TS_fc.fitted.plots$Yr_Comb_rao_FP,
                             TS_tc.fitted.plots$Yr_species_rich,
                             nrow = 3, ncol = 1, legend = "none", labels = "auto")

#save plot
ggsave(TS_LogR_partial, filename = "~/Documents/ClimateExtremes/CleanCode/New_figures/TS_LogR_partial.jpeg", device = "jpeg",
       units = "cm", width = 30, height = 45, dpi = 300)

#put together prediction plots for land-use intensity and day of year

#combine lists of predictions from fc and tc model for land-use intensity and day of year
Fitted_fc_lui_dy_sm <- c(TS_fc.fitted[c(3, 4, 5)], TS_tc.fitted[c(2, 3, 4)])

#check nrow for each data.frame of fitted values
sapply(Fitted_fc_lui_dy_sm, nrow) #all 500

#walk through list of predictions to create a new column with the name of the predictor and rename the first column as 'X'
Fitted_fc_lui_dy_sm <- do.call(rbind, lapply(seq_along(Fitted_fc_lui_dy_sm), function(i) {
  dtf <- Fitted_fc_lui_dy_sm[[i]]
  dtf$var <- colnames(dtf)[1]
  colnames(dtf)[1] <- "X"
  return(dtf)
}))

#create a new column and assign to it a unique value, i.e. 'NO'. This is needed to then combine information about 
#significant relationships for land-use intensity, day of year and soil moisture for fc and tc models
Fitted_fc_lui_dy_sm$biod_sign <- "NO"

#the next 2 lines simply re-map information on significance of the relationships to a more complex coding
#which will allow giving different colours to significant relationships observed in fc or tc models
Fitted_fc_lui_dy_sm[1:1500, "biod_sign"] <- ifelse(Fitted_fc_lui_dy_sm$sign[1:1500] == "NO", "FUN_NO", "FUN_YES")

Fitted_fc_lui_dy_sm[1501:3000, "biod_sign"] <- ifelse(Fitted_fc_lui_dy_sm$sign[1501:3000] == "NO", "TAX_NO", "TAX_YES")

#create prediction plot for land-use intensity
TS_LogR_partial_lui <- ggplot(Fitted_fc_lui_dy_sm[Fitted_fc_lui_dy_sm$var == "Yr_lui", ], aes(x = X_original, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = biod_sign), alpha = 0.2) +
  geom_line(aes(col = biod_sign)) +
  scale_color_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                                "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                     labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                                "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  scale_fill_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                               "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                    labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                               "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  facet_wrap(~ Year, nrow = 2, ncol = 5) +
  ylab("LogR") + xlab(Lab_vars.f[["Yr_lui"]]) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

#create prediction plot for day of the year
TS_LogR_partial_day <- ggplot(Fitted_fc_lui_dy_sm[Fitted_fc_lui_dy_sm$var == "day_of_year", ], aes(x = X_original, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = biod_sign), alpha = 0.2) +
  geom_line(aes(col = biod_sign)) +
  scale_color_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                                "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                     labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                                "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  scale_fill_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                               "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                    labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                               "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  facet_wrap(~ Year, nrow = 2, ncol = 5) +
  ylab("LogR") + xlab(Lab_vars.f[["day_of_year"]]) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

#create prediction plot for soil moisture
TS_LogR_partial_sm <- ggplot(Fitted_fc_lui_dy_sm[Fitted_fc_lui_dy_sm$var == "SoilMoist10", ], aes(x = X_original, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = biod_sign), alpha = 0.2) +
  geom_line(aes(col = biod_sign)) +
  scale_color_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                                "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                     labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                                "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  scale_fill_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                               "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                    labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                               "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  facet_wrap(~ Year, nrow = 2, ncol = 5) +
  ylab("LogR") + xlab(Lab_vars.f[["SoilMoist10"]]) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))


#combine plots
TS_LogR_part_luidaysm <- ggarrange(TS_LogR_partial_lui, TS_LogR_partial_day, TS_LogR_partial_sm,
                                 nrow = 3, ncol = 1, legend = "none", labels = "auto")

#save plot
ggsave(TS_LogR_part_luidaysm, filename = "~/Documents/ClimateExtremes/CleanCode/New_figures/TS_LogR_part_luidaysm.jpeg", device = "jpeg",
       units = "cm", width = 30, height = 45, dpi = 300, bg = "white")



#-------------------------------Analysis LogRref-plot

#check NA
anyNA(TA_biom.LogR_tr_plot)
#check str of TA_biom.LogR_tr_plot
str(TA_biom.LogR_tr_plot)

#save predictors' averages in an object as these will be needed later to transform model predictions
#back to the original scale (models are fitted on centered predictors)
TA_LogRtrpl_sc_fact <- attr(scale(TA_biom.LogR_tr_plot[!colnames(TA_biom.LogR_tr_plot) %in% c("LogR_tr_plot", "Explo", "Year", "Useful_EP_PlotID")], center = T, scale = F), "scaled:center")

#center all numerical predictors
TA_biom.LogR_tr_plot.cent <- data.frame(TA_biom.LogR_tr_plot[c("LogR_tr_plot", "Explo", "Year", "Useful_EP_PlotID")],
                                        scale(TA_biom.LogR_tr_plot[!colnames(TA_biom.LogR_tr_plot) %in% c("LogR_tr_plot", "Explo", "Year", "Useful_EP_PlotID")], center = T, scale = F))

#check all averages are 0
summary(TA_biom.LogR_tr_plot.cent)

#check years in the data.frame
unique(TA_biom.LogR_tr_plot.cent$Year)

#add Year_int (integer version of Year needed by nlme for the temporal autocorrelation structure)
TA_biom.LogR_tr_plot.cent$Year_integer <- as.integer(TA_biom.LogR_tr_plot.cent$Year)

#-------------------------------Functional components 

#fit models for the functional components (functional diversity and slow-fast continuum)

#note that the code below is not thoroughly commented as the procedure is exactly the same as the one implemented
#for the model for LogR including the functional components. So, instructions on each step presented below can be found in that section

#start with a linear model
TS_modLM_logr_plot.fc <- gls(LogR_tr_plot ~ -1 + Year + Explo + Year:(Across_PC1_CWM + Yr_Comb_rao_FP + Yr_lui + day_of_year + SoilMoist10),
                              data = TA_biom.LogR_tr_plot.cent)


#check temporal autocorrelation in residuals
plot(ACF(TS_modLM_logr_plot.fc, form = ~ Year_integer | Useful_EP_PlotID), alpha = 0.05)

Correct_acf_gaps(dtf = data.frame(Resi = resid(TS_modLM_logr_plot.fc, type = "r"),
                                  Year = TA_biom.LogR_tr_plot.cent$Year_integer,
                                  PlotID = TA_biom.LogR_tr_plot.cent$Useful_EP_PlotID),
                 lag = 0:5)

#fit ar(1)
TS_modAR1_logr_plot.fc <- update(TS_modLM_logr_plot.fc, . ~ .,
                                 correlation = corARMA(p = 1, q = 0, form = ~ Year_integer | Useful_EP_PlotID)) 


#fit ar(1) - note that the ar(1) model seems to provide a good fit to the temporal autocor structure
plot(ACF(TS_modAR1_logr_plot.fc, resType = "n"), alpha = 0.05)

#diagnostics
qqnorm(residuals(TS_modAR1_logr_plot.fc, type = "n"))

#r2
round(rr2::R2.lik(TS_modAR1_logr_plot.fc), 2) #0.41

#predictions
TS_fc_logr_refp.fitted <- as.data.frame(predictorEffects(TS_modAR1_logr_plot.fc))

#subset TS_fc_logr_refp.fitted to only keep predictions of interest
TS_fc_logr_refp.fitted <- TS_fc_logr_refp.fitted[c(3, 4, 5, 6, 7)]

#extract table p-values
LogRrefp_fc_tTab <- summary(TS_modAR1_logr_plot.fc)$tTable

#get years with statistically significant years
LogRrefp_fc_pval <- lapply(names(TS_fc_logr_refp.fitted), Extract_sign_yrs, table_p = LogRrefp_fc_tTab)

names(LogRrefp_fc_pval) <- names(TS_fc_logr_refp.fitted)

#transform predictors back to original scale and add column with information on years with significant relationships
TS_fc_logr_refp.fitted <- lapply(TS_fc_logr_refp.fitted, function(x) {
  x$X_original <- (x[[1]] + unname(TA_LogRtrpl_sc_fact[[colnames(x)[1]]]))
  sign_yrs <- LogRrefp_fc_pval[[(colnames(x)[1])]]
  if(isTRUE(length(sign_yrs) == 0)) {
    x$sign <- "NO"
  } else {
    x$sign <- "NO"
    x[x[["Year"]] %in% sign_yrs, "sign"] <- "YES"
  }
  return(x)
})

TS_fc_logr_refp.plots <- lapply(names(TS_fc_logr_refp.fitted), function(nm) {
  ggplot(TS_fc_logr_refp.fitted[[nm]], aes(x = X_original, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = sign), alpha = 0.4) +
    geom_line(aes(col = sign)) +
    scale_color_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    scale_fill_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    facet_wrap(~ Year, nrow = 2, ncol = 6) +
    ylab("LogR ref-plot") + xlab(Lab_vars.f[[nm]]) +
    theme_pubr() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12), strip.text = element_text(size = 12))
})


names(TS_fc_logr_refp.plots) <- names(TS_fc_logr_refp.fitted)

ggarrange(TS_fc_logr_refp.plots[[1]], TS_fc_logr_refp.plots[[2]],
          TS_fc_logr_refp.plots[[3]], TS_fc_logr_refp.plots[[4]],
          TS_fc_logr_refp.plots[[5]])


#-------------------------------Taxonomic diversity 

#fit models for taxonomic diversity (species richness)

#note that the code below is not thoroughly commented as the procedure is exactly the same as the one implemented
#for the model for LogR including the functional components. So, instructions on each step presented below can be found in that section

TS_modLM_logr_plot.tc <- gls(LogR_tr_plot ~ -1 + Year + Explo + Year:(Yr_species_rich + Yr_lui + day_of_year + SoilMoist10),
                              data = TA_biom.LogR_tr_plot.cent)


#check temporal autocorrelation in normalised residuals
plot(ACF(TS_modLM_logr_plot.tc, form = ~ Year_integer | Useful_EP_PlotID), alpha = 0.05)

Correct_acf_gaps(dtf = data.frame(Resi = resid(TS_modLM_logr_plot.tc, type = "r"),
                                  Year = TA_biom.LogR_tr_plot.cent$Year_integer,
                                  PlotID = TA_biom.LogR_tr_plot.cent$Useful_EP_PlotID),
                 lag = 0:5)

#fit ar(1) - note that the ar(1) model seems to provide a good fit to the temporal autocor structure
TS_modAR1_logr_plot.tc <- update(TS_modLM_logr_plot.tc, . ~ .,
                                 correlation = corARMA(p = 1, q = 0, form = ~ Year_integer | Useful_EP_PlotID))

plot(ACF(TS_modAR1_logr_plot.tc, resType = "n"), alpha = 0.05)

#diagnostics
qqnorm(residuals(TS_modAR1_logr_plot.tc, type = "n"))

#r2
round(rr2::R2.lik(TS_modAR1_logr_plot.tc), 2) #0.38

#predictions
TS_tc_logr_refp.fitted <- as.data.frame(predictorEffects(TS_modAR1_logr_plot.tc))

#subset to only keep predictions of interest
TS_tc_logr_refp.fitted <- TS_tc_logr_refp.fitted[c(3, 4, 5, 6)]

#extract summary table
LogRrefp_tc_tTab <- summary(TS_modAR1_logr_plot.tc)$tTable

#get years with significant relationships
LogRrefp_tc_pval <- lapply(names(TS_tc_logr_refp.fitted), Extract_sign_yrs, table_p = LogRrefp_tc_tTab)

names(LogRrefp_tc_pval) <- names(TS_tc_logr_refp.fitted)

TS_tc_logr_refp.fitted <- lapply(TS_tc_logr_refp.fitted, function(x) {
  x$X_original <- (x[[1]] + unname(TA_LogRtrpl_sc_fact[[colnames(x)[1]]]))
  sign_yrs <- LogRrefp_tc_pval[[(colnames(x)[1])]]
  if(isTRUE(length(sign_yrs) == 0)) {
    x$sign <- "NO"
  } else {
    x$sign <- "NO"
    x[x[["Year"]] %in% sign_yrs, "sign"] <- "YES"
  }
  return(x)
})

TS_tc_logr_refp.plots <- lapply(names(TS_tc_logr_refp.fitted), function(nm) {
  ggplot(TS_tc_logr_refp.fitted[[nm]], aes(x = X_original, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = sign), alpha = 0.4) +
    geom_line(aes(col = sign)) +
    scale_color_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    scale_fill_manual(values = c("YES" = "#333333", "NO" = "lightgrey"), name = "Significant") +
    facet_wrap(~ Year, nrow = 2, ncol = 6) +
    ylab("LogR ref-plot") + xlab(Lab_vars.sr[[nm]]) +
    theme_pubr() +
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12), strip.text = element_text(size = 12))
})


names(TS_tc_logr_refp.plots) <- names(TS_tc_logr_refp.fitted)

ggarrange(TS_tc_logr_refp.plots[[1]], TS_tc_logr_refp.plots[[2]],
          TS_tc_logr_refp.plots[[3]], TS_tc_logr_refp.plots[[4]])


#-----------Put together prediction plots for functional and taxonomic components, and for land-use intensity, day of the year and soil moisture


#plots for models including fc or tc
#add + scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) to plot for species ricnhess if nrow = 1 in ggarrange
TS_LogRtrpl_partial <- ggarrange(TS_fc_logr_refp.plots$Across_PC1_CWM, 
                                 TS_fc_logr_refp.plots$Yr_Comb_rao_FP,
                                 TS_tc_logr_refp.plots$Yr_species_rich,
                                 nrow = 3, ncol = 1, legend = "none", labels = "auto")

#save plot
ggsave(TS_LogRtrpl_partial, filename = "~/Documents/ClimateExtremes/CleanCode/New_figures/TS_LogRtrpl_partial.jpeg", device = "jpeg",
       units = "cm", width = 30, height = 45, dpi = 300)

#create prediction plots for land-use intensity, day of the year and soil moisture

Fitted_fc_lui_dy_sm.logrtrpl <- c(TS_fc_logr_refp.fitted[c(3, 4, 5)], TS_tc_logr_refp.fitted[c(2, 3, 4)])

sapply(Fitted_fc_lui_dy_sm.logrtrpl, nrow) #all 550

Fitted_fc_lui_dy_sm.logrtrpl <- do.call(rbind, lapply(seq_along(Fitted_fc_lui_dy_sm.logrtrpl), function(i) {
  dtf <- Fitted_fc_lui_dy_sm.logrtrpl[[i]]
  dtf$var <- colnames(dtf)[1]
  colnames(dtf)[1] <- "X"
  return(dtf)
}))

Fitted_fc_lui_dy_sm.logrtrpl$biod_sign <- "NO"

Fitted_fc_lui_dy_sm.logrtrpl$biod_sign[1:1650] <- ifelse(Fitted_fc_lui_dy_sm.logrtrpl$sign[1:1650] == "NO", "FUN_NO", "FUN_YES")

Fitted_fc_lui_dy_sm.logrtrpl$biod_sign[1651:3300] <- ifelse(Fitted_fc_lui_dy_sm.logrtrpl$sign[1651:3300] == "NO", "TAX_NO", "TAX_YES")

TS_LogRtrpl_partial_lui <- ggplot(Fitted_fc_lui_dy_sm.logrtrpl[Fitted_fc_lui_dy_sm.logrtrpl$var == "Yr_lui", ], aes(x = X_original, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = biod_sign), alpha = 0.2) +
  geom_line(aes(col = biod_sign)) +
  scale_color_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                                "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                     labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                                "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  scale_fill_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                               "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                    labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                               "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  facet_wrap(~ Year, nrow = 2, ncol = 6) +
  ylab("LogR ref-plot") + xlab(Lab_vars.f[["Yr_lui"]]) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

TS_LogRtrpl_partial_day <- ggplot(Fitted_fc_lui_dy_sm.logrtrpl[Fitted_fc_lui_dy_sm.logrtrpl$var == "day_of_year", ], aes(x = X_original, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = biod_sign), alpha = 0.2) +
  geom_line(aes(col = biod_sign)) +
  scale_color_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                                "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                     labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                                "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  scale_fill_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                               "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                    labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                               "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  facet_wrap(~ Year, nrow = 2, ncol = 6) +
  ylab("LogR ref-plot") + xlab(Lab_vars.f[["day_of_year"]]) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))


TS_LogRtrpl_partial_sm <- ggplot(Fitted_fc_lui_dy_sm.logrtrpl[Fitted_fc_lui_dy_sm.logrtrpl$var == "SoilMoist10", ], aes(x = X_original, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = biod_sign), alpha = 0.2) +
  geom_line(aes(col = biod_sign)) +
  scale_color_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                                "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                     labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                                "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  scale_fill_manual(values = c("FUN_YES" = "darkgreen", "TAX_YES" = "purple",
                               "FUN_NO" = "lightgrey", "TAX_NO" = "lightgrey"),
                    labels = c("FUN_YES" = "YES", "TAX_YES" = "YES",
                               "FUN_NO" = "NO", "TAX_NO" = "NO"), name = "Significant") +
  facet_wrap(~ Year, nrow = 2, ncol = 6) +
  ylab("LogR ref-plot") + xlab(Lab_vars.f[["SoilMoist10"]]) +
  theme_pubr() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

TS_LogRtrpl_part_luidaysm <- ggarrange(TS_LogRtrpl_partial_lui, TS_LogRtrpl_partial_day, TS_LogRtrpl_partial_sm,
                                       nrow = 3, ncol = 1, legend = "none", labels = "auto")

#save plot
ggsave(TS_LogRtrpl_part_luidaysm, filename = "~/Documents/ClimateExtremes/CleanCode/New_figures/TS_LogRtrpl_part_luidaysm.jpeg", device = "jpeg",
       units = "cm", width = 30, height = 45, dpi = 300, bg = "white")
