##Code to replicate results presented in the manuscript "Biodiversity promotes resistance but dominant species shape recovery of grasslands under extreme drought" by Bazzichetto et al.

##Data preparation, computation of the log response ratios, exploratory analysis on biomass trend and 
##correlation between biodiversity (functional and taxonomic diversity), slow-fast continuum and land use

library(dplyr)
library(ggplot2)
library(ggpubr)

##----------------------------------------------------------------Import biomass data

#import biomass data
biomass <- read.csv(file = "~/Documents/ClimateExtremes/CleanCode/Data/biomass.csv",
                    header = T, sep = ",", dec = ".")


#create column with acronyms of the regions of the Biodiversity Exploratories (BE): ALB (South-West), HAI (Center), SCH (North-East)

#extract first letter of BE regions' name
biomass$Explo <- substr(biomass$Useful_EP_PlotID, start = 1, stop = 1)

#create key to rename regions
Explo_key <- c(H = "HAI", A = "ALB", S = "SCH")

#overwrite Explo column replacing BE regions' first letters with acronyms
biomass$Explo <- unname(Explo_key[biomass$Explo])

##----------------------------------------------------------------Computation of log ratios

#---------------------compute LogR (annual log ratio)

#the 'Compute_LogR' function is in the 'Functions' script
#the function operates at the (vegetation) plot level to compute yearly change in biomass value 
biomass.logratio.df <- Compute_LogR(biomass)

#join biomass and biomass.logratio.df
biomass <- dplyr::inner_join(biomass, biomass.logratio.df, by = c("EP_PlotID" = "PlotID", "Explo" = "Expl", "Year" = "Year"))

#--------------------compute summaries of biomass

#compute different, region-specific measures summarising biomass trend:
#year- and regional-specific average biomass
#year- and regional-specific median biomass
#year- and regional-specific 1st quartile biomass
#year- and regional-specific 3rd quartile biomass
biomass.avg <- tapply(biomass[[c("biomass")]], INDEX = list(biomass$Explo, biomass$Year), FUN = mean, na.rm = T, simplify = T)
biomass.med <- tapply(biomass[[c("biomass")]], INDEX = list(biomass$Explo, biomass$Year), FUN = median, na.rm = T, simplify = T)
biomass.low_q <- tapply(biomass[[c("biomass")]], INDEX = list(biomass$Explo, biomass$Year), FUN = quantile, probs = 0.25, na.rm = T, simplify = T)
biomass.up_q <- tapply(biomass[[c("biomass")]], INDEX = list(biomass$Explo, biomass$Year), FUN = quantile, probs = 0.75, na.rm = T, simplify = T)

#create a data.frame with the measures of biomass trend
biomass.med.q <- data.frame(Explo = character(36), year = character(36), Avg = numeric(36), Med = numeric(36),
                            L_q = numeric(36), U_q = numeric(36))

biomass.med.q$Explo <- rep(rownames(biomass.med), each = ncol(biomass.med))
biomass.med.q$year <- rep(colnames(biomass.med), times = 3)
biomass.med.q$Avg <- as.numeric(t(biomass.avg))
biomass.med.q$Med <- as.numeric(t(biomass.med))
biomass.med.q$L_q <- as.numeric(t(biomass.low_q))
biomass.med.q$U_q <- as.numeric(t(biomass.up_q))

#compute region-specific median of yearly (median) biomass value
biomass$biom_trend <- ifelse(biomass$Explo == "ALB", median(biomass.med.q[biomass.med.q$Explo == "ALB", "Med"], na.rm = T),
                             ifelse(biomass$Explo == "HAI", median(biomass.med.q[biomass.med.q$Explo == "HAI", "Med"], na.rm = T),
                                    median(biomass.med.q[biomass.med.q$Explo == "SCH", "Med"], na.rm = T)))

#coerce year to integer
biomass.med.q$year <- as.integer(biomass.med.q$year)

#----------------------compute LogRref-plot (plot reference log ratio)

#compute plot-specific median value of biomass over the time-series
biomass.med.plot <- tapply(biomass[["biomass"]], biomass[["Useful_EP_PlotID"]], median, na.rm = T)

biomass.med.plot <- data.frame(LogR_plot = unname(biomass.med.plot),
                               Useful_EP_PlotID = names(biomass.med.plot))

#join plot-specific median value of biomass to biomass data.frame
biomass <- dplyr::inner_join(x = biomass, y = biomass.med.plot, by = c("Useful_EP_PlotID"))

#compute LogRref-plot (plot reference LogR)
biomass$LogR_tr_plot <- with(biomass, log(biomass/LogR_plot))

#get rid of columns that will not be used
biomass <- biomass[!(colnames(biomass) %in% c("EP_PlotID"))]


##--------------------------------------------Import and process data on functional traits and land-use intensity


#----------------------------import (yearly) values of community weighted means (CWMs) of functional traits

##FOR ME - only keep data with no _pld - change file name and directory
load("~/Documents/ClimateExtremes/CleanCode/Data/fundata.RData")


#get rid of functional traits data pooled over the time-series, as the focus is on yearly values of the traits (CWM)
rm(cwm_pooled_gr, dcpl_gr_pld, comb_raophylo_gr_pld)

#update column names
colnames(dcpl_gr)[-c(1:3)] <- paste("Yr", colnames(dcpl_gr)[-c(1:3)], sep = "_")
colnames(cwm_gr)[-c(1:4, 14)] <- paste("Yr", colnames(cwm_gr)[-c(1:4, 14)], sep = "_")
colnames(comb_raophylo_gr)[-c(1:2, 7)] <- paste("Yr", "Comb", colnames(comb_raophylo_gr)[-c(1:2, 7)], sep = "_")

Yr_join_div <- dplyr::inner_join(dcpl_gr, comb_raophylo_gr, by = c("Useful_EP_PlotID", "Year"))
Yr_join_compl <- dplyr::inner_join(Yr_join_div, cwm_gr, by = c("Useful_EP_PlotID", "Year"))

#delete unused columns
Yr_join_compl$Exploratory.y <- NULL
Yr_join_compl$Exploratory.x <- NULL
Yr_join_compl$SSD_log <- NULL

#check NAs
anyNA(Yr_join_compl) #F

#join yearly data on functional traits to the 'biomass' data.frame

#coerce 'Year' to character
biomass$Year <- as.character(biomass$Year)

#FOR ME - here I inserted this line as (probably) biomass had no Exploratory column, which was in Df_grasslands (not imported here)
biomass$Exploratory <- ifelse(biomass$Explo == "ALB", "A", 
                              ifelse(biomass$Explo == "HAI", "H", "S"))

biomass <- dplyr::left_join(x = biomass, y = Yr_join_compl, c("Useful_EP_PlotID", "Year", "Exploratory"))

#get rid of columns that will not be used
biomass$Hab <- NULL
biomass$Yr_clonal_0 <- NULL

#reorder vars
biomass <- biomass[c("Useful_EP_PlotID", "Year", "day_of_year", "Explo", "Exploratory", colnames(biomass)[!(colnames(biomass) %in% c("Useful_EP_PlotID", "Year", "day_of_year", "Explo", "Exploratory"))])]

#----------------------------import data on (yearly) plant species richness and land use

load("~/Documents/ClimateExtremes/CleanCode/Data/rich+lui.RData")


#lui_separate is a data frame including the value of the land-use intensity index
#richness_plot_year is also a data frame, with the value of plant species richness recorded in each plot (per each year)

#coerce lui_separate to a data.frame
lui_separate <- as.data.frame(lui_separate)

#some of the plot names do not match between lui_separate and richness_plot_year. The code below fixes the problem.
ID_not_in_PlotID <- unique(richness_plot_year$Useful_EP_PlotID[which(is.na(match(richness_plot_year$Useful_EP_PlotID, lui_separate$PlotID)))])
PlotID_not_in_ID <- unique(lui_separate$PlotID[which(is.na(match(lui_separate$PlotID, richness_plot_year$Useful_EP_PlotID)))])
ID_to_match <- setNames(ID_not_in_PlotID, PlotID_not_in_ID)

lui_separate$Useful_EP_PlotID <- lui_separate$PlotID

for(i in names(ID_to_match)) {
  lui_separate$Useful_EP_PlotID[lui_separate$Useful_EP_PlotID %in% i] <- ID_to_match[[i]] #use [[]] here
}

#join data.frame with data on plant species richness and land use
Rich_lui_j <- dplyr::left_join(richness_plot_year, lui_separate, by = c("Useful_EP_PlotID", "Year"))

#rename columns 
colnames(Rich_lui_j)[3] <- "Yr_species_rich"
colnames(Rich_lui_j)[9] <- "Yr_lui"

#now join the data.frame on plant species richness AND land use with 'biomass'
biomass <- dplyr::left_join(x = biomass,
                            y = Rich_lui_j[c("Useful_EP_PlotID", "Year", "Yr_species_rich",
                                             "Yr_lui")],
                            by = c("Useful_EP_PlotID", "Year"))

biomass$Yr_species_rich <- unname(biomass$Yr_species_rich)

##--------------------------------------------Compute yearly value of slow-fast continuum using CWMs

#the 'Compute_PCA_CWM' function is in the 'Functions' script
#the function performs a Principal Component Analysis on the (yearly) value of the CWMs
#note that this is equivalent to fitting year-specific PCA
Across_yrs_PC1 <- Compute_PCA_CWM(biomass)

#multiply the 'Across_PC1_CWM' columns by -1 to get the axis on a slow -> fast continuum
#low values of 'Across_PC1_CWM' are associated with slow-growing communities and positive values with fast-growing communities 
Across_yrs_PC1$Across_PC1_CWM <- -1*(Across_yrs_PC1$Across_PC1_CWM)

#extract PCA object using the 'Extract_PCA_CWM', which is in the 'Functions' script
PCA_CWM_obj <- Extract_PCA_CWM(biomass) 

#check amount of variance explained by the first axis: ~ 48%
summary(PCA_CWM_obj)

#join Across_yrs_PC1 with 'biomass'
biomass <- dplyr::left_join(x = biomass, y = Across_yrs_PC1, by = c("Useful_EP_PlotID", "Year"))

##------------------------------Import data on soil moisture and join them to the biomass data.frame

#import daily data on soil moisture (10 cm depth) - data already quality checked for physical range and were not interpolated across space

#import data on soil moisture
SoilMoistData <- read.csv(file = "/Users/manbaz/Documents/ClimateExtremes/CleanCode/Data/SM_10_2009_2019_Marta_no_int/plots.csv", sep = ",", dec = ".", stringsAsFactors = F)

#check years
unique(SoilMoistData$year)

#remove NAs for soil moisture
SoilMoistData <- SoilMoistData[!is.na(SoilMoistData$SM_10), ]

#check
sum(is.na(SoilMoistData$SM_10)) #0
range(SoilMoistData$SM_10) #0.18882 64.43250

#compute plot-specific mean value of soil moisture at 10 cm depth
Plot_spec_SM10 <- tapply(SoilMoistData$SM_10, INDEX = list(SoilMoistData$plotID), mean, na.rm = T)


#join soil moisture data to biomass
biomass <- dplyr::left_join(x = biomass, y = data.frame(SoilMoist10 = unname(Plot_spec_SM10),
                                                        PlotID = names(Plot_spec_SM10)),
                            by = c("Useful_EP_PlotID" = "PlotID"))

##--------------------------------------------Have a look at the trend of plant biomass

#create a copy. of biomass.med.q for plotting
biomass.med.q.plot <- biomass.med.q

#add NAs for 2008 to match SPEI plot (which has data from 2008 to 2018)
biomass.med.q.plot <- rbind(biomass.med.q.plot, 
                            data.frame(Explo = c("ALB", "HAI", "SCH"), year = 2008L, Avg = NA, Med = NA, L_q = NA, U_q = NA))


#modified on 29th Nov 2023 to increase text size
Biomass_plot_trend <- ggplot(biomass.med.q.plot[biomass.med.q.plot$year != 2020, ], aes(x = as.factor(year), y = Med)) + 
  geom_hline(data = biomass[biomass$Explo == "SCH" & biomass$Year != "2020", ], aes(yintercept = biom_trend), col = "#BC1EE5", lty = "dotdash") +
  geom_hline(data = biomass[biomass$Explo == "HAI" & biomass$Year != "2020", ], aes(yintercept = biom_trend), col = "#BC1EE5", lty = "dotdash") +
  geom_hline(data = biomass[biomass$Explo == "ALB" & biomass$Year != "2020", ], aes(yintercept = biom_trend), col = "#BC1EE5", lty = "dotdash") +
  geom_point(col = "#139680") +
  geom_errorbar(aes(x = as.factor(year), y = Med, ymin = L_q, ymax = U_q),
                width = .3) +
  facet_wrap(~ Explo, labeller = labeller(Explo = c(ALB = "South-West", HAI = "Center", SCH = "North-East"))) + 
  ylab("Biomass") + xlab(NULL) +
  theme_pubr() +
  theme(legend.position = "none", axis.title.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 11, angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(size = 12))


##--Have a look at the correlation between biomass, species richness, functional and tax. diversity, slow-fast continuum and land use

#these 2 packages are needed to use plant icons in figures
library(cowplot)
library(magick)

#compute average (across years) correlation among species richness, functional and tax. diversity, slow-fast continuum and land use
Cor_Sr_all <- colMeans(do.call(rbind, lapply(unique(biomass$Year)[-12], function(yr) {
  df <- biomass[biomass$Year == yr, c("Yr_species_rich", "Across_PC1_CWM", "Yr_Comb_rao_FP", "Yr_lui")]
  c(Sr_sf = with(df, cor(Yr_species_rich, Across_PC1_CWM, use = "pairwise.complete.obs")),
    Sr_fd = with(df, cor(Yr_species_rich, Yr_Comb_rao_FP, use = "pairwise.complete.obs")),
    Sr_lui = with(df, cor(Yr_species_rich, Yr_lui, use = "pairwise.complete.obs")))
})))

round(Cor_Sr_all, digits = 2)


#create objects to store plant icons (.png)
#Festuca spp.
Icon_fesspp <- "~/Documents/ClimateExtremes/CleanCode/Conceptual_fig/Plants_icon/festuca-spp-fescue.png" #Check attribution at https://ian.umces.edu/media-library/
#Euphorbia esula
Icon_eupesu <- "~/Documents/ClimateExtremes/CleanCode/Conceptual_fig/Plants_icon/euphorbia-esula-leafy-spurge.png" #Check attribution at https://ian.umces.edu/media-library/

#following, a list of plots showing the relationships among all biodiversity-related components + land use
#and box-plots showing the distribution of the values of these variables in the different regions

#biomass vs slow-fast continuum
Biomass_vs_slowfast <- ggplot(biomass[biomass$Year != "2020", ], aes(x = Across_PC1_CWM, y = biomass)) +
  geom_point(aes(col = Explo), alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, col = "black", lwd = 2, lty = "dotdash") +
  geom_smooth(aes(col = Explo), method = "lm", formula = y ~ x, se = F, lwd = 1.2) +
  scale_color_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                     labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Biomass") + xlab("Slow-fast continuum") +
  #annotation_raster(Icon_poapra, xmin = 3, ymin = -40, xmax = 4.5, ymax = 140) +
  #annotation_raster(Icon_censol, xmin = -6, ymin = -60, xmax = -4.5, ymax = 120) +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
        legend.position = "none") 


Biomass_vs_slowfast <- ggdraw(Biomass_vs_slowfast) +
  draw_image(image = Icon_eupesu, x = 1, hjust = 1, valign = 0, width = 0.1, height = 0.13) +
  draw_image(image = Icon_fesspp, x = 0.11, hjust = 0, valign = 0, width = 0.15, height = 0.14)

#box-plots of slow-fast continuum
Slowfast_vs_reg <- ggplot(biomass[biomass$Year != "2020", ], aes(y = Across_PC1_CWM, x = Explo, fill = Explo)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                    labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Slow-fast continuum") + xlab(NULL) +
  #annotation_raster(Icon_poapra, xmin = 0.4, ymin = 3, xmax = 0.8, ymax = 4.5) +
  #annotation_raster(Icon_censol, xmin = 0.43, ymin = -6, xmax = 0.8, ymax = -4.5) +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        axis.text.x.bottom = element_blank(), axis.ticks = element_blank(),
        legend.position = "none") 


Slowfast_vs_reg <- ggdraw(Slowfast_vs_reg) +
  draw_image(image = Icon_eupesu, x = 0.14, y = 0.85, valign = 0, width = 0.1, height = 0.13) +
  draw_image(image = Icon_fesspp, x = 0.14, y = 0.02, valign = 0, width = 0.15, height = 0.15)


#biomass vs functional diversity
Biomass_vs_fd <- ggplot(biomass[biomass$Year != "2020", ], aes(x = Yr_Comb_rao_FP, y = biomass)) +
  geom_point(aes(col = Explo), alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, col = "black", lwd = 2, lty = "dotdash") +
  geom_smooth(aes(col = Explo), method = "lm", formula = y ~ x, se = F, lwd = 1.2) +
  scale_color_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                     labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Biomass") + xlab("Functional diversity") +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))

#box-plots of functional diversity
Fd_vs_reg <- ggplot(biomass[biomass$Year != "2020", ], aes(y = Yr_Comb_rao_FP, x = Explo, fill = Explo)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                    labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Functional diversity") + xlab(NULL) +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        axis.text.x.bottom = element_blank(), axis.ticks = element_blank())

#biomass vs species richness
Biomass_vs_sr <- ggplot(biomass[biomass$Year != "2020", ], aes(x = Yr_species_rich, y = biomass)) +
  geom_point(aes(col = Explo), alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, col = "black", lwd = 2, lty = "dotdash") +
  geom_smooth(aes(col = Explo), method = "lm", formula = y ~ x, se = F, lwd = 1.2) +
  scale_color_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                     labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Biomass") + xlab("Species richness") +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))

#box-plots of species richness
Sr_vs_reg <- ggplot(biomass[biomass$Year != "2020", ], aes(y = Yr_species_rich, x = Explo, fill = Explo)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                    labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Species richness") + xlab(NULL) +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        axis.text.x.bottom = element_blank(), axis.ticks = element_blank())

#biomass vs land-use intensity
Biomass_vs_lui <- ggplot(biomass[biomass$Year != "2020", ], aes(x = Yr_lui, y = biomass)) +
  geom_point(aes(col = Explo), alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, col = "black", lwd = 2, lty = "dotdash") +
  geom_smooth(aes(col = Explo), method = "lm", formula = y ~ x, se = F, lwd = 1.2) +
  scale_color_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                     labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Biomass") + xlab("Land-use intensity") +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))

#box-plots of land-use intensity
Lui_vs_reg <- ggplot(biomass[biomass$Year != "2020", ], aes(y = Yr_lui, x = Explo, fill = Explo)) +
  geom_boxplot() +
  scale_fill_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                    labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  ylab("Land-use intensity") + xlab(NULL) +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        axis.text.x.bottom = element_blank(), axis.ticks = element_blank())


#all correlation and box-plots together
Biom_vs_all_plots <- ggarrange(plotlist = list(Biomass_vs_slowfast, Slowfast_vs_reg,
                                               Biomass_vs_fd, Fd_vs_reg,
                                               Biomass_vs_sr, Sr_vs_reg,
                                               Biomass_vs_lui, Lui_vs_reg), nrow = 4, ncol = 2,
                               widths = c(1.5, 1), common.legend = T)

#unique plot with only relationships 
Biom_vs_all_corr <- ggarrange(plotlist = list(Biomass_vs_slowfast,
                                              Biomass_vs_fd,
                                              Biomass_vs_sr,
                                              Biomass_vs_lui), nrow = 2, ncol = 2, common.legend = T)

#unique plot with only box-plots 
Biom_vs_all_bx <- ggarrange(plotlist = list(Slowfast_vs_reg,
                                            Fd_vs_reg,
                                            Sr_vs_reg,
                                            Lui_vs_reg), nrow = 2, ncol = 2, common.legend = T)


#save plot
ggsave(plot = Biom_vs_all_corr, device = "jpeg", filename = "~/Documents/ClimateExtremes/CleanCode/New_figures/Biomass_vs_all_cor.jpeg",
       dpi = 500, units = "cm", width = 22, height = 22, bg = "white")

#save plot
ggsave(plot = Biom_vs_all_bx, device = "jpeg", filename = "~/Documents/ClimateExtremes/CleanCode/New_figures/Biomass_vs_all_bx.jpeg",
       dpi = 500, units = "cm", width = 22, height = 22, bg = "white")


#land-use intensity vs slow-fast continuum
Lui_vs_slowfast <- ggplot(biomass[biomass$Year != "2020", ], aes(y = Across_PC1_CWM, x = Yr_lui)) +
  geom_point(aes(col = Explo), alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, col = "black", lwd = 2, lty = "dotdash") +
  geom_smooth(aes(col = Explo), method = "lm", formula = y ~ x, se = F, lwd = 1.2) +
  scale_color_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                     labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  xlab("Land-use intensity") + ylab("Slow-fast continuum") +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        legend.position = "none")

Lui_vs_slowfast <- ggdraw(Lui_vs_slowfast) +
  draw_image(image = Icon_eupesu, x = 0.14, y = 0.85, valign = 0, width = 0.1, height = 0.13) +
  draw_image(image = Icon_fesspp, x = 0.14, y = 0.13, valign = 0, width = 0.15, height = 0.14)

#land-use intensity vs functional diversity
Lui_vs_fd <- ggplot(biomass[biomass$Year != "2020", ], aes(y = Yr_Comb_rao_FP, x = Yr_lui)) +
  geom_point(aes(col = Explo), alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, col = "black", lwd = 2, lty = "dotdash") +
  geom_smooth(aes(col = Explo), method = "lm", formula = y ~ x, se = F, lwd = 1.2) +
  scale_color_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                     labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  xlab("Land-use intensity") + ylab("Functional diversity") +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))

#land-use intensity vs species richness
Lui_vs_sr <- ggplot(biomass[biomass$Year != "2020", ], aes(y = Yr_species_rich, x = Yr_lui)) +
  geom_point(aes(col = Explo), alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = F, col = "black", lwd = 2, lty = "dotdash") +
  geom_smooth(aes(col = Explo), method = "lm", formula = y ~ x, se = F, lwd = 1.2) +
  scale_color_manual(values = c("ALB" = "#009E73", "HAI" = "#CC79A7", "SCH" = "#0072B2"),
                     labels = c("ALB" = "South-West", "HAI" = "Center", "SCH" = "North-East"), name = "Region:") +
  xlab("Land-use intensity") + ylab("Species richness") +
  theme_pubr() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 12),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14))

#all land-use intensity plots together
Lui_vs_all_corr <- ggarrange(Lui_vs_slowfast, Lui_vs_fd, Lui_vs_sr, nrow = 2, ncol = 2, common.legend = T)

#save plot
ggsave(plot = Lui_vs_all_corr, device = "jpeg", filename = "~/Documents/ClimateExtremes/CleanCode/New_figures/Lui_vs_all_cor.jpeg",
       dpi = 500, units = "cm", width = 22, height = 22, bg = "white")


##-----------------------Compute descriptive statistics on period of biomass harvest


#compute earliest day of biomass harvest per year and region
Min_day_harv <- with(biomass, tapply(X = day_of_year, INDEX = list(Explo, Year), min, na.rm = T))
#round results to integer
round(rowMeans(Min_day_harv[, -12]), digits = 0)

#compute latest day of biomass harvest per year and region
Max_day_harv <- with(biomass, tapply(X = day_of_year, INDEX = list(Explo, Year), max, na.rm = T))
round(rowMeans(Max_day_harv[, -12]), digits = 0)

#compute median of day of biomass harvest per year and region
Med_day_harv <- with(biomass, tapply(X = day_of_year, INDEX = list(Explo, Year), median, na.rm = T))
round(rowMeans(Med_day_harv[, -12]), digits = 0)

#compute lower interquartile bound (prob = 0.25) of day of biomass harvest per year and region
Low_iqr_day_harv <- with(biomass, tapply(X = day_of_year, INDEX = list(Explo, Year), quantile, probs = 0.25, na.rm = T))
round(rowMeans(Low_iqr_day_harv[, -12]), digits = 0)

#compute upper interquartile bound (prob = 0.75) of day of biomass harvest per year and region
Up_iqr_day_harv <- with(biomass, tapply(X = day_of_year, INDEX = list(Explo, Year), quantile, probs = 0.75, na.rm = T))
round(rowMeans(Up_iqr_day_harv[, -12]), digits = 0)



