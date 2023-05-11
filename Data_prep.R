##Code to replicate results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought" by Bazzichetto et al.
##Data preparation

library(dplyr)
library(ggplot2)
library(ggpubr)

#import biomass data
biomass <- read.csv(file = "~/Documents/ClimateExtremes/CleanCode/Data/biomass.csv",
                        header = T, sep = ",", dec = ".")


#keep only years used in the analyses
#biomass <- biomass[biomass$Year %in% 2009:2020, ]

#create column with region acronym: ALB (South-West), HAI (Center), SCH (North-East)
biomass$Explo <- substr(biomass$Useful_EP_PlotID, start = 1, stop = 1)

#create key to rename regions
Explo_key <- c(H = "HAI", A = "ALB", S = "SCH")

biomass$Explo <- unname(Explo_key[biomass$Explo])

##----------------------------------------------------------------Computation of log ratios

#compute LogR (annual log-ratio)
biomass.logratio.df <- Compute_LogR(biomass)

#join biomass and biomass.logratio.df
biomass <- dplyr::inner_join(biomass, biomass.logratio.df, by = c("EP_PlotID" = "PlotID", "Explo" = "Expl", "Year" = "Year"))

#compute LogR_reg_ref (regional-reference LogR)

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

#compute LogR_reg_ref
biomass$LogR_trend <- with(biomass, log(biomass/biom_trend))

#coerce year to integer
biomass.med.q$year <- as.integer(biomass.med.q$year)

#compute LogR_plot_ref (plot-reference LogR)

#compute plot-specific median value of biomass over the time-series
biomass.med.plot <- tapply(biomass[["biomass"]], biomass[["Useful_EP_PlotID"]], median, na.rm = T)

biomass.med.plot <- data.frame(LogR_plot = unname(biomass.med.plot),
                               Useful_EP_PlotID = names(biomass.med.plot))

#join plot-specific median value of biomass to biomass data.frame
biomass <- dplyr::inner_join(x = biomass, y = biomass.med.plot, by = c("Useful_EP_PlotID"))

#compute LogR_plot_ref (plot-reference LogR)
biomass$LogR_tr_plot <- with(biomass, log(biomass/LogR_plot))

#get rid of unused column
biomass <- biomass[!(colnames(biomass) %in% c("EP_PlotID"))]

##--------------------------------------------Import and process data on bio and non-biological variables

#Import data on soil humidity
HumData <- read.csv(file = "/Users/manbaz/Documents/ClimateExtremes/SoilHumidity/plots.csv",
                    sep = ",", dec = ".", stringsAsFactors = F)

#select observations for years used in the analyses
HumData <- HumData[(HumData$datetime %in% seq(2009, 2020, 1)), ]

#compute plot-specific mean value of soil humidity at 10 cm depth
Plot_spec_SH10 <- tapply(HumData$SM_10, INDEX = list(HumData$plotID), mean, na.rm = T)

#import (yearly data on functional traits)
##only keep data with no _pld
load("~/Documents/ClimateExtremes/CleanCode/Data/fundata.RData")

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

#join yearly data on functional traits to biomass data.frame
biomass$Year <- as.character(biomass$Year)

#FOR ME - here I inserted this line as (probably) biomass had no Exploratory column, which was in Df_grasslands (not imported here)
biomass$Exploratory <- ifelse(biomass$Explo == "ALB", "A", 
                              ifelse(biomass$Explo == "HAI", "H", "S"))

biomass <- dplyr::left_join(x = biomass, y = Yr_join_compl, c("Useful_EP_PlotID", "Year", "Exploratory"))

#remove useless vars
biomass$Hab <- NULL
biomass$Yr_clonal_0 <- NULL

#reorder vars
biomass <- biomass[c("Useful_EP_PlotID", "Year", "day_of_year", "Explo", "Exploratory", colnames(biomass)[!(colnames(biomass) %in% c("Useful_EP_PlotID", "Year", "day_of_year", "Explo", "Exploratory"))])]

#import data on yearly species richness and land use index (lui)
load("~/Documents/ClimateExtremes/Func_Phylo_Div/rich+lui.RData") 

#coerce lui_separate to data.frame
lui_separate <- as.data.frame(lui_separate)

#harmonize plot names
ID_not_in_PlotID <- unique(richness_plot_year$Useful_EP_PlotID[which(is.na(match(richness_plot_year$Useful_EP_PlotID, lui_separate$PlotID)))])
PlotID_not_in_ID <- unique(lui_separate$PlotID[which(is.na(match(lui_separate$PlotID, richness_plot_year$Useful_EP_PlotID)))])
ID_to_match <- setNames(ID_not_in_PlotID, PlotID_not_in_ID)

lui_separate$Useful_EP_PlotID <- lui_separate$PlotID

for(i in names(ID_to_match)) {
  lui_separate$Useful_EP_PlotID[lui_separate$Useful_EP_PlotID %in% i] <- ID_to_match[[i]] #use [[]] here
}

#join data.frame richness and lui
Rich_lui_j <- dplyr::left_join(richness_plot_year, lui_separate, by = c("Useful_EP_PlotID", "Year"))

#rename columns 
colnames(Rich_lui_j)[3] <- "Yr_species_rich"
colnames(Rich_lui_j)[9] <- "Yr_lui"

#join data.frame lui, species richness and biomass
biomass <- dplyr::left_join(x = biomass,
                            y = Rich_lui_j[c("Useful_EP_PlotID", "Year", "Yr_species_rich",
                                             "Yr_lui")],
                            by = c("Useful_EP_PlotID", "Year"))

biomass$Yr_species_rich <- unname(biomass$Yr_species_rich)

##--------------------------------------------Compute yearly value of slow-fast continuum using CWM

#PCA across years
Across_yrs_PC1 <- Compute_PCA_CWM(biomass)

#multiply Across_PC1_CWM by -1 to have slow-fast axis
Across_yrs_PC1$Across_PC1_CWM <- -1*(Across_yrs_PC1$Across_PC1_CWM)

#join result with biomass data.frame
biomass <- dplyr::left_join(x = biomass, y = Across_yrs_PC1, by = c("Useful_EP_PlotID", "Year"))

##--------------------------------------------Join Soil Humidity
biomass <- dplyr::left_join(x = biomass, y = data.frame(SH10 = unname(Plot_spec_SH10),
                                                            PlotID = names(Plot_spec_SH10)),
                            by = c("Useful_EP_PlotID" = "PlotID"))


##--------------------------------------------Have a look at biomass trend

#biomass.med.q for plotting
biomass.med.q.plot <- biomass.med.q

#add NAs for 2008 to match SPEI plot
biomass.med.q.plot <- rbind(biomass.med.q.plot, 
                            data.frame(Explo = c("ALB", "HAI", "SCH"), year = 2008L, Avg = NA, Med = NA, L_q = NA, U_q = NA))

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
  theme(legend.position = "none", text = element_text(size = 10),
        axis.text.x = element_text(size = 9, angle = 45, vjust = 1, hjust = 1))


