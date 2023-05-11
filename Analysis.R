##Code to replicate results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought" by Bazzichetto et al.
##Statistical analyses - year-by-year models

#for plotting
library(ggplot2)
library(ggpubr)
library(patchwork)
#library(GGally)
#for modelling
library(car)

##----------------------------------------------------------------Fit yearly statistical models

##----------------------------------------------------------------Fit yearly models considering functional composition and diversity

#labels for predictors functional models
Lab_vars.f <- c('Across_PC1_CWM' = "Slow-fast continuum",
                    'Yr_Comb_rao_FP' = "Functional diversity",
                    'Yr_lui' = "Land-use intensity",
                    'Yr_Myco_intensity' = "Mycorrhizae int.", 
                    "SH10" = "Soil humidity")

#fit models considering functional composition and diversity
Across_yrs_regression <- lapply(c("LogR", "LogR_tr_plot", "LogR_trend"), function(resp) {
  lapply(switch(resp, "LogR_trend" = unique(biomass$Year)[-12],
                "LogR" = unique(biomass$Year)[-c(1, 12)],
                "LogR_tr_plot" = unique(biomass$Year)[-12]), function(i) {
                  reslts <- Yr_regression(biomass[biomass$Year == i, ], y = resp,
                                          Xs = c("Explo", "Yr_lui",
                                                 "Yr_Comb_rao_FP",
                                                 "Yr_Myco_intensity",
                                                 "Across_PC1_CWM", 
                                                 "SH10"), std = T)
                  return(reslts)
                })
  })

names(Across_yrs_regression) <- c("LogR", "LogR_tr_plot", "LogR_trend")
names(Across_yrs_regression$LogR) <- paste0("Yr", seq(2010, 2019))
names(Across_yrs_regression$LogR_tr_plot) <- paste0("Yr", seq(2009, 2019))
names(Across_yrs_regression$LogR_trend) <- paste0("Yr", seq(2009, 2019))


#check multi-collinearity
lapply(Across_yrs_regression, function(i) {
  lapply(i, function(.) {
    vif(.[[3]])
  })
})


#create matrix with confidence intervals for yearly models
CI_mat_yearsp <- do.call(rbind, lapply(names(Across_yrs_regression), function(nm) {
  mat <- do.call(rbind, lapply(Across_yrs_regression[[nm]], '[[', 1))
  mat$Resp <- nm
  return(mat)
})
)

CI_mat_yearsp$Coef_name[grep(pattern = "(Intercept)", x = CI_mat_yearsp$Coef_name)] <- "ExploALB"

CI_mat_yearsp$Sign <- ifelse(CI_mat_yearsp$Upper_endpoint > 0 & CI_mat_yearsp$Lower_endpoint < 0, "NO", "YES")

#plot models' results

#change order in CI_mat_yearsp
CI_mat_yearsp$Coef_name_f <- factor(CI_mat_yearsp$Coef_name, levels = c("Across_PC1_CWM", "Yr_Comb_rao_FP",
                                                                        "Yr_lui", "Yr_Myco_intensity", "SH10",
                                                                        "ExploALB", "ExploHAI", "ExploSCH"))

CI_mat_yearsp$Resp_f <- factor(CI_mat_yearsp$Resp, levels = c("LogR", "LogR_tr_plot", "LogR_trend"))


CIplot_across <- ggplot(CI_mat_yearsp[(!(CI_mat_yearsp$Coef_name %in% c("ExploALB", "ExploHAI", "ExploSCH"))), ],
                        aes(y = Estimate, x = as.factor(Year), col = Sign)) +
  #try to change background of panel for LogR_trend, which will be used as control
  #from: https://stackoverflow.com/questions/9847559/conditionally-change-panel-background-with-facet-grid
  geom_rect(data = CI_mat_yearsp[(!(CI_mat_yearsp$Coef_name %in% c("ExploALB", "ExploHAI", "ExploSCH"))) & CI_mat_yearsp$Resp == "LogR_trend", ],
            col = "black", fill = "grey70", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_rect(data = CI_mat_yearsp[(!(CI_mat_yearsp$Coef_name %in% c("ExploALB", "ExploHAI", "ExploSCH"))) & CI_mat_yearsp$Resp != "LogR_trend", ],
            col = "black", fill = "grey94", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, lty = "solid", col = "grey70") +
  geom_point() +
  geom_errorbar(aes(ymin = Lower_endpoint, ymax = Upper_endpoint), width = .3) +
  ylab(NULL) + xlab(NULL) +
  scale_color_discrete(type = c("NO" = "grey54", "YES" = "blue3")) +
  facet_grid(Coef_name_f ~ Resp_f, labeller = labeller(Coef_name_f = Lab_vars.f, 
                                                       Resp_f = c("LogR" = "LogR",
                                                                  "LogR_trend" = "LogR ref-reg",
                                                                  "LogR_tr_plot" = "LogR ref-plot"))) +
  theme_pubr() + theme(legend.position = "none", axis.text.y = element_text(size = 12), 
                       strip.text.x = element_text(size = 12, face = "bold"),
                       strip.text.y = element_text(size = 12, face = "bold"),
                       axis.text.x.bottom = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                       axis.text.y.left = element_text(size = 9), strip.background = element_blank())


##----------------------------------------------------------------Fit yearly models considering species richness


#labels for predictors species richness models
Lab_vars.sr <- c('Yr_species_rich' = "Species richness",
                'Yr_lui' = "Land-use intensity",
                'Yr_Myco_intensity' = "Mycorrhizae int.", 
                "SH10" = "Soil humidity")

#fit models considering species richness
Across_yrs_regression.sr <- lapply(c("LogR", "LogR_tr_plot", "LogR_trend"), function(resp) {
  lapply(switch(resp, "LogR_trend" = unique(biomass$Year)[-12],
                "LogR" = unique(biomass$Year)[-c(1, 12)],
                "LogR_tr_plot" = unique(biomass$Year)[-12]), function(i) {
                  reslts <- Yr_regression(biomass[biomass$Year == i, ], y = resp,
                                          Xs = c("Explo", "Yr_lui",
                                                 "Yr_species_rich",
                                                 "Yr_Myco_intensity",
                                                 "SH10"), std = T)
                  return(reslts)
                })
})

names(Across_yrs_regression.sr) <- c("LogR", "LogR_tr_plot", "LogR_trend")
names(Across_yrs_regression.sr$LogR) <- paste0("Yr", seq(2010, 2019))
names(Across_yrs_regression.sr$LogR_tr_plot) <- paste0("Yr", seq(2009, 2019))
names(Across_yrs_regression.sr$LogR_trend) <- paste0("Yr", seq(2009, 2019))


#check multi-collinearity
lapply(Across_yrs_regression.sr, function(i) {
  lapply(i, function(.) {
    vif(.[[3]])
  })
})


#create matrix with confidence intervals for yearly models
CI_mat_yearsp.sr <- do.call(rbind, lapply(names(Across_yrs_regression.sr), function(nm) {
  mat <- do.call(rbind, lapply(Across_yrs_regression.sr[[nm]], '[[', 1))
  mat$Resp <- nm
  return(mat)
})
)

CI_mat_yearsp.sr$Coef_name[grep(pattern = "(Intercept)", x = CI_mat_yearsp.sr$Coef_name)] <- "ExploALB"

CI_mat_yearsp.sr$Sign <- ifelse(CI_mat_yearsp.sr$Upper_endpoint > 0 & CI_mat_yearsp.sr$Lower_endpoint < 0, "NO", "YES")

#plot models' results

#change order in CI_mat_yearsp
CI_mat_yearsp.sr$Coef_name_f <- factor(CI_mat_yearsp.sr$Coef_name, levels = c("Yr_species_rich", "Yr_lui", "Yr_Myco_intensity", "SH10",
                                                                        "ExploALB", "ExploHAI", "ExploSCH"))

CI_mat_yearsp.sr$Resp_f <- factor(CI_mat_yearsp.sr$Resp, levels = c("LogR", "LogR_tr_plot", "LogR_trend"))

#plot
CIplot_across.sr <- ggplot(CI_mat_yearsp.sr[!(CI_mat_yearsp.sr$Coef_name %in% c("ExploALB", "ExploHAI", "ExploSCH")), ],
                           aes(y = Estimate, x = as.factor(Year), col = Sign)) +
  geom_rect(data = CI_mat_yearsp.sr[(!(CI_mat_yearsp.sr$Coef_name %in% c("ExploALB", "ExploHAI", "ExploSCH"))) &
                                      CI_mat_yearsp.sr$Resp == "LogR_trend", ],
            col = "black", fill = "grey70", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_rect(data = CI_mat_yearsp.sr[(!(CI_mat_yearsp.sr$Coef_name %in% c("ExploALB", "ExploHAI", "ExploSCH"))) &
                                      CI_mat_yearsp.sr$Resp != "LogR_trend", ],
            col = "black", fill = "grey94", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, lty = "solid", col = "grey70") +
  geom_point() +
  geom_errorbar(aes(ymin = Lower_endpoint, ymax = Upper_endpoint), width = .3) +
  ylab(NULL) + xlab(NULL) +
  scale_color_discrete(type = c("NO" = "grey54", "YES" = "blue3", "HCCM" = "black")) +
  facet_grid(Coef_name_f ~ Resp, labeller = labeller(Coef_name_f = Lab_vars.sr, 
                                                      Resp = c("LogR" = "LogR",
                                                                  "LogR_trend" = "LogR ref-reg",
                                                                  "LogR_tr_plot" = "LogR ref-plot"))) +
  theme_pubr() + theme(legend.position = "none", axis.text.y = element_text(size = 12), 
                       strip.text.x = element_text(size = 12, face = "bold"),
                       strip.text.y = element_text(size = 12, face = "bold"),
                       axis.text.x.bottom = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                       axis.text.y.left = element_text(size = 9), strip.background = element_blank())




#plot merging results for functional composition and diversity, and species richness


CI_mat.combined <- rbind(data.frame(CI_mat_yearsp[CI_mat_yearsp$Coef_name %in% c("Yr_Comb_rao_FP", "Across_PC1_CWM"), ], Mod = "F"),
                         data.frame(CI_mat_yearsp.sr[CI_mat_yearsp.sr$Coef_name %in% c("Yr_species_rich"), ], Mod = "SR"))


CI_mat.combined$Coef_name_f <- factor(CI_mat.combined$Coef_name_f, levels = c("Across_PC1_CWM", "Yr_Comb_rao_FP", "Yr_species_rich"))

CI_mat.combined$Color <- with(CI_mat.combined, paste(Sign, Mod, sep = "_"))

F_SR_combined <- ggplot(CI_mat.combined,
       aes(y = Estimate, x = as.factor(Year), col = Color)) +
  geom_rect(data = CI_mat.combined[CI_mat.combined$Resp == "LogR_trend", ],
            col = "black", fill = "grey70", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_rect(data = CI_mat.combined[CI_mat.combined$Resp != "LogR_trend", ],
            col = "black", fill = "grey94", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, lty = "solid", col = "grey70") +
  geom_point() +
  geom_errorbar(aes(ymin = Lower_endpoint, ymax = Upper_endpoint), width = 0) +
  ylab(NULL) + xlab(NULL) +
  scale_color_discrete(type = c("NO_F" = "grey54", "NO_SR" = "grey54", "YES_F" = "blue3", "YES_SR" = "red")) +
  facet_grid(Coef_name_f ~ Resp, labeller = labeller(Coef_name_f = c(Lab_vars.f, "Yr_species_rich" = "Species richness"), 
                                                     Resp = c("LogR" = "LogR",
                                                              "LogR_trend" = "LogR ref-reg",
                                                              "LogR_tr_plot" = "LogR ref-plot"))) +
  theme_pubr() + theme(legend.position = "none", axis.text.y = element_text(size = 12), 
                       strip.text.x = element_text(size = 12, face = "bold"),
                       strip.text.y = element_text(size = 12, face = "bold"),
                       axis.text.x.bottom = element_blank(),
                       axis.text.y.left = element_text(size = 9), strip.background = element_blank())


CI_lui.combined <- rbind(data.frame(CI_mat_yearsp[CI_mat_yearsp$Coef_name %in% c("Yr_lui"), ], Mod = "F"),
                         data.frame(CI_mat_yearsp.sr[CI_mat_yearsp.sr$Coef_name %in% c("Yr_lui"), ], Mod = "SR"))


CI_lui.combined$Coef_name_f <- factor(CI_lui.combined$Coef_name_f, levels = c("Yr_lui"))

CI_lui.combined$Color <- with(CI_lui.combined, paste(Sign, Mod, sep = "_"))

LUI_combined <- ggplot(CI_lui.combined,
       aes(y = Estimate, x = as.factor(Year), col = Color)) +
  geom_rect(data = CI_lui.combined[CI_lui.combined$Resp == "LogR_trend", ],
            col = "black", fill = "grey70", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_rect(data = CI_lui.combined[CI_lui.combined$Resp != "LogR_trend", ],
            col = "black", fill = "grey94", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, lty = "solid", col = "grey70") +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = Lower_endpoint, ymax = Upper_endpoint), width = 0, position = position_dodge(width = 1)) +
  ylab(NULL) + xlab(NULL) +
  scale_color_discrete(type = c("NO_F" = "grey54", "NO_SR" = "grey54", "YES_F" = "blue3", "YES_SR" = "red")) +
  facet_grid(Coef_name_f ~ Resp, labeller = labeller(Coef_name_f = c("Yr_lui" = "Land-use intensity"), 
                                                     Resp = c("LogR" = "LogR",
                                                              "LogR_trend" = "LogR ref-reg",
                                                              "LogR_tr_plot" = "LogR ref-plot"))) +
  theme_pubr() + theme(legend.position = "none", axis.text.y = element_text(size = 12), 
                       strip.text.x = element_blank(),
                       strip.text.y = element_text(size = 12, face = "bold"),
                       axis.text.x.bottom = element_text(size = 9, angle = 45, vjust = 1, hjust = 1),
                       axis.text.y.left = element_text(size = 9), strip.background = element_blank())


Yr_models_res <- F_SR_combined / LUI_combined + plot_layout(heights = c(2.5, 1))

#figure summarising model results, biomass and spei trend

Summary_fig <- ((Biomass_plot_trend / MainSPEI_plot) | Yr_models_res) + plot_annotation(tag_levels = "a")

ggsave(plot = Summary_fig, "~/Documents/ClimateExtremes/CleanCode/Summary_res.jpeg", device = "jpeg", dpi = 300, width = 32, height = 25, units = "cm")


#plot merging results for intensity of mycorrhizae and soil humidity for manuscript appendix

CI_myco_sh <- rbind(data.frame(CI_mat_yearsp[CI_mat_yearsp$Coef_name %in% c("Yr_Myco_intensity", "SH10"), ],
                               Mod = "F"),
                    data.frame(CI_mat_yearsp.sr[CI_mat_yearsp.sr$Coef_name %in% c("Yr_Myco_intensity", "SH10"), ],
                               Mod = "SR"))

CI_myco_sh$Color <- with(CI_myco_sh, paste(Sign, Mod, sep = "_"))

CI_myco_sh$Ungr_facet <- with(CI_myco_sh, paste(Coef_name, Mod, sep = "_"))

CI_myco_sh$Ungr_facet <- factor(CI_myco_sh$Ungr_facet,
                                 levels = c("Yr_Myco_intensity_F", "Yr_Myco_intensity_SR",
                                            "SH10_F", "SH10_SR"))

Myco_SH_combined <- ggplot(CI_myco_sh, aes(y = Estimate, x = as.factor(Year), col = Color)) +
  geom_rect(data = CI_myco_sh[CI_myco_sh$Resp == "LogR_trend", ],
            col = "black", fill = "grey70", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_rect(data = CI_myco_sh[CI_myco_sh$Resp != "LogR_trend", ],
            col = "black", fill = "grey94", xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_hline(yintercept = 0, lty = "solid", col = "grey70") +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = Lower_endpoint, ymax = Upper_endpoint), width = 0, position = position_dodge(width = 1)) +
  ylab(NULL) + xlab(NULL) +
  scale_color_discrete(type = c("NO_F" = "grey54", "NO_SR" = "grey54", "YES_F" = "blue3", "YES_SR" = "red")) +
  facet_grid(Ungr_facet ~ Resp, labeller = labeller(Ungr_facet = c("Yr_Myco_intensity_F" = "Mycorrhizae int.",
                                                                   "Yr_Myco_intensity_SR" = "Mycorrhizae int.",
                                                                   "SH10_F" = "Soil humidity",
                                                                   "SH10_SR" = "Soil humidity"), 
                                                     Resp = c("LogR" = "LogR",
                                                              "LogR_trend" = "LogR ref-reg",
                                                              "LogR_tr_plot" = "LogR ref-plot"))) +
  theme_pubr() + theme(legend.position = "none", axis.text.y = element_text(size = 12), 
                       strip.text.x = element_text(size = 12, face = "bold"),
                       strip.text.y = element_text(size = 12, face = "bold"),
                       axis.text.x.bottom = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
                       axis.text.y.left = element_text(size = 9), strip.background = element_blank())

ggsave(plot = Myco_SH_combined, "~/Documents/ClimateExtremes/CleanCode/Myco_SH_app.jpeg", device = "jpeg",
       dpi = 300, width = 18, height = 16, units = "cm")
