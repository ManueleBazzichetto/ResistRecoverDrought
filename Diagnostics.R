##Code to replicate results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought" by Bazzichetto et al.
##Diagnostics of year-by-year models

library(ggplot2)
library(ggpubr)
library(car)
library(performance)


##-------------------------------------------------Across years correlation between land use and biodiversity

lapply(unique(biomass$Year), function(yr) {
  df <- biomass[biomass$Year == yr, ]
  lapply(unique(df$Explo), function(.) {
    df <- df[df$Explo == ., c("Yr_lui", "Yr_Comb_rao_FP", "Across_PC1_CWM", "Yr_species_rich")]
    with(df, cor(Filter(is.numeric, df), use = "pairwise.complete.obs"))
  })
})

##-------------------Across years multi-collinearity between species richness, land use and functional components

#fit models considering functional composition and diversity, and species richness (all)
Across_yrs_regression.all <- lapply(c("LogR", "LogR_tr_plot", "LogR_trend"), function(resp) {
  lapply(switch(resp, "LogR_trend" = unique(biomass$Year)[-12],
                "LogR" = unique(biomass$Year)[-c(1, 12)],
                "LogR_tr_plot" = unique(biomass$Year)[-12]), function(i) {
                  reslts <- Yr_regression(biomass[biomass$Year == i, ], y = resp,
                                          Xs = c("Explo", "Yr_lui",
                                                 "Yr_species_rich",
                                                 "Yr_Comb_rao_FP",
                                                 "day_of_year",
                                                 "Across_PC1_CWM", 
                                                 "SH10"), std = T)
                  return(reslts)
                })
  })

names(Across_yrs_regression.all) <- c("LogR", "LogR_tr_plot", "LogR_trend")
names(Across_yrs_regression.all$LogR) <- paste0("Yr", seq(2010, 2019))
names(Across_yrs_regression.all$LogR_tr_plot) <- paste0("Yr", seq(2009, 2019))
names(Across_yrs_regression.all$LogR_trend) <- paste0("Yr", seq(2009, 2019))


#check multi-collinearity
lapply(Across_yrs_regression.all, function(i) {
  lapply(i, function(.) {
    vif(.[[3]])
  })
})


##-------------------Average (across years) correlation between species richness, land use and functional components

Cor_Sr_all <- colMeans(do.call(rbind, lapply(unique(biomass$Year)[-12], function(yr) {
  df <- biomass[biomass$Year == yr, c("Yr_species_rich", "Across_PC1_CWM", "Yr_Comb_rao_FP", "Yr_lui")]
  c(Sr_sf = with(df, cor(Yr_species_rich, Across_PC1_CWM, use = "pairwise.complete.obs")),
    Sr_fd = with(df, cor(Yr_species_rich, Yr_Comb_rao_FP, use = "pairwise.complete.obs")),
    Sr_lui = with(df, cor(Yr_species_rich, Yr_lui, use = "pairwise.complete.obs")))
  })))

round(Cor_Sr_all, digits = 2)

##-------------------------------------------------Diagnostics of year-by-year models

##----------------------------------------------------------------------------LogR

lapply(Across_yrs_regression$LogR, function(i) check_model(i[[3]]))

##----------------------------------------------------------------------------LogR ref-plot

lapply(Across_yrs_regression$LogR_tr_plot, function(i) check_model(i[[3]]))

##----------------------------------------------------------------------------LogR ref-reg

lapply(Across_yrs_regression$LogR_trend, function(i) check_model(i[[3]]))


##-------------------------------------------------Marginal relationships between resistance, biodiversity and drought intensity

##SPEI-3

#FD
ggplot(TA_LogR_trpl_SPEI3, aes(x = Yr_Comb_rao_FP, y = LogR_tr_plot, col = Year)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(Explo ~ SPEI_final_cat)

ggplot(TA_LogR_trpl_SPEI3, aes(x = Yr_Comb_rao_FP, y = LogR_tr_plot)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(Explo ~ SPEI_final_cat)

#Slow-fast
ggplot(TA_LogR_trpl_SPEI3, aes(x = Across_PC1_CWM, y = LogR_tr_plot)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(Explo ~ SPEI_final_cat)

#SR
ggplot(TA_LogR_trpl_SPEI3, aes(x = Yr_species_rich, y = LogR_tr_plot)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(Explo ~ SPEI_final_cat)


##SPEI-12

#FD
ggplot(TA_LogR_trpl_SPEI12, aes(x = Yr_Comb_rao_FP, y = LogR_tr_plot)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(Explo ~ SPEI_final_cat)

#Slow-fast
ggplot(TA_LogR_trpl_SPEI12, aes(x = Across_PC1_CWM, y = LogR_tr_plot)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(Explo ~ SPEI_final_cat)

#SR
ggplot(TA_LogR_trpl_SPEI12, aes(x = Yr_species_rich, y = LogR_tr_plot)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(Explo ~ SPEI_final_cat)

