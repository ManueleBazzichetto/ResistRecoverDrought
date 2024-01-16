##Code to replicate results presented in the manuscript "Biodiversity promotes resistance but dominant species shape recovery of grasslands under extreme drought" by Bazzichetto et al.
##Functions

##----------------------------Functions used in the analyses

#function to compute LogR (annual log ratio)

Compute_LogR <- function(x) {
  df <- x[c("EP_PlotID", "Explo", "Year", "biomass")]
  df.full <- do.call(rbind, lapply(unique(df$Explo), function(.) {
    df.expl <- df[df$Explo == ., ]
    ratio.df <- do.call(rbind, lapply(unique(df.expl$EP_PlotID), function(plot) { 
      df.expl.plot <- df.expl[df.expl$EP_PlotID == plot, ]
      df.expl.plot <- df.expl.plot[order(df.expl.plot$Year), ]
      year <- df.expl.plot[["Year"]]
      log.r <- log((df.expl.plot[["biomass"]][2:(length(year))])/(df.expl.plot[["biomass"]][1:(length(year)-1)]))
      log.r <- c(NA_real_, log.r)
      #add another measure of biomass change, weighted difference biomass between t1 - t0
      w.dif <- ((df.expl.plot[["biomass"]][2:(length(year))]) - (df.expl.plot[["biomass"]][1:(length(year)-1)]))/(df.expl.plot[["biomass"]][1:(length(year)-1)])
      w.dif <- c(NA_real_, w.dif)
      df.res <- data.frame(LogR = log.r, Wdif = w.dif, Year = year, PlotID = plot, Expl = .)
    }))
    return(ratio.df)
  }))
  return(df.full)
  }


#function to perform PCA on CWM, and extract 1st axis

Compute_PCA_CWM <- function(dt) {
  dt <- dt[c("Yr_LDMC", "Yr_SLA", "Yr_LeafP",
             "Yr_Height_cmlog", "Yr_LeafN_log",
             "Yr_Seed_mass_log", "Useful_EP_PlotID", "Year")]
  NAs <- sort(unique(which(is.na(dt), arr.ind = T)[, 1]))
  ID_notNAs <- dt[-NAs, c("Useful_EP_PlotID", "Year")]
  dt <- na.omit(dt)
  PCA_CWM <- prcomp(x = dt[!(colnames(dt) %in% c("Useful_EP_PlotID", "Year"))],
                    center = T, scale. = T)
  biplot(PCA_CWM) #check whether first axis has to be multiplied by -1 
  PC1 <- PCA_CWM$x[, 1]
  PC1.df <- data.frame(Across_PC1_CWM = PC1, Useful_EP_PlotID = ID_notNAs[[1]],
                       Year = ID_notNAs[[2]])
  return(PC1.df)
  }

#modified from Compute_PCA_CWM - the function below simply re-run the same PCA
#but returns the prcomp object

Extract_PCA_CWM <- function(dt) {
  dt <- dt[c("Yr_LDMC", "Yr_SLA", "Yr_LeafP",
             "Yr_Height_cmlog", "Yr_LeafN_log",
             "Yr_Seed_mass_log", "Useful_EP_PlotID", "Year")]
  dt <- na.omit(dt)
  PCA_CWM <- prcomp(x = dt[!(colnames(dt) %in% c("Useful_EP_PlotID", "Year"))],
                    center = T, scale. = T)
  return(PCA_CWM)
  }

#function to find years in which to test recovery

Find_recovery_yrs <- function(x) {
  reslt <- do.call(rbind, lapply(unique(x$Explo), function(expl) {
    df <- x[x$Explo == expl, ]
    df$Year <- as.numeric(df$Year)
    Years_for_recovery <- vector(mode = "numeric")
    for(i in (unique(df$Year)[-c(11)])) {
      if(isTRUE(df[df$Year == (i), "SPEI_final_cat"] %in% c("ModerateDry", "ExtremeDry")) && isTRUE(df[df$Year == (i+1), "SPEI_final_cat"] %in% c("Normal", "ModerateWet", "ExtremeWet"))) {
        Years_for_recovery <- c(Years_for_recovery, i+1)
        } else {
          next
        }
      }
    Years_for_recovery <- as.character(Years_for_recovery)
    Years_for_recovery <- data.frame(Yrs = Years_for_recovery, Explo = expl)
    return(Years_for_recovery)
    }))
  return(reslt)
  }


#--------------

#the 2 functions below are used to plot the temporal autocorrelation function of residuals from
#the time-series models presented in the script "Time_series_log_ratios"
#this is can be achieved by nlme::ACF, although ACF does not consider gaps in the time-series
#the functions below allow computing the temporal ACF excluding gaps in the time-series
#this avoids estimating the ACF with non-consecutive time-steps

#WARNING! The functions below are somehow experimental and shouldn't be used outside this project
#the functions were used to check how much the ACF computed by nlme would deviate from that
#accounting for gaps in the time-series
#The datasets used in this project have only very few gaps in the time-series and the 
#almost totality of these gaps are of 1 time-stpe (i.e., 1 year)
#so, I am assuming the ACF estimated using nlme is reasonably reliable


#Compute_ACF is internally called by Correct_acf_gaps (see below)
#this function computes group-specific components of the ACF formula
#these components are then aggregated among groups (here, plot identifier)
#the function should provide sensible results only in case of few gaps in the time-series

Compute_ACF <- function(gr, lag) {
  gr_yrs <- gr[["Year"]]
  gr_yrs_nm <- setNames(gr[["Resi"]], as.character(gr_yrs))
  #create full seq with gaps
  sq_yrs <- seq(min(gr_yrs), max(gr_yrs), 1)
  sq_yrs_nm <- setNames(rep(NA, length(sq_yrs)), as.character(sq_yrs))
  #create seq to use
  sq_yrs_nm[names(gr_yrs_nm)] <- unname(gr_yrs_nm)
  res <- do.call(rbind, lapply(lag, function(i) {
    #prod of residuals at a given lag
    prod_num <- unname(sq_yrs_nm[1:(length(sq_yrs_nm)-i)] * sq_yrs_nm[(i + 1):(length(sq_yrs_nm))])
    n_num <- length(prod_num[!is.na(prod_num)])
    sum_num <- sum(prod_num, na.rm = T)
    prod_den <- sum(sq_yrs_nm^2, na.rm = T)
    n_denom <- sum(!is.na(sq_yrs_nm))
    compon_form <- data.frame(Nrt = sum_num, Dnm = prod_den, N_nrt = n_num, N_dnm = n_denom, lag = i)
    return(compon_form)
  }))
  return(res)
}

#Correct_acf_gaps compute the ACF for a vector of time lags
#basically, the function is designed to get a data.frame of residuals associated with groups (here, plot identifier),
#split residuals by group, and pass them to Compute_ACF to return the ACF value for each lag
#also, the function produces a plot showing the temporal ACF (as done by nlme::ACF)

Correct_acf_gaps <- function(dtf, lag, alpha_lev = 0.05) {
  require(ggplot2)
  require(ggpubr)
  #PlotID_col <- dtf[["PlotID"]]
  list_plots <- lapply(unique(dtf[["PlotID"]]), function(id) dtf[dtf[["PlotID"]] == id, ])
  acf_components <- do.call(rbind, lapply(list_plots, Compute_ACF, lag = lag))
  acf_vals <- do.call(rbind, lapply(unique(acf_components[["lag"]]), function(lg) {
    acf_lg <- acf_components[acf_components[["lag"]] == lg, c(1, 2, 3, 4)]
    acf_lg <- colSums(acf_lg)
    #compute acf
    acf_value <- (acf_lg["Nrt"]/acf_lg["N_nrt"])/(acf_lg["Dnm"]/acf_lg["N_dnm"])
    acf_lim <- qnorm(p = (1-(alpha_lev/2)), mean = 0, sd = 1, lower.tail = T, log.p = F)
    rslt <- data.frame(lag = lg, ACF = acf_value, N_nrt = acf_lg[["N_nrt"]],
                       ACF_low = -1*(acf_lim/sqrt(acf_lg[["N_nrt"]])),
                       ACF_up = (acf_lim/sqrt(acf_lg[["N_nrt"]])))
    return(rslt)
  }))
  #plot ACF
  Plot_ACF <- ggplot(acf_vals, aes(x = lag, y = ACF)) +
    geom_ribbon(aes(ymin = ACF_low, ymax = ACF_up), alpha = 0.3) +
    geom_point(size = 3, alpha = .6) +
    geom_line() +
    ylab("Temporal ACF") + xlab("Lag") +
    theme_pubclean() +
    theme(text = element_text(size = 14))
  plot(Plot_ACF)
  return(acf_vals)
}


