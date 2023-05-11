##Code to replicate results presented in the manuscript "Functional components of biodiversity mediate stability of grasslands under extreme drought" by Bazzichetto et al.
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


#function to fit linear regression and extract different info (vocariance matrix of coefficient estimators, confidence intervals for parameters and so on)

Yr_regression <- function(df, y, Xs, std = FALSE) {
  vars <- c(y, Xs)
  df <- df[c(vars, "Year", "Useful_EP_PlotID")]
  message(paste("Deleted", length(unique(which(is.na(df), arr.ind = T)[, 1])), "NAs", sep = " "))
  df <- na.omit(df)
  if(std) {
    #only standardize predictors
    df.std <- cbind(df[y], scale(Filter(is.numeric, df[!(colnames(df) %in% y)])), Filter(is.character, df))
  }
  form <- as.formula(paste(y, "~", paste(Xs, collapse = " + "), collapse = " ")) #only additive effect here
  if(std) {
    mod_std <- lm(form, data = df.std)
  }
  mod <- lm(form, data = df)
  if(std) {
    ConfInt <- confint.lm(mod_std)
    Coef <- coef(mod_std)
  } else {
    ConfInt <- confint.lm(mod)
    Coef <- coef(mod)
  }
  rownames(ConfInt) <- seq_len(nrow(ConfInt))
  ConfInt <- data.frame(unname(Coef), ConfInt,
                        Coef_name = names(Coef), Year = as.numeric(unique(df[["Year"]])))
  colnames(ConfInt)[1:3] <- c("Estimate", "Lower_endpoint", "Upper_endpoint")
  return(list(ConfInt = ConfInt, Vcov = vcov(mod), Mod = mod, Df = df))
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


