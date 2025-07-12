# Updated functions to Calculate partitionR() from comstab

# Author(s): Reilly O'Connor
# Version: 2024-06-18

# Pkgs
library(tidyverse)
library(zoo)

##### Code #####

#function to calcualte standard error
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

#Function to record warnings in dataframe from partition R
capture_warnings <- function(expr) {
  warnings <- character()
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings)
}

partitionR <- function (z, ny = 1) {
  if (!is.matrix(z)) 
    stop("Error: z is not a matrix")
  if (!is.numeric(z)) 
    stop("Error: non-numerical values in z")
  if (any(z < 0)) 
    stop("Error: negative values in z")
  if (dim(z)[1] == 1) 
    stop("Error: single-row matrix")
  if (!is.numeric(ny)) 
    stop("ny must be numeric")
  z[is.na(z)] <- 0
  z <- z[, colSums(z) > 0, drop = FALSE]
  z <- z[, apply(X = z, MARGIN = 2, FUN = min) != apply(X = z, MARGIN = 2, FUN = max), drop = FALSE]
  nyi <- apply(X = z, MARGIN = 2, FUN = function(x) sum(x > 0))
  z <- z[, nyi > ny, drop = FALSE]
  n <- ncol(z)
  varsum <- stats::var(rowSums(z))
  meansum <- mean(rowSums(z))
  CV <- sqrt(varsum)/meansum
  if (CV == 0) 
    stop("The community CV is zero. This analysis does not apply to \n                   perfectly stable communities.")
  if (dim(z)[2] == 1) {
    warning("This analysis is not relevant for single-species communities. \n            All stabilizing effects were fixed to 1.")
    TPLs <- stats::setNames(object = c(1, 1, 1, 1), nm = c("b","b_intercept", "TPL", "TPL_intercept"))
    CVs <- stats::setNames(object = c(CV, CV, CV, CV), nm = c("CVe", 
                                                              "CVtilde", "CVa", "CVc"))
    Stabilization <- stats::setNames(object = c(1, 1, 1, 
                                                1), nm = c("tau", "Delta", "Psi", "omega"))
    Relative <- stats::setNames(object = rep(NA, 3), nm = c("Delta_cont", 
                                                            "Psi_cont", "omega_cont"))
    res <- list(CVs = CVs, Stabilization = Stabilization, 
                Relative = Relative)
    class(res) <- "comstab"
    return(res)
  }
  if (dim(z)[2] > 1) {
    vari <- apply(X = z, MARGIN = 2, FUN = stats::var)
    meani <- colMeans(z)
    CVi <- sqrt(vari)/meani
    if (any(CVi == 0)) 
      warning("Non-fluctuating species found in the data.")
    CV0 <- which(CVi > 0)
    TPL <- stats::coef(stats::lm(log10(CVi[CV0]) ~ log10(meani[CV0])))
    b <- stats::coef(stats::lm(log10(vari) ~ log10(meani)))
    CVe <- 10^TPL[1] * (meansum/n)^TPL[2]
    if (sum(CV0) > 5) {
      testcor <- stats::cor.test(log10(vari[CV0]), log10(meani[CV0]))$p.value > 
        0.05
      if (testcor) 
        warning("No significant power law between species variance and mean abundances.")
    }
    else {
      warning("Low number of species. The power law between species variance and mean abundances cannot be tested.")
    }
    sumsd <- sum(sqrt(vari))
    CVtilde <- sumsd/meansum
    Delta <- CVtilde/CVe
    if (Delta > 1) 
      warning("Destabilizing effect of dominants. Relative effects cannot be computed.")
    sdsum <- sqrt(varsum)
    rootPhi <- sdsum/sumsd
    sumvar <- sum(vari)
    beta <- log10(1/2)/(log10(sumvar/(sumsd^2)))
    Psi <- rootPhi^beta
    omega <- rootPhi/Psi
    if (omega > 1) 
      warning("Community diversity is lower than the null diversity. Relative effects cannot be computed.")
    tau <- Delta * Psi * omega
    TPLs <- stats::setNames(object = c(b[2], b[1],TPL[2], TPL[1]), nm = c("b","b_intercept", "TPL", "TPL_intercept"))
    CVs <- stats::setNames(object = c(CVe, CVtilde, CVtilde * Psi, CV), 
                           nm = c("CVe", "CVtilde", "CVa", "CVc"))
    Stabilization <- stats::setNames(object = c(tau, Delta, Psi, omega), 
                                     nm = c("tau", "Delta", "Psi", "omega"))
    if (any(Stabilization > 1)) {
      Relative <- stats::setNames(object = rep(NA, 3), 
                                  nm = c("Delta_cont", "Psi_cont", "omega_cont"))
    }
    else {
      Relative <- stats::setNames(object = c(log10(Delta)/log10(tau), 
                                             log10(Psi)/log10(tau), log10(omega)/log10(tau)), 
                                  nm = c("Delta_cont", "Psi_cont", "omega_cont"))
    }
    res <- list(CVs = CVs, Stabilization = Stabilization, 
                Relative = Relative, TPLs = TPLs)
    class(res) <- "comstab"
    return(res)
  }
}


partitionR_window <- function(window_matrix) {
  
  #Convert to data frame to access column names
  window_df <- as.data.frame(window_matrix)
  
  start_year <- min(window_df$Year)
  end_year <- max(window_df$Year)
  
  community_matrix <- window_df %>% select(-Year) %>% as.matrix()
  
  captured <- capture_warnings(partitionR(community_matrix, ny = 1))
  res <- captured$result
  
  data.frame(
    start_year = start_year,
    end_year = end_year,
    CVe = as.numeric(res$CVs[["CVe"]]),
    CVtilde = as.numeric(res$CVs[["CVtilde"]]),
    CVa = as.numeric(res$CVs[["CVa"]]),
    CVc = as.numeric(res$CVs[["CVc"]]),
    tau = as.numeric(res$Stabilization[["tau"]]),
    Delta = as.numeric(res$Stabilization[["Delta"]]),
    Psi = as.numeric(res$Stabilization[["Psi"]]),
    omega = as.numeric(res$Stabilization[["omega"]]),
    Delta_cont = as.numeric(res$Relative[["Delta_cont"]]),
    Psi_cont = as.numeric(res$Relative[["Psi_cont"]]),
    omega_cont = as.numeric(res$Relative[["omega_cont"]]),
    vari_mean_slope = as.numeric(res$TPLs[["b"]]),
    cvi_mean_slope = as.numeric(res$TPLs[["TPL"]]),
    vari_mean_intercept = as.numeric(res$TPLs[["b_intercept"]]),
    cvi_mean_intercept = as.numeric(res$TPLs[["TPL_intercept"]]),
    warnings = paste(captured$warnings, collapse = "; ")
  )
}



plot_TPL <- function(df, window_size, ny = 1, form = c("vari", "CVi")) {
  form <- match.arg(form)
  cmat <- as.matrix(df %>% dplyr::select(-Year))
  
  start_year <- 1
  
  #Loop over each starting point of the window
  while (start_year <= (nrow(df) - window_size + 1)) {
    end_year <- start_year + window_size - 1
    
    #Subset the community matrix for the current window
    cmat_window <- cmat[start_year:end_year, ]
    
    #Initial validation and preprocessing steps from partitionR
    if (!is.matrix(cmat_window)) 
      stop("Error: cmat_window is not a matrix")
    if (!is.numeric(cmat_window)) 
      stop("Error: non-numerical values in cmat_window")
    if (any(cmat_window < 0)) 
      stop("Error: negative values in cmat_window")
    if (dim(cmat_window)[1] == 1) 
      stop("Error: single-row matrix")
    if (!is.numeric(ny)) 
      stop("ny must be numeric")
    cmat_window[is.na(cmat_window)] <- 0
    cmat_window <- cmat_window[, colSums(cmat_window) > 0, drop = FALSE]
    cmat_window <- cmat_window[, apply(X = cmat_window, MARGIN = 2, FUN = min) != apply(X = cmat_window, MARGIN = 2, FUN = max), drop = FALSE]
    nyi <- apply(X = cmat_window, MARGIN = 2, FUN = function(x) sum(x > 0))
    cmat_window <- cmat_window[, nyi > ny, drop = FALSE]
    
    if (dim(cmat_window)[2] > 1) {
      vari <- apply(X = cmat_window, MARGIN = 2, FUN = stats::var)
      meani <- colMeans(cmat_window)
      CVi <- sqrt(vari)/meani
      if (any(CVi == 0)) 
        warning("Non-fluctuating species found in the data.")
      CV0 <- which(CVi > 0)
      
      #Compute necessary statistics
      TPL <- stats::coef(stats::lm(log10(CVi[CV0]) ~ log10(meani[CV0])))
      b <- stats::coef(stats::lm(log10(vari) ~ log10(meani)))
      p_value_CVi <- stats::cor.test(log10(CVi[CV0]), log10(meani[CV0]))$p.value
      p_value_Vari <- stats::cor.test(log10(vari), log10(meani))$p.value
      corr_CVi <- stats::cor(log10(CVi[CV0]), log10(meani[CV0]))
      corr_Vari <- stats::cor(log10(vari), log10(meani))
      
      #Unlogged slopes
      slope_CVi <- TPL[2]
      slope_vari <- b[2]
      
      #Plotting based on selected form
      if (form == "vari") {
        plot(log10(meani), log10(vari), 
             main = paste("Taylor's Power Law: Variance vs Mean\nTime Window:", start_year, "-", end_year),
             xlab = "log10(Mean i)", ylab = "log10(Var i)",
             pch = 19, col = "blue")
        abline(lm(log10(vari) ~ log10(meani)), col = "red")
        abline(a = 0, b = 2, col = "black", lty = 2) # Add black dotted line with slope 2
        text(x = min(log10(meani)), y = max(log10(vari)) - 0.5)#, labels = paste("Slope:", round(slope_vari, 2)), pos = 4)
        legend("topright", legend = paste("p-value:", round(p_value_Vari, 5)), col = "black", cex = 0.8)
        text(x = min(log10(meani)), y = max(log10(vari)), 
             labels = paste("Slope:", round(slope_vari, 2)), pos = 4)
      } else if (form == "CVi") {
        plot(log10(meani[CV0]), log10(CVi[CV0]), 
             main = paste("Taylor's Power Law: CV vs Mean\nTime Window:", start_year, "-", end_year),
             xlab = "log10(Mean i)", ylab = "log10(CV i)",
             pch = 19, col = "blue")
        abline(lm(log10(CVi[CV0]) ~ log10(meani[CV0])), col = "red")
        text(x = min(log10(meani[CV0])), y = max(log10(CVi[CV0])) - 0.5) #, labels = paste("Slope:", round(slope_CVi, 2)), pos = 4)
        legend("topright", legend = paste("p-value:", round(p_value_CVi, 5)), col = "black", cex = 0.8)
        text(x = min(log10(meani[CV0])), y = max(log10(CVi[CV0])), 
             labels = paste("Slope:", round(slope_CVi, 2)), pos = 4)
      }
    } else {
      stop("The matrix must contain more than one column for this analysis.")
    }
    
    # Wait for user input to proceed to the next plot, handle "Esc" to exit
    tryCatch({
      readline(prompt = "Press [Enter] to see the next plot or [Esc] to exit...")
    }, interrupt = function(ex) {
      cat("\nPlotting interrupted by user.\n")
      break
    })
    
    start_year <- start_year + 1
  }
}


#Synchrony decomposition function from Zhao et al 2022
zhao_decomp <- function(X){
  
  X1 <- X[,colSums(X, na.rm=T)>0]  #remove sp not showing
  #if all sp. in this year are NAs, delete this year from tot 
  i.not.na <- !(apply(X1, 1, function(xx) all(is.na(xx))))
  X2 <- X1[i.not.na,]
  
  
  tot <- rowSums(X2,na.rm=T)
  mu <- mean(tot)
  stab.comm <- mu / sd(tot, na.rm=T)       #the reverse of CVcom
  p <- colSums(X2,na.rm=T) / sum(tot)
  alpha.inv.simpson <- 1/sum(p^2)      #Simpson
  alpha.shannon <- -sum(p * log(p))  #Shannon
  nn <- apply(X2, 1, function(xx) sum(xx>0))
  richness <- mean(nn)  #species richness
  evenness.shannon <- alpha.shannon / log(length(p))  #evenness of the distribution of density across sp.
  
  stab.sp <- mu / sum(apply(X2, 2, sd, na.rm=T))   #the reverse of CVpop
  asyn <- stab.comm / stab.sp      #based on LdM
  
  compensatory <- sqrt(sum(apply(X2, 2, var, na.rm=T)) / var(tot))
  statistical <- sum(apply(X2, 2, sd, na.rm=T), na.rm=T) / sqrt(sum(apply(X2, 2, var, na.rm=T), na.rm=T))
  
  ### asyn = compensatory * statistical
  
  
  return(c(stab.comm=stab.comm,  stab.pop=stab.sp, asyn=asyn, 
           compensatory=compensatory, statistical=statistical,
           richness=richness, shannon=alpha.shannon,
           inv.simpson=alpha.inv.simpson, evenness=evenness.shannon))
}


zhao_decomp_w <- function(df) {
  
  df_zd <- tibble()
  df_zd_yr <- tibble()
  
  cmat <- as.matrix(df %>% dplyr::select(-Year))
  
  for (window_size in 5:nrow(df)) {
    window_results <- tibble()
    
    for (start_year in 1:(nrow(df) - window_size + 1)) {
      
      
      end_year <- start_year + window_size - 1
      
      cmat_window <- cmat[start_year:end_year, ]
      
      df_zd_sub <- zhao_decomp(cmat_window)
      
      df_zd_sub$start_year <- start_year
      df_zd_sub$end_year <- end_year
      df_zd_sub$window_size <- window_size
      
      
      df_zd_yr <- rbind(df_zd_yr, df_zd_sub)
      
    }
    
    df_zd <- rbind(df_zd, df_zd_yr)
    
  }
  return(df_zd)
}


#Turn count data into density
count_to_density <- function(N, C_D, C_P){
  result <- (N/(50 * pi * 400^2)) * ((400/C_D)^2) * C_P * (1.5e7/1)
  return(result)
} 


CalcZeroInfGeomDens = function(x){
  prob_obs = sum(x > 0)/length(x)
  geom_mean_dens = ifelse(prob_obs > 0, exp(mean(log(x[x > 0]))), 0)
  return(geom_mean_dens * prob_obs)
}


