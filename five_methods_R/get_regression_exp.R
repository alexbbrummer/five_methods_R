get_regression_exp <- function(treeDF, plot = TRUE){
  ## This function calculates the scaling exponents a and b using the regression method.  Output is provided in the following order: a value, a R-squared fit, a Confidence Interval at 95%, b value, b R-squared fir, and b confidence interval at 95%.  Finally, by setting plot = TRUE (the default), plots of the regression fit will be automatically produced.  To remove such plots, set plot = FALSE
  
  ## Check that distal tips have been calculated.
  if(!("distal_tips" %in% names(treeDF))){
    treeDF <- get_distal_tips(treeDF)
  }
  
  ## Need to identify bad rows.  These are rows that are either terminal branches (distal_tips = 0), or where the diameter or length are equal to 0.
  bad_rows <- unique(c(which(treeDF$distal_tips == 0), which(treeDF$diameter == 0), which(treeDF$length == 0)))
  
  library("lmodel2")
  
  log_tips <- log(treeDF$distal_tips[-bad_rows])
  log_diameter <- log(treeDF$diameter[-bad_rows])
  log_length <- log(treeDF$length[-bad_rows])
  
  ## Calculate a
  tryCatch({
    
    treeDF_lm <- lmodel2(log_diameter~log_tips) 
    
    a_rsq <- treeDF_lm$rsquare
    a_slp <- treeDF_lm$regression.results[3,3]
    a_ci <- treeDF_lm$confidence.intervals[3,5] - treeDF_lm$confidence.intervals[3,4]
    a_int <- treeDF_lm$regression.results[3,2]
    
    rp = vector('expression', 3)
    rp[1] = substitute(expression(a == SLP),
                       list(SLP = format(a_slp, dig = 3)))[2]
    rp[2] = substitute(expression('95%' == CI),
                       list(CI = format(a_ci, dig = 3)))[2]
    rp[3] = substitute(expression(italic(R)^2 == RSQ),
                       list(RSQ = format(a_rsq, dig = 3)))[2]
    
    ## Comment the call to png() and dev.off() to plot in RStudio.  Otherwise, check your directory is appropriate.
    plot(log_tips, log_diameter, main = NULL, xlab = "", ylab = "", xlim = c(-1, max(log_tips)))
    title(xlab = expression(paste("ln(Number of Distal Terminal Branches)")), ylab = expression(paste("ln(Diameter)")), mgp = c(2.25, 1, 0))
    abline(a = a_int, b = a_slp, col = "red")
    legend("topleft", legend = rp, bty = "n")
    
  }, error = function(err) {
    
    print(paste("Error occurred, returning NA's", err))
    a_rsq <- NA
    a_slp <- NA
    a_ci <- NA
    a_int <- NA
  })
  
  ## Calculate b
  tryCatch({
    
    treeDF_lm <- lmodel2(log_length~log_tips) 
    
    b_rsq <- treeDF_lm$rsquare
    b_slp <- treeDF_lm$regression.results[3,3]
    b_ci <- treeDF_lm$confidence.intervals[3,5] - treeDF_lm$confidence.intervals[3,4]
    b_int <- treeDF_lm$regression.results[3,2]
    
    rp = vector('expression', 3)
    rp[1] = substitute(expression(b == SLP),
                       list(SLP = format(b_slp, dig = 3)))[2]
    rp[2] = substitute(expression('95%' == CI),
                       list(CI = format(b_ci, dig = 3)))[2]
    rp[3] = substitute(expression(italic(R)^2 == RSQ),
                       list(RSQ = format(b_rsq, dig = 3)))[2]
    
    ## Comment the call to png() and dev.off() to plot in RStudio.  Otherwise, check your directory is appropriate.
    plot(log_tips, log_length, main = NULL, xlab = "", ylab = "", xlim = c(-1, max(log_tips)))
    title(xlab = expression(paste("ln(Number of Distal Terminal Branches)")), ylab = expression(paste("ln(Length)")), mgp = c(2.25, 1, 0))
    abline(a = b_int, b = b_slp, col = "red")
    legend("topleft", legend = rp, bty = "n")
    
  }, error = function(err) {
    
    print(paste("Error occurred, returning NA's", err))
    a_rsq <- NA
    a_slp <- NA
    a_ci <- NA
    a_int <- NA
  })
return(c(a_slp, a_rsq, a_ci, b_slp, b_rsq, b_ci))
}