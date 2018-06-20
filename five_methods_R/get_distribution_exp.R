get_distribution_exp <- function(treeDF, binning_type, plot = TRUE){
  ## This function calculates the scaling exponents a and b using the distribution method.  Binning_type determines the whether to use linear (1) or logarithmic (2) binning.  Output is provided in the following order: a value, a R-squared fit, a Confidence Interval at 95%, b value, b R-squared fir, and b confidence interval at 95%.  Finally, by setting plot = TRUE (the default), plots of the regression fit will be automatically produced.  To remove such plots, set plot = FALSE
  
  library("lmodel2")
  
  ## Bin size determination
  
  n_bin_dia <- ceiling(sqrt(length(treeDF$diameter)))
  
  n_bin_len <- ceiling(sqrt(length(treeDF$length)))
  
  if(binning_type == 1){
    ## Linear binning
    # Calculate a
    
    a_hist <- hist(treeDF$diameter, breaks = n_bin_dia, plot = FALSE)
    empty_bins <- which(a_hist$counts == 0)
    counts <- a_hist$counts[-empty_bins]
    centers <- a_hist$mids[-empty_bins]
    
    relativefreq <- counts/length(treeDF$diameter)
    
    if(length(centers) <= 1){
      a_dist <- NaN
    }else{
      
      x_log <- log(centers)
      y_log <- log(relativefreq)
      
      tryCatch({
        treeDF_lm <- lmodel2(y_log~x_log) 
        
        a_rsq <- treeDF_lm$rsquare
        a_slp <- treeDF_lm$regression.results[3,3]
        a_ci <- treeDF_lm$confidence.intervals[3,5] - treeDF_lm$confidence.intervals[3,4]
        a_int <- treeDF_lm$regression.results[3,2]
        a_dist <- -1/a_slp
        
        if(plot == TRUE){
          # Format nicely for legend
          rp = vector('expression', 3)
          rp[1] = substitute(expression(a == SLP),
                             list(SLP = format(a_dist, dig = 3)))[2]
          rp[2] = substitute(expression('95%' == CI),
                             list(CI = format(a_ci, dig = 3)))[2]
          rp[3] = substitute(expression(italic(R)^2 == RSQ),
                             list(RSQ = format(a_rsq, dig = 3)))[2]
          # Make graph
          plot(x_log, y_log, main = NULL, xlab = "", ylab = "")
          title(ylab= expression(paste("ln(Relative Freq.)")), xlab = expression(paste("ln(Radius (mm))")), mgp = c(2.25, 1, 0))
          abline(a = a_int, b = a_slp, col = "red")
          legend("topright", legend = rp, bty = "n")
        }
        
      }, error = function(err) {
        
        print(paste("Error occurred, returning NA's", err))
        a_rsq <- NA
        a_slp <- NA
        a_ci <- NA
        a_int <- NA
        a_dist <- NA
      })
    }
      
      # Calculate b
      
      b_hist <- hist(treeDF$length, breaks = n_bin_dia, plot = FALSE)
      empty_bins <- which(b_hist$counts == 0)
      counts <- b_hist$counts[-empty_bins]
      centers <- b_hist$mids[-empty_bins]
      
      relativefreq <- counts/length(treeDF$length)
      
      if(length(centers) <= 1){
        b_dist <- NaN
      }else{
        
        x_log <- log(centers)
        y_log <- log(relativefreq)
        
        tryCatch({
          treeDF_lm <- lmodel2(y_log~x_log) 
          
          b_rsq <- treeDF_lm$rsquare
          b_slp <- treeDF_lm$regression.results[3,3]
          b_ci <- treeDF_lm$confidence.intervals[3,5] - treeDF_lm$confidence.intervals[3,4]
          b_int <- treeDF_lm$regression.results[3,2]
          b_dist <- -1/b_slp
          
          if(plot == TRUE){
            # Format nicely for legend
            rp = vector('expression', 3)
            rp[1] = substitute(expression(b == SLP),
                               list(SLP = format(b_dist, dig = 3)))[2]
            rp[2] = substitute(expression('95%' == CI),
                               list(CI = format(b_ci, dig = 3)))[2]
            rp[3] = substitute(expression(italic(R)^2 == RSQ),
                               list(RSQ = format(b_rsq, dig = 3)))[2]
            # Make graph
            plot(x_log, y_log, main = NULL, xlab = "", ylab = "")
            title(ylab= expression(paste("ln(Relative Freq.)")), xlab = expression(paste("ln(Length)")), mgp = c(2.25, 1, 0))
            abline(a = b_int, b = b_slp, col = "red")
            legend("topright", legend = rp, bty = "n")
          }
          
        }, error = function(err) {
          
          print(paste("Error occurred, returning NA's", err))
          b_rsq <- NA
          b_slp <- NA
          b_ci <- NA
          b_int <- NA
          b_dist <- NA
        })
      }
  }else if(binning_type == 2){
    ## Transform radius and length to log scale
    dia_log <- log(treeDF$diameter)
    len_log <- log(treeDF$length)
    
    # Calculate a
    
    a_hist <- hist(dia_log, breaks = n_bin_dia, plot = FALSE)
    empty_bins <- which(a_hist$counts == 0)
    counts <- a_hist$counts[-empty_bins]
    centers <- exp(a_hist$mids[-empty_bins])
    
    relativefreq <- counts/length(dia_log)
    
    if(length(centers) <= 1){
      a_dist <- NaN
    }else{
      
      x_log <- log(centers)
      y_log <- log(relativefreq)
      
      tryCatch({
        treeDF_lm <- lmodel2(y_log~x_log) 
        
        a_rsq <- treeDF_lm$rsquare
        a_slp <- treeDF_lm$regression.results[3,3]
        a_ci <- treeDF_lm$confidence.intervals[3,5] - treeDF_lm$confidence.intervals[3,4]
        a_int <- treeDF_lm$regression.results[3,2]
        a_dist <- -1/a_slp
        
        if(plot == TRUE){
          # Format nicely for legend
          rp = vector('expression', 3)
          rp[1] = substitute(expression(a == SLP),
                             list(SLP = format(a_dist, dig = 3)))[2]
          rp[2] = substitute(expression('95%' == CI),
                             list(CI = format(a_ci, dig = 3)))[2]
          rp[3] = substitute(expression(italic(R)^2 == RSQ),
                             list(RSQ = format(a_rsq, dig = 3)))[2]
          # Make graph
          plot(x_log, y_log, main = NULL, xlab = "", ylab = "")
          title(ylab= expression(paste("ln(Relative Freq.)")), xlab = expression(paste("ln(Radius)")), mgp = c(2.25, 1, 0))
          abline(a = a_int, b = a_slp, col = "red")
          legend("topright", legend = rp, bty = "n")
        }
        
        
      }, error = function(err) {
        
        print(paste("Error occurred, returning NA's", err))
        a_rsq <- NA
        a_slp <- NA
        a_ci <- NA
        a_int <- NA
        a_dist <- NA
      })
    }
    
    # Calculate b
    
    b_hist <- hist(len_log, breaks = n_bin_dia, plot = FALSE)
    empty_bins <- which(b_hist$counts == 0)
    counts <- b_hist$counts[-empty_bins]
    centers <- exp(b_hist$mids[-empty_bins])
    
    relativefreq <- counts/length(len_log)
    
    if(length(centers) <= 1){
      b_dist <- NaN
    }else{
      
      x_log <- log(centers)
      y_log <- log(relativefreq)
      
      tryCatch({
        treeDF_lm <- lmodel2(y_log~x_log) 
        
        b_rsq <- treeDF_lm$rsquare
        b_slp <- treeDF_lm$regression.results[3,3]
        b_ci <- treeDF_lm$confidence.intervals[3,5] - treeDF_lm$confidence.intervals[3,4]
        b_int <- treeDF_lm$regression.results[3,2]
        b_dist <- -1/b_slp
        
        
        if(plot == TRUE){
          # Format nicely for legend
          rp = vector('expression', 3)
          rp[1] = substitute(expression(b == SLP),
                             list(SLP = format(b_dist, dig = 3)))[2]
          rp[2] = substitute(expression('95%' == CI),
                             list(CI = format(b_ci, dig = 3)))[2]
          rp[3] = substitute(expression(italic(R)^2 == RSQ),
                             list(RSQ = format(b_rsq, dig = 3)))[2]
          # Make graph
          plot(x_log, y_log, main = NULL, xlab = "", ylab = "")
          title(ylab= expression(paste("ln(Relative Freq.)")), xlab = expression(paste("ln(Length)")), mgp = c(2.25, 1, 0))
          abline(a = b_int, b = b_slp, col = "red")
          legend("topright", legend = rp, bty = "n")
        }
        
      }, error = function(err) {
        
        print(paste("Error occurred, returning NA's", err))
        b_rsq <- NA
        b_slp <- NA
        b_ci <- NA
        b_int <- NA
        b_dist <- NA
      })
    }
  }
  return(c(a_dist, a_rsq, a_ci, b_dist, b_rsq, b_ci))
}