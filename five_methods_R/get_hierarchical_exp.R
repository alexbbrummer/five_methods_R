get_hierarchical_exp <- function(treeDF, runs = 1000, plot = TRUE, seed = 1234){
  ## This method produces only one number for a collection of data/network (it is a global metric).  However, this method is based on averaging over the nodal (locally calculated) branching ratios, or scale factors, beta and gamma.  The idea is that we are sorting the beta and gamma distributions based on branch size (radii for beta and lengths for gamma), then binning (either linear or logarithmic) to then average across the size classes.  This approach could in principle be applied to the conservation exponent as well, but implementation has not yet occured.  In order to generate error bars, the method of boot-strapping is employed.  It is from the boot-strapping distribution that we also extract the mean value for the hierarchically averaged result. By setting plot = TRUE (the default), plots of the boot strapped distribution statistic are generated for use in checking robustness of approach in terms of distiribution normality and continuity.  To remove such plots, set plot = FALSE.  Also, the option to toggle the number of bootstrap simulations is available through runs.  The default is 1000, which is OK for basic testing.  For presenting results internally use 5000, for sharing results use 10000, for publishing use 100000.
  
  
  ## Load library "mltools" to use the bin_data() function later.
  library("mltools")
  
  ## Remove bad_rows as defined by rows where beta or gamma equal infinity.
  bad_rows <- unique(c(which(treeDF$beta == Inf), which(treeDF$gamma == Inf)))
  treeDF <- treeDF[-bad_rows,]
  
  ## First we must define the function that we will be bootstrapping with later.
  hierarchical_ave <- function(treeDF, d){
    ## This might be the most mysterious line present.  The boot strapping function requires a variable index for future use, ence the presence of the undefined "d" index.
    treeDF_sort <- treeDF[d,]
    
    # Defining number of bins
    dia_bin_number <- ceiling(sqrt(length(treeDF_sort$diameter)))
    len_bin_number <- ceiling(sqrt(length(treeDF_sort$length)))
    
    # Define temporary data.frame for sorting/binning.  Here we sue the bin_data() function to bin the diameter and length bins by the number of bins defined above, and identify which bin each diameter or length variable sits inside of.
    treeDF_sort <- data.frame(treeDF_sort, diameter_bins = bin_data(treeDF_sort$diameter, bins = dia_bin_number), length_bins = bin_data(treeDF_sort$length, bins = len_bin_number))
    
    # Initialize vectors for bin means
    beta_bin_means <- c()
    gamma_bin_means <- c()
    
    # Loop through diameter and length bins to search for corresponding scale factors and caluclate means within each bin
    for(i in 1:length(levels(treeDF_sort$diameter_bins))){
      dia_bin_rows <- which(treeDF_sort$diameter_bins == levels(treeDF_sort$diameter_bins)[i])
      beta_bin_means[i] <- mean(treeDF_sort$beta[dia_bin_rows], na.rm = T)
    }
    
    for(i in 1:length(levels(treeDF_sort$length_bins))){
      len_bin_rows <- which(treeDF_sort$length_bins == levels(treeDF_sort$length_bins)[i])
      gamma_bin_means[i] <- mean(treeDF_sort$gamma[len_bin_rows], na.rm = T)
    }
    # Average over bin means.
    beta_mean <- mean(beta_bin_means, na.rm = T)
    gamma_mean <- mean(gamma_bin_means, na.rm = T)
    
    return(c(beta_mean, gamma_mean))
  }
  
  
  # This section runs the bootstrapping
  library("boot")
  
  # For reproducibility, set seed for random number generator
  set.seed(seed)
  
  # Perform bootsampling.
  boot_test <- boot(data = treeDF, statistic = hierarchical_ave, R = runs)
  
  
  # Estimate bootsample mean of beta
  beta_boot <- mean(boot_test[[2]][,1])
  
  # Estimate bootsample mean of gamma
  gamma_boot <- mean(boot_test[[2]][,2])
  
  # Estimate upper and lower confidence levels for 95% of beta
  beta_boot_err_low <- boot.ci(boot_test, type = "perc", index = c(1))[[4]][4]
  beta_boot_err_high <- boot.ci(boot_test, type = "perc", index = c(1))[[4]][5]
  
  # Estimate upper and lower confidence levels for 95% of b
  gamma_boot_err_low <- boot.ci(boot_test, type = "perc", index = c(2))[[4]][4]
  gamma_boot_err_high <- boot.ci(boot_test, type = "perc", index = c(2))[[4]][5]
  
  
  if(plot){
    plot(boot_test, index = c(1))
    plot(boot_test, index = c(2))
    # hist(boot_test[[2]][,1], breaks = ceiling(sqrt(nrow(boot_test[[2]]))), xlab = "", main = "")
    # title(main = sprintf("Histogram of bootsampled beta from hierarchical averaging \n R = %i", runs), xlab = "Bootsampled beta")
    # 
    # hist(boot_test[[2]][,2], breaks = ceiling(sqrt(nrow(boot_test[[2]]))), xlab = "", main = "")
    # title(main = sprintf("Histogram of bootsampled gamma from hierarchical averaging \n R = %i", runs), xlab = "Bootsampled gamma")
  }
  return(c(beta_boot, beta_boot_err_low, beta_boot_err_high, gamma_boot, gamma_boot_err_low, gamma_boot_err_high))
}