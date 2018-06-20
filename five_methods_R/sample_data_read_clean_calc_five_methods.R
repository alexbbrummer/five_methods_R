####  This code serves to clean up .tsv files for use with other R functions, specifically those needed for performing the 5 methods.  Afterward, sample code is presented for running the 5 methods, as well as plotting some of the basic results.


########## READ AND CLEAN DATA TABLE ###########

## First set the working directory to the location of your individual .tsv file.

setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/angicart_cplusplus_multiples/new_version_plus_comparison_methods/output_comparison_minSurfTests/171211/hht/")

## Check directory to examine what files are present.

dir()

## Because Angicart outputs .tsv files with variable column number, we need a fancy package for reading data into R so that we can specific the number of columns we want to read in.

install.packages("data.table")
library("data.table")

## Now read in the individual .tsv file to be analyzed.  Be sure to read in the "withRoots" versions of the .tsv files as they have the parent to child lineage informatiopn

## Also, take note of the units that you are working in as they are removed from the data table for easier integration with existing code.
treeDF <- fread(input = "hht_01_voxdim_1400x1400x1680_th_28_withRoots.tsv", col.names = c("nodeid", "vol", "length", "radius_vol", "radius", "parent", "n"), sep = "\t", select = c(1:7), fill = TRUE, na.strings = "N/A")


## All of the pre-written functions are based on diameter measurements as opposed to radius measurements, so we need to quickly convert the radii to diameter.

treeDF[ , c(4,5)] <- treeDF[ , c(4,5)]*2
names(treeDF)[c(4,5)] <- c("diameter_vol", "diameter")


## Now we can start using preexisting functions to add needed columns to the data table for further analysis.

########## BRANCHING RATIOS ###########
## This method produces one number for every vessel. thus it is a local network metric.  As such, extra columns are added to the data table upon using the function get_ratio().  See the R script get_ratio.R to examine the function itself.  Afterward are presented examples of plotting and calculating statistics (means and variance).

##Calculating and adding branching ratios beta and gamma to the data table.  We'll do this using the get_ratio() function.

beta_vol <- get_ratio(treeDF, treeDF$diameter_vol)
beta <- get_ratio(treeDF, treeDF$diameter)
gamma <- get_ratio(treeDF, treeDF$length)

## Redefine data table and add new columns
treeDF <- data.frame(treeDF, "beta_vol" = beta_vol, "beta" = beta, "gamma" = gamma)

## From here one can plot histograms of distributions of beta and gamma, and calculate scaling exponents based on the mean values of the ratios.  Sample code is provided below.

#### Scaling exponent a from beta

## First remove some bad rows and/or outliers
if(length(which(treeDF$beta == Inf)) > 0){
  bad_rows <- which(treeDF$beta == Inf)
  beta_vec <- treeDF$beta[-bad_rows]
}else{
  beta_vec <- treeDF$beta
}
## Now plot histogram, calculate mean, median, and variance in beta, and caluclate associated scaling exponent a
hist(beta_vec, breaks = ceiling(sqrt(length(beta_vec))))

mean(beta_vec, na.rm = T)
median(beta_vec, na.rm = T)
var(beta_vec, na.rm = T)

a <- -log(mean(beta_vec, na.rm = T))/log(2)
## Check this formula for calculating function of variance of a disribution.  It is not necessarily correct.
a_var <- -log(var(beta_vec, na.rm = T))/log(2)

#### Scaling exponent b from gamma

## First remove some bad rows and/or outliers
if(length(which(treeDF$gamma >= 10)) > 0){
  bad_rows <- which(treeDF$gamma >= 10)
  gamma_vec <- treeDF$gamma[-bad_rows]
}else{
  gamma_vec <- treeDF$gamma
}
## Now plot histogram, calculate mean, median, and variance in gamma, and calculate associated scaling exponent b
hist(gamma_vec, breaks = ceiling(sqrt(length(gamma_vec))))

mean(gamma_vec, na.rm = T)
median(gamma_vec, na.rm = T)
var(gamma_vec, na.rm = T)

b <- -log(median(gamma_vec, na.rm = T))/log(2)
## Check this formula for calculating function of variance of a disribution.  It is not necessarily correct.
b_var <- -log(var(gamma_vec, na.rm = T))/log(2)




########## CONSERVATAION BASED SCALING EXPONENTS ###########
## This method produces one number for every vessel that satisfies assumptions (bifurcuation), thus it is a local network metric.  As such, extra columns are added to the data table upon using the function get_conservation().  See the R script get_conservation.R to examine the function itself.
treeDF <- get_conservation(treeDF)

#### Now plot histogram, and calculate mean a.
## First remove bad rows and/or outliers.  Note that we are temporarily focusing on positive values of q, but the negative values can also be of interest.  To examine them, vary the logical statements below
if(length(which(treeDF$q >= 0 & treeDF$q <= 1000))){
  bad_rows <- which(treeDF$q >= 0 & treeDF$q <= 1000)
  q_vec <- treeDF$q[which(treeDF$q >= 0 & treeDF$q <= 1000)]
}else{
  q_vec <- treeDF$q
}
## Plotting the histogram
hist(q_vec, breaks = ceiling(sqrt(length(q_vec))))
## Calculate mean, variance, and median.
mean(q_vec)
var(q_vec)
median(q_vec)


#### Now plot histogram, and calculate mean b.
## First remove bad rows and/or outliers.  Note that we are temporarily focusing on positive values of q, but the negative values can also be of interest.  To examine them, vary the logical statements below
if(length(which(treeDF$s >= 0 & treeDF$s <= 1000))){
  bad_rows <- which(treeDF$s >= 0 & treeDF$s <= 1000)
  s_vec <- treeDF$s[which(treeDF$s >= 0 & treeDF$s <= 1000)]
}else{
  s_vec <- treeDF$s
}

## Plotting the histogram
hist(s_vec, breaks = ceiling(sqrt(length(s_vec))))
## Calculate mean, variance, and median.
mean(s_vec)
var(s_vec)
median(s_vec)




########## DISTRIBUTRION BASED SCALING EXPONENT ###########
## This method produces only one number for a collection of data/network (it is a global network metric).  Thus, all we provide here is the calling of the method to produce the exponents, uncertainty, goodness of fit, and graphs if desired.  As such, no extra columns are added to the data table.  See the .R script to examine the function itself.  Binning_type determines the whether to use linear (1) or logarithmic (2) binning.  Output is provided in the following order: a value, a R-squared fit, a Confidence Interval at 95%, b value, b R-squared fit, and b confidence interval at 95%.  Finally, by setting plot = TRUE (the default), plots of the regression fit will be automatically produced.  To remove such plots, set plot = FALSE
get_distribution_exp(treeDF = treeDF, binning_type = 1, plot = TRUE)




########## REGRESSION BASED SCALING EXPONENT ###########
## This method produces only one number for a collection of data/network (it is a global network metric).  However, this method does require knowing the number of distal tips to any given branch.  Thus, we provide here the calling of the methods to calculate the distal tips as well as produce the exponents, uncertainty, goodness of fit, and graphs if desired.  As such, no extra columns are added to the data table.  See the .R script to examine the function itself.  Output is provided in the following order: a value, a R-squared fit, a Confidence Interval at 95%, b value, b R-squared fit, and b confidence interval at 95%.  Finally, by setting plot = TRUE (the default), plots of the regression fit will be automatically produced.  To remove such plots, set plot = FALSE
treeDF <- get_distal_tips(treeDF)

get_regression_exp(treeDF = treeDF, plot = TRUE)



########## HIERARCHICAL AVERAGING BASED SCALING EXPONENT ###########
## This method produces only one number for a collection of data/network (it is a global metric).  However, this method is based on averaging over the nodal (locally calculated) branching ratios, or scale factors, beta and gamma.  The idea is that we are sorting the beta and gamma distributions based on branch size (radii for beta and lengths for gamma), then binning (either linear or logarithmic) to then average across the size classes.  This approach could in principle be applied to the conservation exponent as well, but implementation has not yet occured.  In order to generate error bars, the method of boot-strapping is employed.  It is from the boot-strapping distribution that we also extract the mean value for the hierarchically averaged result. By setting plot = TRUE (the default), plots of the boot strapped distribution statistic are generated for use in checking robustness of approach in terms of distiribution normality and continuity.  To remove such plots, set plot = FALSE.  Also, the option to toggle the number of bootstrap simulations is available through runs.  The default is 1000, which is OK for basic testing.  For presenting results internally use 5000, for sharing results use 10000, for publishing use 100000.  Output is provided in the following order: a value, a Confidence Interval at 95% lower, a confidence interavl at 95% upper, b value, b confidence interval at 95% lower, b confidence interval at 95% upper.  For goodness of method, one must inspect the resulting bootstrapped distributions for normality.  They are presented in the order of beta then gamma.  

## Still need to transform mean beta values, and 95% confidence intervals from beta space to exponent space.  Same for gamma.

get_hierarchical_exp(treeDF = treeDF, plot = TRUE, runs = 500, seed = 1234)


