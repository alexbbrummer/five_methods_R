get_conservation <- function(treeDF){
  ### This function calculates the radial and length scalin exponents a and b from the conservation based approach.  The method employed is to identify roots of the conservation based equations.  An alternative method is possible where the distribution of sibling scale factors are fit with a non-linear fitting algorithm.  This approach has not been tested extensively and is worth further investigation at a later time.
  
  ## First we must define the two equations we are going to be solving.
  
  betafun <- function(x, beta1, beta2){
    return((beta1)**(x) + (beta2)**(x) - 1)
  }
  
  gammafun <- function(y, gamma1, gamma2){
    return((gamma1)**(y) + (gamma2)**(y) - 1)
  }
  
  ## Initialize exponent vectors.
  radial_exp <- mat.or.vec(length(treeDF$beta), 1)
  length_exp <- mat.or.vec(length(treeDF$gamma), 1)
  
  ## Loop through data table.
  for(i in 1:nrow(treeDF)){
    # Check for bifurcation
    if(treeDF$n[i] == 2){
      # Build subtable of siblings that are child to this branch
      subtab <- treeDF[which(treeDF$parent == treeDF$nodeid[i]), ]
      # Define beta1, beta2, gamma1, and gamma2 from subtable for use in root finder below
      
      beta1 <- subtab$beta[1]
      beta2 <- subtab$beta[2]
      gamma1 <- subtab$gamma[1]
      gamma2 <- subtab$gamma[2]
      
    # Now we can use the root finters to calculate the exponents.
      radial_exp[i] <- tryCatch(uniroot(betafun, lower = -1000, upper = 1000, beta1 = beta1, beta2 = beta2)$root,
                                error = function(error_message){
                                  return(NaN)
                                })
      length_exp[i] <- tryCatch(uniroot(gammafun, lower = -1000, upper = 1000, gamma1 = gamma1, gamma2 = gamma2)$root,
                                error = function(error_message){
                                  return(NaN)
                                })
    }else{
      radial_exp[i] <- NA
      length_exp[i] <- NA
    }
  }
  # Check for which exponents were defaulted to maximal values.  These are asymptotically erroneous values to be thrown out.
  radial_exp[which(abs(radial_exp) >= 1000)] <- NaN
  length_exp[which(abs(length_exp) >= 1000)] <- NaN
  return(data.frame(treeDF, "q" = 1/radial_exp, "s" = 1/length_exp))
}
