get_distal_tips <- function(treeDF){
  ## This function calculates the number of distal tips to a given branch in a vascular network.  No assumptions regarding branchin type are made.
  ## Intialize distal_tips vector and insert into treeDF.
  distal_tips <- mat.or.vec(length(treeDF$n),1)
  treeDF <- data.frame(treeDF, "distal_tips" = distal_tips)
  ## Must break treeDF into subtrees by roots.
  roots <- which(is.na(treeDF$parent))
  
  for(j in 1:(length(roots)-1)){
    subtree <- treeDF[roots[j]:(roots[j+1]-1),]
    
    ## Identify tips in subtree
    if(length(which(is.na(subtree$n))) > 0){
      tips <- which(is.na(subtree$n))
    } else if(length(which(is.na(subtree$n))) == 0){
      tips <- which(subtree$n == 0)
    }
    ## Seed tips as having distal_tip values of 1
    subtree$distal_tips[tips] <- 1
    while(subtree$distal_tips[1] == 0){
      new_tips <- c()
      for(i in tips){
        sibs <- which(subtree$parent == subtree$parent[i])
        parent <- which(subtree$nodeid == subtree$parent[i])
        if((subtree$distal_tips[parent] == 0) & (min(subtree$distal_tips[sibs]) > 0)){
          subtree$distal_tips[parent] <- sum(subtree$distal_tips[sibs])
          new_tips <- c(new_tips, parent)
        } else if(min(subtree$distal_tips[sibs]) == 0){
          new_tips <- c(new_tips, i)
        }
      }
      tips <- new_tips
    }
    treeDF$distal_tips[roots[j]:(roots[j+1] - 1)] <- subtree$distal_tips
  }
  

  # treeDF$distal_tips[which(is.na(treeDF$n))] <- 0
  # treeDF$distal_tips[which(treeDF$n == 0)] <- 0
  return(treeDF)
}