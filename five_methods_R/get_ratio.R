get_ratio <- function(treeDF, treelist){
  
  #This function returns a list of the geometric ratios in a
  #given tree data.frame.  The treeDF data.frame must contain the
  #columns "parent" and the submitted "treelist" (i.e. 
  #treeDF$diameter or treeDF$length).
  
  ratiolist <- mat.or.vec(nrow(treeDF), 1)
  for (i in 1:length(treelist)){
    ## Check that parent branch exists (is branch i a root or not?)
    if(!is.na(treeDF$parent[i])){
      ratiolist[i] <- treelist[i]/treelist[which(treeDF$nodeid==treeDF$parent[i])]
    } else{
      ratiolist[i] <- NA
    }
  }
  return(ratiolist)
}