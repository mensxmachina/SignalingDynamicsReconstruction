run_scorpius <- function(D, activator, transf, drMethod, randomStart = F, write2disk = T, wd = NULL){
  
# 1st step => dim reduction
dr_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_",drMethod,".txt")

  # if a file with this dim red was already created, dont recreate
if ((file.exists(paste0("results/",dr_file_name)) == FALSE) | (randomStart == T)){
  if (drMethod == "distPear"){ 
    space <- SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "pearson")
  }
  if (drMethod == "distSpear"){ 
    space <- SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "spearman")
  }
  if (drMethod == "distCos"){ 
    space <- SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "cosine")
  }
  if (drMethod == "distEucl"){ 
    space <- SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "euclidean")
  }
  if (drMethod == "distManh"){ 
    space <- SCORPIUS::reduce_dimensionality(as.matrix(D$expr), ndim = 2, dist = "manhattan")
  }
  if (write2disk) {
    write.table(space, file = file.path(wd, "results",dr_file_name))
  }
} else {space <- read.table(file = paste0("results/",dr_file_name))}

# 2nd step => trajectory 
ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_scorpius_",drMethod,"_output.txt")

# if a file with this analysis is already created, dont reapply 
if ((file.exists(paste0("results/",ti_file_name)) == FALSE) | (randomStart == T)){
  
  traj <- SCORPIUS::infer_trajectory(space)
  pseudot<- traj$time # this is your pseudot
  # 3rd step => Save a table with the results
  if (write2disk) {
    algs_matrix = cbind(D$expr, D$timepoint, pseudot,traj$path)
    write.table(algs_matrix, file = paste0("results/",ti_file_name))
    print(paste0("Finished with SCORPIUS with, ", drMethod," and created ", ti_file_name ,"  :) \n " ))
  } else {
    to_return <- list()
    to_return[["dimRed"]] <- space
    to_return[["pseudotime"]] <- pseudot
    to_return[["trajectory"]] <- traj$path
    to_return[["lineages"]] <- 1
    
    return(to_return)
  }
}

} 