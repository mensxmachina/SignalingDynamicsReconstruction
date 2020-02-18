run_pseudogp <- function(D, activator, transf, drMethod, randomStart = F, write2disk = T, wd = NULL){
  
  # 1st step => dim reduction
  dr_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_",drMethod,".txt") 
  
  # if a file with the requested dimensionality reduction analysis exists, do not repeat 
  if (file.exists(paste0("results/",dr_file_name)) == FALSE){
    if (drMethod == "tSNE"){ 
      set.seed(24)
      tsne_out <- Rtsne.multicore::Rtsne.multicore(D$expr, dims = 2)
      data_drMethod <- tsne_out$Y
      write.table(data_drMethod, file = file.path(wd, "results",dr_file_name))
      rm(tsne_out)
    }
    if (drMethod == "PCA"){
      pca <- prcomp(D$expr, scale. = FALSE)
      data_drMethod <- pca$x[,1:2]
      write.table(data_drMethod, file = file.path(wd, "results",dr_file_name))
      rm(pca)
    } 
  } else { 
    data_drMethod <- read.table(file = paste0("results/",dr_file_name))
  }
  
  # 2nd step => trajectory inference
  ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_pseudoGP_",drMethod,"_output.txt")
  
  # if a file with the inferred trajectory is already created, do not repeat
  if ((file.exists(paste0("results/",ti_file_name)) == FALSE) | (randomStart == T)){
    #  It acts on a reduced dimension representation of the data (two-dimensional) in which it fits a probabilistic curve, 
    # allowing posterior pseudotime uncertainty to be quantified. => aka again you apply any kind of dim reduction!
    # subsample <- c(sample(1:(dim(data_norm))[1], 20, replace=F))
    data_drMethod <- as.matrix(data_drMethod)
    # set.seed(24)
    s <- sample(1:sum(D$N),500)
    le_fit <- pseudogp::fitPseudotime(data_drMethod[s,], smoothing_alpha = 10, smoothing_beta = 3, iter = 100, chains = 1)
    pst <- rstan::extract(le_fit, pars="t")$t
    pseudot <- pst[dim(pst)[1],]
    
    # 3rd step => Save a table with Normalized_values [1:limit], Experimental time, Pseudotime
    algs_matrix = cbind(D$expr[s,], D$timepoint[s], pseudot)
    
    write.table(algs_matrix, file = paste0("results/",ti_file_name))
    print(paste0("Finished with", drMethod, " and created", ti_file_name ,"  :) \n " ))
  } 
  
  if ( (file.exists(dr_file_name) == T) & (file.exists(ti_file_name) == T) ){ 
    print(paste0(alg," already applied :) "))
  }
  
} # END pseudoGP_with_tSNE/PCA