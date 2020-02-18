run_prinCurves <- function(D, activator, transf, drMethod, randomStart = F, write2disk = T, wd = NULL){
  
  # 1st step => dim reduction
  dr_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_",drMethod,".txt") 
  
  # if a file with the requested dimensionality reduction analysis exists, do not repeat 
  if ((file.exists(paste0("results/",dr_file_name)) == FALSE) | (randomStart == T)){
    if (drMethod == "tSNE"){ 
      set.seed(24)
      tsne_out <- Rtsne.multicore::Rtsne.multicore(D$expr, dims = 2)
      data_drMethod <- tsne_out$Y
    }
    if (drMethod == "diffMaps"){
      DifMap_output <- destiny::DiffusionMap(D$expr)
      data_drMethod <-data.frame(cbind(DifMap_output$DC1, DifMap_output$DC2))
    } 
    if (write2disk) {
      write.table(data_drMethod, file = file.path(wd, "results",dr_file_name))
    }
  } else { 
    data_drMethod <- read.table(file = paste0("results/",dr_file_name))
  }
  
  # 2nd step => trajectory
  ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_prinCurves_",drMethod,"_output.txt")  
  
  if ((file.exists(paste0("results/",ti_file_name)) == FALSE) | (randomStart == T)){
    # if we want to assign a starting position, a minimal trajectory must be set (two points)
    if (randomStart) { 
      startTraj = t(as.matrix(data_drMethod[sample(1:nrow(data_drMethod),2),]))
    } else {
      startTraj = NULL
    }
    fit1 <- princurve::principal_curve(x = as.matrix(data_drMethod), 
                                       start = startTraj,
                                       maxit = 100,  
                                       smoother = "smooth_spline")
    pseudot <- fit1$lambda
    # 3rd step => Save a table with the Abundances, Experimental time, Pseudotime values and Pseudotime trajectory
    if (write2disk) {
      algs_matrix = cbind(D$expr, D$timepoint, pseudot, fit1$s)
      write.table(algs_matrix, file = paste0("results/",ti_file_name))
      print(paste0("Finished with Principal Curves with ", drMethod, " and created", ti_file_name ,"  :) \n " ))
    } else {
      to_return <- list()
      to_return[["dimRed"]] <- data_drMethod
      to_return[["pseudotime"]] <- pseudot
      to_return[["trajectory"]] <- fit1$s
      to_return[["lineages"]] <- 1
      
      return(to_return)
    }
    
  }

}