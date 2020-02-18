### Calculate Robustness
# nruns: define how many different datasets will be contrasted
# cell.subset : how many cells should be used (we used 20% )

Robustness <-function(resultFiles, nruns = 100, cell.subset = 0.2){

  reswd = sub("/int.*","",resultFiles[1])
  wd = sub("/results","",reswd)
  metwd = paste0(reswd,"/metrics")
  dir.create(file.path(metwd), showWarnings = FALSE)
  
  Robustness_result <- data.frame()
  for (f in 1:length(resultFiles)){
    
    # load the data
    input_matrix <- read.table (file = resultFiles[f])
    
    npar = which(colnames(input_matrix) == "D.timepoint") - 1
    ncells = nrow(input_matrix)
    expr_matrix = input_matrix[,1:npar]
    timepoint.annot = input_matrix[, npar + 1]

    #subset_order will have the order of each cell, after the TI method is applied on each subset.
    nsub.cells <- round(cell.subset*ncells)
    initial_order <- subset_order <- as.data.frame(matrix(nrow =  nsub.cells, ncol = nruns))
    
    ####
    #### ~ ~ ~ apply alg on subsampled data ~ ~ ~
    ####
    
    # first decode the TI method used (cr: curve reconstruction, dr: dimensionality reduction, act: activation, transf: transformation)
    experiment = sub(".*results/","",sub("_output.*","",resultFiles[f]))
    experiment = unlist(strsplit(experiment,"_"))
    act = experiment[2]
    transf = experiment[4]
    cr.method = experiment[5]
    dr.method = experiment[6]
    
    lineages <- sum(sapply("curve", grepl, colnames(input_matrix))) # if cr.method == slingshot
    if (lineages > 1){
        spear <- 0
        kend <- 0
    } else {
      FUN = paste0("run_",cr.method)
      #setup parallel backend to use many processors
      cores=detectCores()
      cl <- makeCluster(cores[1]-1) #do not overload the computer
      registerDoParallel(cl)
      
      finalMatrix <- foreach(j=1:nruns, .export = FUN, .verbose = T) %dopar% {
        # Subsampling from the initial dataset, with stratification
        subIndex <- sample(1:ncells, nsub.cells)
        #print(subIndex[1:10])
        subMatrix <- list(expr = expr_matrix[subIndex, ], 
                          ntimepoints = length(unique(input_matrix$D.timepoint)),
                          timepoint = timepoint.annot[subIndex])
        
        tempMatrix = do.call(FUN, list(subMatrix, act, transf, dr.method, randomStart = T, write2disk = F))
        tempMatrix[["subIndex"]] = subIndex
        
        tempMatrix
      }
      #stop cluster
      stopCluster(cl)
      
      spear <- kend <- c()
      # # Compare the initial ordering with those returned by the subsets
      for (n in 1:nruns){
        if (finalMatrix[[n]]$lineages == 1){
          initial_order[,n] <- input_matrix$pseudot[finalMatrix[[n]]$subIndex]
          subset_order[,n] <- finalMatrix[[n]]$pseudotime
        
          spear <- c(spear, cor(initial_order[,n],subset_order[,n], use="pairwise.complete.obs", method = "spearman"))
          kend <- c(kend, cor(initial_order[,n],subset_order[,n], use="pairwise.complete.obs", method = "kendall"))
        } 
      }
      if (is.null(spear)) { spear <- 0 }
      if (is.null(kend)) { kend <- 0 }
    }

    Robustness_result <- rbind(Robustness_result, c(mean(abs(spear)), mean(abs(kend))))
    rownames(Robustness_result)[f] <- paste0(cr.method,"_",dr.method)
  }
  
  colnames(Robustness_result) <- c("spearman", "kendall")
  return(Robustness_result)  
}


