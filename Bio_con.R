### Calculate Biological Consistency

Bio_con <- function(D, Pathway_Hierarchy_file, resultFiles, nruns = 100){

  reswd = sub("/int.*","",resultFiles[1])
  wd = sub("/results","",reswd)
  metwd = paste0(reswd,"/metrics")
  dir.create(file.path(metwd), showWarnings = FALSE)
  
  # Pathway_Hierarchy_file : the file with the pathway info. Put with path, extension
	Pathway_Hierarchy <- as.data.frame(read_delim(file = Pathway_Hierarchy_file, ";", escape_double = FALSE, trim_ws = TRUE))
	
	npar <- ncol(D$expr)
	
	method_score = list()
	
	method_score[["naive"]] <- check_pairs(dat = D$expr, t = D$timepoint, known_pairs = Pathway_Hierarchy)

  for(i in 1:length(resultFiles)){	
    
    # decode the TI method used (cr: curve reconstruction, dr: dimensionality reduction, act: activation, transf: transformation)
    experiment = sub(".*results/","",sub("_output.*","",resultFiles[i]))
    experiment = unlist(strsplit(experiment,"_"))
    act = experiment[2]
    transf = experiment[4]
    cr.method = experiment[5]
    dr.method = experiment[6]
    
    FUN = paste0("run_",cr.method)
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #do not overload the computer
    registerDoParallel(cl)
    
    finalMatrix <- foreach(j=1:nruns, .export = FUN, .verbose = T) %dopar% {
      tempMatrix = do.call(FUN, list(D, act, transf, dr.method, randomStart = T, write2disk = F))
      tempMatrix
    }
    #stop cluster
    stopCluster(cl)
    
    temp_score = rep(0,length(finalMatrix))
    for (k in 1:length(finalMatrix)){

    	# sort the plot values in ascending pseudotime order
    	ix = order(finalMatrix[[k]]$pseudotime)
    	dat = D$expr[ix,]
    	Pseudotime = finalMatrix[[k]]$pseudotime[ix]
    	if (finalMatrix[[k]]$lineages > 1 ){
    	  print(paste0("The algorithm ", cr.method," with ", dr.method ," returned a branching trajectory. Cannot calculate the metric !"))
    	  next
    	}
    	
    	# rescale pseudotime to [0,1] for the figures between algorithms to be comparable
    	Pseudotime_zscore <- (Pseudotime-min(Pseudotime,na.rm=T))/(max(Pseudotime,na.rm=T)-min(Pseudotime,na.rm=T))
    	
    	# check also the reverse direction
    	Pseudotime_zscore_rev <- reverse_pseudotime(Pseudotime_zscore)
    	
    	temp = check_pairs(dat,Pseudotime_zscore, known_pairs = Pathway_Hierarchy)
    	temp_rev = check_pairs(dat,Pseudotime_zscore_rev, known_pairs = Pathway_Hierarchy)
    	
    	temp_score[k] <- max(temp,temp_rev)
    	
    	# validate some intermediate results
    	# plot(finalMatrix[[k]]$dimRed)
    	# points(finalMatrix[[k]]$trajectory, col = "red")
    	# plot(Pseudotime_zscore_rev, dat$pPlcg2, col = D$timepoint[ix])
    	# spline.df = calc_spline(dat,Pseudotime_zscore)
    	# points(spline.df$time, spline.df$pPlcg2, col = "black")
    }
    method_score[[paste0(cr.method,"_",dr.method)]] = c(score = max(temp_score), output = finalMatrix[[which.max(temp_score)]])
  }
	
	return(method_score)
}
