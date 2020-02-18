## This function calculates the bias corrected distance correlation 
## between the extracted pseudotime and the input data. If there are 
## NAs in the data or the pseudotime vector (in case of a branching 
## trajectory) this metric cannot be calculated.

Dist_cor<-function(resultFiles){
  
  reswd = sub("/int.*","",resultFiles[1])
  wd = sub("/results","",reswd)
  metwd = paste0(reswd,"/metrics")
  dir.create(file.path(metwd), showWarnings = FALSE)
  
  distCorMatrix <- data.frame(distCor = rep(0,length(resultFiles)))
  for(i in 1:length(resultFiles)){
    
    # load the data
    input_matrix <- read.table (file = resultFiles[i])
    
    # extract parameters
    npar = which(colnames(input_matrix) == "D.timepoint") - 1
    dat <- input_matrix[,1:npar]
    Pseudotime <- input_matrix[,npar + 2]
    
  	# calculate the bias corrected distance correlation
    d <- energy::dcorT.test(Pseudotime, dat)
    distCorMatrix[i,1] <- d$estimate
  }
	return(distCorMatrix)
}