### Calculate Roughness

Rough <- function(resultFiles){
  
  reswd = sub("/int.*","",resultFiles[1])
  wd = sub("/results","",reswd)
  metwd = paste0(reswd,"/metrics")
  dir.create(file.path(metwd), showWarnings = FALSE)
  
  roughMatrix <- data.frame(roughness = rep(0,length(resultFiles)), pvalue = rep(0,length(resultFiles)))
  for(i in 1:length(resultFiles)){
    
    # load the data
    input_matrix <- read.table (file = resultFiles[i])
    
    # extract parameters
    npar = which(colnames(input_matrix) == "D.timepoint") - 1
    dat <- input_matrix[,1:npar]
    Timepoints = input_matrix[, npar + 1]
    Pseudotime <- input_matrix[,npar + 2]
	

  	AlgRoughMatrix <- data.frame(matrix(NA, nrow = npar, ncol = 1))
  	r_perm_ofprotein <-  data.frame(matrix(NA, nrow = 100, ncol = 1))
  	observedRoughness<-c()
  	listProteinRoughs<-rep(0,npar)
  
  	myOrder<-c(rank(Pseudotime, ties.method = "first"))
  	start_time <-Sys.time()
  	for (my_protein in 1:npar){
  	  p = dat[myOrder, my_protein]
  	  a<-p[myOrder]
  	  a<-a[!is.na(a)]
  	  observedRoughness[my_protein] <- DeLorean:::calc.roughness(a)
  	  for (perms in 1:100){
  	    a_perm <- sample(a) 
  	    r_perm_ofprotein[perms, my_protein] <- DeLorean:::calc.roughness(a_perm)
  	  }
  	}      
  	PermMeanMatrix <- rowMeans(r_perm_ofprotein, na.rm = FALSE, dims = 1)
  	AlgRoughMatrix<-data.frame(observedRoughness)
  	permRoughList <- PermMeanMatrix[1:100]
  
  	TEST <- t.test(x = observedRoughness, y = permRoughList, alternative = "less")
  	
  	roughMatrix[i,] <- c(mean(observedRoughness), TEST$p.value)
  }
  return(roughMatrix)
}


