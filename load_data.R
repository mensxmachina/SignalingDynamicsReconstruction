load.Data <- function(files.paths, index, measurement.time, transf = 'asinh'){
  
  data.raw = read.flowSet(files = dataFiles)
  
  data.info <- pData(parameters(data.raw[[1]]))
  
  # make a flowSet with only the abundances of the markers of choice. (optional: arcsinh transform)
  data <- fsApply(data.raw, function(x, cofactor=5){
    colnames(x) <- data.info$desc
    expr <- flowCore::exprs(x)
    if (transf == 'asinh') {
      expr <- asinh(expr[, index]/ cofactor)  
    } else if (transf == 'log') {
      ## Automatically estimate the logicle transformation based on the data
      lgcl <- estimateLogicle(samp, channels = index)
      ## transform  parameters using the estimated logicle transformation
      expr <- transform(expr[, index], lgcl) 
    } else {
      expr <- expr[, index]
    }
    exprs(x) <- expr
    x
  })
  data
  
  # how many cells are in each dataset
  N <- fsApply(data, function(x){ x <- dim(flowCore::exprs(x))[1]})
  
  # make a matrix with the abundances and a vector with the time
  expr <- as.data.frame(fsApply(data, flowCore::exprs))
  timepoint <- unlist(sapply(c(1:length(N)), function(i){ rep(measurement.time[i],N[i])} ))
  
  # find the mean at each timepoint
  exp_means <- fsApply(data,each_col,mean)
  
  D <- list()
  D$N <- N
  D$expr <- expr
  D$exp_means <- exp_means
  D$timepoint <- timepoint
  D$ntimepoints <- length(N)
  D$info <- data.info
  
  return (D)
}