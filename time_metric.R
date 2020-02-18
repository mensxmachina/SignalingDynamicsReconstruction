time_metric <- function( resultFiles, nruns = 100) {
                         
  reswd = sub("/int.*","",resultFiles[1])
  wd = sub("/results","",reswd)
  metwd = paste0(reswd,"/metrics")
  dir.create(file.path(metwd), showWarnings = FALSE)
  
  time_metric_result <- data.frame()
  for (f in 1:length(resultFiles)){
                           
    # load the data
    input_matrix <- read.table (file = resultFiles[f])
    
    npar = which(colnames(input_matrix) == "D.timepoint") - 1
    exp.time = input_matrix[, npar + 1]
    pseudo.time <- input_matrix[,npar + 2]
    
    # repeat for the forward and reverce trajectory direction
    P = c(0,0)
    for (d in 1:2){
      if (d == 2){
        pseudo.time <- reverse_pseudotime(pseudo.time)
      }
  
      # decode the experiment
      experiment = sub(".*results/","",sub("_output.*","",resultFiles[f]))
      experiment = unlist(strsplit(experiment,"_"))
      cr.method = experiment[5]
      dr.method = experiment[6]
      
      # skip if the pseudo.time includes branches (slingshot)
      if (sum(is.na(pseudo.time)) > 1) {
        time_metric_result <- rbind(time_metric_result, 0)
        rownames(time_metric_result)[f] <- paste0(cr.method,"_",dr.method)
        next
      }
      
      # create all possible combinations of experimental timepoints
      measurement.times = unique(exp.time)
      exp.time.combns = combn(length(measurement.times),2)
      
      # initialization of P( PT[i]<PT[j] | ET[i]<ET[j], for i != j)
      hits = 0
      
      for ( c in 1:dim(exp.time.combns)[2] ){
        # Create all pairs from a sample of n cells where ET[i] < ET[j]
        i = sample(which(exp.time == measurement.times[exp.time.combns[1,c]]),nruns, replace = T)
        j = sample(which(exp.time == measurement.times[exp.time.combns[2,c]]),nruns, replace = T)
        pairs = rbind( i[rep(1:nruns,1,each = nruns)], j[rep(1:nruns,nruns)])
        
        # Check for every pair whether also PT[i] < PT[j]
        for (k in 1:nruns^2){
          if (pseudo.time[pairs[1,k]] <= pseudo.time[pairs[2,k]]) {  hits <- hits + 1  }
        }
      }
  
    (P[d] = hits / (dim(exp.time.combns)[2]*nruns^2))
    }
    
    time_metric_result <- rbind(time_metric_result, max(P))
    rownames(time_metric_result)[f] <- paste0(cr.method,"_",dr.method)
  }
  colnames(time_metric_result) <- "time.metric"
  return(time_metric_result)
}