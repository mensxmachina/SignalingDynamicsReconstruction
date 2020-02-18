run_monocle <- function(D, activator, transf, drMethod, randomStart = F, write2disk = T, wd = NULL){
library(monocle)
dr_file_name1 = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_MST_",drMethod,".txt")
dr_file_name2 = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_MST_",drMethod,".Rdata")

# if a file with the requested dimensionality reduction analysis exists, do not repeat 
if (file.exists(paste0("results/",dr_file_name2)) == FALSE | (randomStart == T)){
  # set.seed(24)
  exprs_matrix = as.matrix(t(D$expr))
  pd <- new("AnnotatedDataFrame", data = as.data.frame(D$timepoint))
  rownames(pd@data) = colnames(exprs_matrix)
  cds <- monocle::newCellDataSet(exprs_matrix, phenoData = pd, expressionFamily = VGAM::uninormal())
  # 1st step => dim reduction
  if (drMethod == "ICA"){ 
    my_cds <- monocle::reduceDimension(cds, max_components = 2, reduction_method = "ICA", norm_method = "none")
  }
  if (drMethod == "DDRTree"){
    library('DDRTree')
    my_cds <- monocle::reduceDimension(cds, reduction_method = "DDRTree", norm_method = "none")
  }
  if (write2disk){
    save(my_cds, file = file.path(wd, "results",dr_file_name2))
    write.table(t(my_cds@reducedDimS), file = file.path(wd, "results",dr_file_name1))
  }
} else {load(file = paste0("results/",dr_file_name2))}

# 2nd step => trajectory 
ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_monocle_",drMethod,"_output.txt")

# if a file with this analysis is already created, dont reapply 
if (file.exists(paste0("results/",ti_file_name)) == FALSE | (randomStart == T)){
  # find the traj
  my_cds <- monocle::orderCells(my_cds)
  nstates = length(unique(pData(my_cds)$State))
  if (randomStart){
    my_cds <- monocle::orderCells(my_cds, root_state = as.character(sample(1:nstates,1)))
  }
  pseudot <- pData(my_cds)$Pseudotime
  
  if ( nstates > 1){
    print(paste0("Monocle with ", drMethod, " found more than one lineages !"))
  }
  
  # 3rd step => Save a table with the results
  temp=approx(t(my_cds@reducedDimK)[,1],t(my_cds@reducedDimK)[,2],n = nrow(D$expr))
  if (write2disk){
    algs_matrix = cbind(D$expr, D$timepoint, pseudot,temp$x,temp$y)
    write.table(algs_matrix, file = paste0("results/",ti_file_name))
  } else {
    to_return <- list()
    to_return[["dimRed"]] <- t(my_cds@reducedDimS)
    to_return[["pseudotime"]] <- pseudot
    to_return[["trajectory"]] <- cbind(temp$x,temp$y)
    to_return[["lineages"]] <- nstates
    
    return(to_return)
  }
  print(paste0("Finished with Monocle with ", drMethod, " and created ", ti_file_name ,"  :) \n " ))
} 

if (file.exists(ti_file_name) == T) { 
  print(paste0(alg," already applied :) "))
}
  
} # END Monocle with ICA/DDRTree