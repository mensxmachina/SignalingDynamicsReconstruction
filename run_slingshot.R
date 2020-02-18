run_slingshot <- function(D, activator, transf, drMethod, randomStart = F, write2disk = T, wd = NULL){

# 1st step => dim reduction
dr_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_",drMethod,".txt") 

# if a file with the requested dimensionality reduction analysis exists, do not repeat 
if ((file.exists(paste0("results/",dr_file_name)) == FALSE) | (randomStart == T)){
  if (drMethod == "tSNE"){ 
    set.seed(24)
    tsne_out <- Rtsne.multicore::Rtsne.multicore(D$expr, dims = 2)
    data_drMethod <- tsne_out$Y
  }
  if (drMethod == "PCA"){
    pca <- prcomp(D$expr, scale. = FALSE)
    data_drMethod <- pca$x[,1:2]
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

# 2nd step => trajectory inference
ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_slingshot_",drMethod,"_output.txt")

# if a file with the inferred trajectory is already created, do not repeat
if ((file.exists(paste0("results/",ti_file_name)) == FALSE)  | (randomStart == T)){
  
  #Clustering
  cl1 <- kmeans(data_drMethod, centers = 5)$cluster
  #Mapping
  if (randomStart){
    lin1 <- slingshot::getLineages(data_drMethod, cl1, start.clus = as.character(sample(1:5,1)))
  } else {
    lin1 <- slingshot::getLineages(data_drMethod, cl1)
  }
  
  #Construct smooth curves 
  crv1 <- slingshot::getCurves(lin1)
  # this is your pseudotime
  pseudot <- slingshot::slingPseudotime(crv1) 

  # 3rd step => Save a table with Normalized_values [1:limit], Experimental time, Pseudotime
  if (length(crv1@curves) > 1) {
    print(paste0("Slingshot with ", drMethod, " found more than one lineages !"))
  }
  
  if (write2disk){
    algs_matrix = cbind(D$expr, D$timepoint, pseudot, crv1@curves$curve1$s)
    write.table(algs_matrix, file = paste0("results/",ti_file_name))
  } else {
    to_return <- list()
    to_return[["dimRed"]] <- data_drMethod
    to_return[["pseudotime"]] <- pseudot
    to_return[["trajectory"]] <- crv1@curves$curve1$s
    to_return[["lineages"]] <- length(crv1@curves)
    
    return(to_return)
  }
  
  print(paste0("Finished with Slingshot with ", drMethod, " and created ", ti_file_name ,"  :) " ))
}

}


# how to plot in the case of bifurcating trajectory
# cl = factor(D$timepoint)
# levels(cl) = 1:D$ntimepoints
# 
# plot(data_drMethod, col = RColorBrewer::brewer.pal(D$ntimepoints,"Set1")[cl])
# slingshot::lines(crv1,show.constraints = TRUE)
# legend("bottomright", legend=unique(D$timepoint),
#        col=RColorBrewer::brewer.pal(D$ntimepoints,"Set1"), pch=19, cex=0.8)
# plot(algs_matrix$curve1,D$expr[,2], col = RColorBrewer::brewer.pal(D$ntimepoints,"Set1")[cl])

