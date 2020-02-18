run_DeLorean <- function(D, activator, transf, drMethod=NULL, randomStart = F, write2disk = T, wd = NULL){
  
# 1st step => dim reduction

ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_DeLorean_output.txt")

inner_datanorm <- D$expr

guo.expr = t(inner_datanorm)
guo.gene.meta = as.data.frame(colnames(inner_datanorm))
names(guo.gene.meta) <- "gene"
cell = as.factor(1:dim(inner_datanorm)[1])
capture = as.factor(D$timepoint)
obstime = as.factor(D$timepoint)
# yes, you need to do this 2 times.
guo.cell.meta1 <- cbind(cell, capture, obstime)
guo.cell.meta <- as.data.frame(guo.cell.meta1)
guo.cell.meta$cell <- as.factor(guo.cell.meta$cell)
guo.cell.meta$capture <- as.factor(guo.cell.meta$capture)
guo.cell.meta$obstime <- as.numeric(guo.cell.meta$obstime)
colnames(guo.expr) <- guo.cell.meta$cell

dl <- DeLorean::de.lorean(guo.expr, guo.gene.meta, guo.cell.meta)

dl <- DeLorean::estimate.hyper(dl, model.name='lowrank')
# by placing sigma.tau=0.5, length.scale=1.5 => error... now (7/5) try without any prior.
# sigma.tau = 0.5
# length.scale = NULL (Length scale for stationary GP covariance function. Defaults to the range of the observed capture times.)
# exact': The model without a low rank approximation that does not estimate the cell sizes.
# adjust.cell.sizes:  Adjust by the cell sizes for better estimates of the hyperparameters

dl <- DeLorean::fit.dl(dl, method = "vb") # through method you can choose how many cores you want. https://www.rdocumentation.org/packages/DeLorean/versions/1.2.5/topics/fit.dl
dl <- DeLorean::examine.convergence(dl)
pseudot <- dl$cell.map$S.hat

# 3rd step => Save a table with Normalized_values [1:limit], Experimental time, Pseudotime
if (write2disk){
  algs_matrix = cbind(inner_datanorm, D$timepoint, pseudot)
  write.table(algs_matrix, file = paste0("results/",ti_file_name))
  print(paste0("Finished with DeLorean and created", ti_file_name ,"  :) \n " ))
} else {
  to_return <- list( dimRed = inner_datanorm, pseudotime = pseudot, trajectory = NULL, lineages = 1)
  print(paste0("Finished with DeLorean and created", ti_file_name ,"  :) \n " ))
  return(to_return)
}


}
