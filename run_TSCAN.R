run_TSCAN <- function(D, activator, transf, drMethod = NULL, randomStart = F, write2disk = T, wd = NULL){

  # 1st step => dim reduction
  dr_file_name1 = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_PCA_MST.txt")
  dr_file_name2 = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_PCA_MST.Rdata")

  # if a file with this dim red was already created, dont recreate
  if (file.exists(paste0("results/",dr_file_name2)) == FALSE | (randomStart == T)){
    mclust <- TSCAN::exprmclust(as.data.frame(t(D$expr)))
    
    # TSCAN::plotmclust(mclust,show_cell_names = F)
    if (write2disk){
      save(mclust, file = file.path(wd, "results",dr_file_name2))
      write.table(mclust$pcareduceres, file = file.path(wd, "results",dr_file_name1))
    }
  } else {load(file = paste0("results/",dr_file_name2))}
  
  # 2nd step => trajectory 
  ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_TSCAN_output.txt")
  
  # if a file with this analysis is already created, dont reapply 
  if ((file.exists(paste0("results/",ti_file_name)) == FALSE) | (randomStart == T)){
    
    if (randomStart){
      ordering <- TSCAN::TSCANorder(mclust, MSTorder = sample(unique(mclust$clusterid)), orderonly = F, flip = F, listbranch = F)
    } else {
      ordering <- TSCAN::TSCANorder(mclust, MSTorder = NULL, orderonly = F, flip = F, listbranch = F)
    }
     
    TSCAN::plotmclust(mclust,show_cell_names = F, MSTorder = sample(unique(mclust$clusterid)))
    
    ordering$sample_name <- as.integer(sub("V","",ordering$sample_name))
    colnames(ordering)[1] <- "cell_ID"
    i = sort(ordering$cell_ID, index.return = T)$ix
    
    # 3rd step => Save a table with the results
    if (nrow(ordering) == nrow(D$expr)){
      pseudot <- ordering$Pseudotime[i]
      algs_matrix = cbind(D$expr, D$timepoint, pseudot)
      lineages <- 1
    } else { 
      algs_matrix = cbind(D$expr, D$timepoint, c(1:nrow(D$expr))) 
      lineages <- 2
    }
      
    if (write2disk){
      write.table(algs_matrix, file = paste0("results/",ti_file_name))
    } else {
      to_return <- list()
      to_return[["dimRed"]] <- mclust$pcareduceres
      to_return[["pseudotime"]] <- pseudot
      to_return[["trajectory"]] <- mclust$MSTtree
      to_return[["lineages"]] <- lineages
      
      return(to_return)
    }
    print(paste0("Finished with TSCAN and created ", ti_file_name ,"  :) \n " ))
  }

}
