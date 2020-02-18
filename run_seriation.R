run_seriation <- function(D, activator, transf, drMethod = NULL, randomStart = F, write2disk = T){
  
  # 2nd step => trajectory 
  ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_seriation_output.txt")
  
  # if a file with this analysis is already created, dont reapply 
  if ((file.exists(paste0("results/",ti_file_name)) == FALSE) | (randomStart == T)){
    
    test.ser = seriation::seriate(as.matrix(D$expr + 10),method="PCA_angle")
    ord.x = unlist(test.ser[[1]])
    ord.y = unlist(test.ser[[2]])
    
    pseudot = seriation::get_rank(ord.x)
    
    # 3rd step => Save a table with the results
    if (write2disk){
      algs_matrix = cbind(D$expr, D$timepoint, pseudot)
      write.table(algs_matrix, file = paste0("results/",ti_file_name))
    } else {
      to_return <- list()
      to_return[["dimRed"]] <- NULL
      to_return[["pseudotime"]] <- pseudot
      to_return[["trajectory"]] <- NULL
      to_return[["lineages"]] <- 1
      
      return(to_return)
    }
    print(paste0("Finished with Seriation with PCA_angle and created ", ti_file_name ,"  :) \n " ))
  }

}
