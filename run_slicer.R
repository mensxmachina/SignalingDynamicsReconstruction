run_slicer <- function(D, activator, transf, drMethod, randomStart = F, write2disk = T, wd = NULL){
  
  # 1st step => dim reduction
  dr_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_",drMethod,".txt")

  # if a file is already created, dont reapply
  if (file.exists(paste0("results/",dr_file_name)) == FALSE | (randomStart == T)){
    # choose best k ( because this can be very slow we subsample the data to 1/5 of the cells)
    s = sample(1:sum(D$N),sum(D$N)*0.2)
    k = SLICER::select_k(D$expr[s,], kmin=5, kmax = 50, by = 5) #10 ?
    # Apply LLE
    traj_lle = lle::lle(D$expr, m=2, k)$Y

    if (write2disk){
      write.table(traj_lle, file = file.path(wd, "results",dr_file_name))
    }
  } else {
    traj_lle <- read.table(file = paste0("results/",dr_file_name))
  }
  
  # 2nd step => trajectory 
  ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_SLICER_",drMethod,"_output.txt")
  
  # if a file with this analysis is already created, dont reapply 
  if ((file.exists(ti_file_name) == FALSE) | (randomStart == T)){
    traj_graph = SLICER::conn_knn_graph(traj_lle,5)
    # traj_graph : An igraph object corresponding to the k-NN graph
    ends = SLICER::find_extreme_cells(traj_graph, traj_lle)
    
    to_return <- list()
    # returns more than one ends most times.
    for (i in 1:length(ends)){
      start = ends[i]
      pseudot = SLICER::cell_order(traj_graph, start)
      ti_file_name = paste0("int_", activator,"_",D$ntimepoints,"tp_",transf,"_slicer_",drMethod,"_output_start_",start,".txt")
      # 3rd step => Save a table with Normalized_values [1:limit], Experimental time, Pseudotime
      if (write2disk){
        algs_matrix = cbind(D$expr, D$timepoint, pseudot)
        write.table(algs_matrix, file = paste0("results/",ti_file_name))
      } else {
        to_return[["dimRed"]] <- traj_lle
        to_return[[paste0("pseudotime_start_",start)]] <- pseudot
        to_return[[paste0("trajectory_start_",start)]] <- SLICER::process_distance(traj_graph,start)
        to_return[[paste0("lineages_start_",start)]] <- 2 # BUG: length(unique(SLICER::assign_branches(traj_graph, start, min_branch_len = 100)))
      }
      print(paste0("Finished with SLICER_",start," and created ", ti_file_name ,"  :) \n " ))
    }
    if (!write2disk){
      return(to_return)
    }
  }

} # END SLICER WITH LLE