length(unique(pData(my_cds)$State))

plot_cell_trajectory(my_cds, color_by = "State")


cell_metadata = data.frame(cell = c(1:sum(D$N)), time.point = D$timepoint)
row.names(cell_metadata) = paste0("cell",c(1:sum(D$N)))

protein_metadata = data.frame(gene_short_name = names(D$expr))
row.names(protein_metadata) = names(D$expr)

my_cds2 <- monocle3::new_cell_data_set(as(t(D$expr), "sparseMatrix"),#as.matrix(t(D$expr)),
                              cell_metadata = cell_metadata,
                               gene_metadata = protein_metadata)

my_cds2 <- monocle3::preprocess_cds(my_cds2, norm_method = 'none', num_dim = 3)
my_cds2 <- monocle3::reduce_dimension(my_cds2,reduction_method = "UMAP")
my_cds2 <- monocle3::cluster_cells(my_cds2,reduction_method = "UMAP")
my_cds2 <- monocle3::learn_graph(my_cds2)
my_cds2 <- monocle3::order_cells(my_cds2, root_cells = paste0("cell",sample(1:sum(D$N),1)))

X = monocle3::pseudotime(my_cds2)
head(X)

X=igraph::V(monocle3::principal_graph(my_cds2))
X=principal_graph(my_cds2)[["UMAP"]]

plot_cells(my_cds2, color_cells_by = "pseudotime")

monocle3::plot_cells(my_cds2)

X=as(t(D$expr), "sparseMatrix")


