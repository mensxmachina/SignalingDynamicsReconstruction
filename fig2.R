fig2 <- function(activator, transf, resultFiles, shift.start = T, circular = F){
  
  reswd = sub("/int.*","",resultFiles[1])
  wd = sub("/results","",reswd)
  figwd = paste0(reswd,"/figures")
  dir.create(file.path(figwd), showWarnings = FALSE)
  setwd(figwd)
  
  for(i in 1:length(resultFiles)){
    
    # identify the TI method used (cr: curve reconstruction, dr: dimensionality reduction)
    alg = sub("_output.*","",sub(paste0(".*",transf,"_"),"",resultFiles[i]))
    cr.method = unlist(strsplit(alg,"_"))[1]
    dr.method = unlist(strsplit(alg,"_"))[2]
    
    # create folder to save the figures
    fp = file.path(figwd, alg)
    dir.create(fp, showWarnings = FALSE)
    
    # load the data
    input_matrix <- read.table (file = resultFiles[i])
    
    # load the dimensionality reduction method results
    dr.File <- list.files(file.path(reswd), pattern=paste0(".*",dr.method,".txt"), full.names=TRUE)
    if (length(dr.File) == 0) {
      dr.File <- list.files(file.path(reswd), pattern=paste0(".*",dr.method,".Rdata"), full.names=TRUE)
      load(file = dr.File)
    } else {
      dr_matrix <- read.table (file = dr.File)
      colnames(dr_matrix) = c("V1","V2")
    }
    
    
    # calculate plot parameters
    npar = which(colnames(input_matrix) == "D.timepoint") - 1
    ntimepoints = length(unique(input_matrix$D.timepoint))
    
    # extract data to plot
    dat <- input_matrix[,1:npar]
    Timepoints = input_matrix[, npar + 1]
    Pseudotime <- input_matrix[,npar + 2]
    trajectory = input_matrix[, c(npar + 3,npar + 4)]
    
    # Find axis limits => you need common limits of asinh or log for y axis
    y_min <- min(dat)
    y_max <- max(dat)

    # sort the plot values in ascending pseudotime order
    ix = order(Pseudotime)
    dat = dat[ix,]
    Timepoints = Timepoints[ix]
    Pseudotime = Pseudotime[ix]
    if (sum(is.na(Pseudotime)) > 0 ){
      print(paste0("The algorithm ", alg," returned a branching trajectory. Cannot visualize !"))
      next
    }
    
    trajectory = trajectory[ix,]
    colnames(trajectory) = c("V1","V2")

    # rescale pseudotime to [0,1] for the figures between algorithms to be comparable
    Pseudotime_zscore <- (Pseudotime-min(Pseudotime,na.rm=T))/(max(Pseudotime,na.rm=T)-min(Pseudotime,na.rm=T))
    
    # reposition the pseudotime values so that cells in experimental time 0 are close to the origin
    if (shift.start){
      Pseudotime_zscore = shift_start(exp.time = Timepoints, pseudotime = Pseudotime_zscore, ps.ordered = T, circular = circular, randomize_t0 = T)
    }

    # add the spline fit
    spline.df = calc_spline(dat,Pseudotime_zscore)
    dfm <- reshape2::melt(spline.df, id.vars='time')
    
    cols <- c(RColorBrewer::brewer.pal(12,"Paired"), RColorBrewer::brewer.pal(9,"Set1")) 
    
    param.plot <- list()
    # repeat for each cell parameter
    for (j in 1:npar){
      plotvar <- dat[,j]
      to_plot = data.frame(pst=Pseudotime_zscore,y=plotvar,t = as.factor(Timepoints))

      # make the plot
      param.plot[[j]] <- ggplot(to_plot,aes(pst, y)) + geom_point(aes(color = t), alpha=.5, size = 0.8) +
                            scale_color_manual(values = cols, name = "time") + ylim(y_min,y_max) + xlim(0,1) + 
                            geom_line(data = subset(dfm, variable == names(dat)[j]), aes(x = time, y = value), size = 1) +
                            # ggtitle(paste0(alg,"_", colnames(dat)[j],"_lineage_",l)) + 
                            ggtitle(colnames(dat)[j]) +
                            # xlab('scaled pseudotime') + ylab(paste0('abundance (',transf,")")) + 
                            theme(legend.position="none", legend.title=element_text(size=8), 
                                  plot.title = element_text(hjust = 0.5, size = 10),
                                  axis.text.x=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                  axis.ticks.x=element_blank(), #axis.line.y = element_line(colour = "black"),
                                  panel.background = element_blank()) +
                            guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2),nrow = 1,title="experimental time"))
      
      # check the density of the reordered experimental time
      # ggplot(to_plot,aes(pst, y)) + geom_point(aes(color = t), alpha=.5, size = 0.8) +
      #   scale_color_manual(values = cols, name = "time") + ylim(y_min,y_max) + xlim(0,1) +
      #   geom_density_2d(aes(color = t), alpha = 0.5) + facet_wrap(~t)
      
      
      # save each plot independently
      # ggsave(paste0(fp,"/int_", activator,"_",ntimepoints,"tp_",transf,"_",alg,"_", colnames(dat)[j], "_lineage_", l ,".png"))
    }
    
    y_label <- ggdraw() + draw_label("protein abundance (asinh)", size = 10, angle = 90, hjust = 0.7)
    x_label <- ggdraw() + draw_label("standardized pseudotime", size = 10, hjust = 0.4) #+ theme(plot.margin=margin(t=0.5,unit="cm"))
    
    # draw the parameters
    p0 = plot_grid(plotlist = param.plot, nrow = 4)
    
    # add the legend
    leg <- get_legend(param.plot[[1]]+ theme(legend.position = 'bottom'))
    # p0 = p0 + draw_grob(leg, 0.6,-0.1, 1/4, 0.5)
    
    # create the final subplot figure
    p1 = plot_grid(y_label, p0, rel_widths = c(0.1,1))
    p2 = plot_grid(p1, x_label, rel_heights = c(1,0.05), nrow = 2)
    
    # create the reduced manifold figure
    dr.df = cbind(dr_matrix[ix,], Timepoints)
    dr.plot <- ggplot(dr.df, aes(x=V1,y=V2)) + geom_point(aes(color = as.factor(Timepoints)), alpha=.5, size = 0.8) +
                  scale_color_manual(values = cols, name = "time") +
                  xlab("reduced dimension 1") + ylab("reduced dimension 2") +
                  theme(legend.position="none", panel.background = element_blank(),
                        axis.text=element_blank(), axis.ticks = element_blank(),
                        axis.title = element_text(size = 10)) + 
                  geom_point(data = trajectory, aes(x = V1, y = V2),size = 0.5)
    
    # combine in one figure the manifold and the reconstructed dynamics per parameter
    p3 = plot_grid(dr.plot,p2, labels = c('A','B'), nrow = 1, axis = "t")
    p3 = plot_grid(p3,leg, ncol = 1, rel_heights = c(1,0.2))
    
    ggsave(paste0(fp,"/fig2.tiff"),plot = p3,device = "tiff",dpi = 300,width = 30,height = 18, units ="cm")
  
  }
  setwd(wd)
}
