fig1 <- function(dat,t, fpath){

  # create folder to save the figure
  reswd = paste0(fpath,"/results")
  figwd = paste0(reswd,"/figures")
  dir.create(figwd, showWarnings = FALSE)
  setwd(figwd)
  
  npar = ncol(dat)
  
  y_min = min(dat)
  y_max = max(dat)
  x_min = min(t)
  x_max = max(t)
  
  spline.df = calc_spline(dat,t)
  dfm <- reshape2::melt(spline.df, id.vars='time')
  
  cols <- c(RColorBrewer::brewer.pal(12,"Paired"), RColorBrewer::brewer.pal(9,"Set1")) 
  
  param.plot <- list()
  # repeat for each cell parameter
  for (j in 1:npar){
    plotvar <- dat[,j]
    to_plot = data.frame(pst=t,y=plotvar,t = as.factor(t))
    
    # make the plot
    param.plot[[j]] <- ggplot(to_plot,aes(pst, y)) + geom_point(aes(color = t), alpha=.5, size = 0.8) +
                          scale_color_manual(values = cols, name = "time") + ylim(y_min,y_max) + xlim(x_min,x_max) + 
                          geom_line(data = subset(dfm, variable == names(dat)[j]), aes(x = time, y = value), size = 1) +
                          ggtitle(colnames(dat)[j]) + scale_x_continuous(breaks=c(x_min,x_max)) +
                          theme(legend.position="none", legend.title=element_text(size=8), legend.title.align=0.5,
                                plot.title = element_text(hjust = 0.5, size = 10),
                                axis.title = element_blank(),
                                axis.line = element_line(colour = "black"),
                                panel.background = element_blank()) +
                          guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2),nrow = 2, title="experimental\ntime"))

    # save each plot independently
    # ggsave(paste0(fp,"/int_", activator,"_",ntimepoints,"tp_",transf,"_",alg,"_", colnames(dat)[j], "_lineage_", l ,".png"))
  }
  # draw the parameters
  p0 = plot_grid(plotlist = param.plot, nrow = 4)
  
  # add the legend
  leg <- get_legend(param.plot[[1]]+ theme(legend.position = 'bottom'))
  p0 = p0 + draw_grob(leg, 0.6,-0.1, 1/4, 0.5)
  
  ggsave("fig1.tiff",plot = p0,device = "tiff",dpi = 300,width = 15,height = 18, units ="cm")
  
  setwd(fpath)
}