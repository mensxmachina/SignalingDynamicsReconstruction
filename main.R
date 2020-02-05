# Script for running all algorithms in a single pipeline

library(flowCore)
library(tidyverse)
library(cowplot)

# setwd("~/Pipeline_results")
wd <- "E:/Projects/CAUSALPATH/Papers/Pseudotime/Code/R/analysisCode/"
setwd(wd)

### ~ ~ ~ Functions ~ ~ ~ ###
source("load_data.R")
source("plot_markers.R")
source("plot_dr.R")
source("shift_start.R")
source("reverse_pseudotime.R")
source("calc_spline.R")
source("cycle_pseudotime.R")
source("fig2.R")

source("run_slingshot.R")
source("run_pseudogp.R")
source("run_slicer.R")
source("run_scorpius.R")
source("run_monocle.R")
source("run_prinCurves.R")
source("run_TSCAN.R")
source("run_DeLorean.R")
source("run_seriation.R")



#** Bodenmiller 2012 Dataset **##

# Choose right files depending to interverntion and number of timepoints
# activators = {01:'pVO4', 02:'IL3',  03:'IL2',  04:'IL12', 05:'GCSF', 06:'GMCSF', 
#               07:'BCR',  08:'IFNg', 09:'IFNa', 10:'LPS',  11:'PMA-IONO'};
act = '01'
dataFolder = "E:/Projects/CAUSALPATH/Data/Bodenmiller, Bendall Data/_fcs/BM/8_TimePoints/"
dataFiles = list.files(dataFolder, pattern=paste0("cd4.*",act,".fcs"), full.names=TRUE)

t = as.matrix(c(0,1,5,15,30,60,120,240))

param.func = c(7, 8, 13, 15, 16, 17, 18, 19, 22, 24, 25, 26, 28, 30)
param.line = c(3, 4, 9, 11, 12, 14, 21, 29, 31, 33)


#** Bodenmiller 2017 Dataset **##
act = 'EGF'
dataFolder = "E:/Projects/CAUSALPATH/Data/ProteinOverexpression/_fcs/replicate_1/"
dataFiles = list.files(dataFolder, pattern="cell_Empty", full.names=TRUE)

t = as.matrix(c(60,0,5,15,30))

param.func = c(12, 14:39, 41:47) 
param.line = c(3, 4, 9, 11, 12, 14, 21, 29, 31, 33)


#** Krishnaswamy 2014 Dataset **##
act = 'costim'
dataFolder = "E:/Projects/CAUSALPATH/Data/DREMI-DREVI/_fcs/CD4_naive_series2_CD4_costim/"
dataFiles = list.files(dataFolder, pattern=".fcs", full.names=TRUE)

t = as.matrix(c(0,10,160,1,20,2,0.5,3,40,4,5,6,80,8)) # cd3-cd28
t = as.matrix(c(0,1,20,2,0.5,3,4,6,80,8)) # cd3-cd28-costim

ii = c(1,7,4,6,10,13)
dataFiles <- dataFiles[ii]
t = t[ii]

param.func = c(23,29,28,16) # c(5, 7, 10, 13:23, 26, 28:30)
param.line = c(3, 4, 6, 8, 9, 11, 24, 25, 27)


#############
#** BEGIN **#
#############
dir.create(file.path(wd, "results"), showWarnings = FALSE)

measurement.time = t
index = param.func

# Load data
transf = 'asinh'
proteins = "Phosph"

D = load.Data(dataFiles, index, measurement.time, transf)

subsample = F

if (subsample) {
  set.seed(123)
  nss = 1e4
  
  tD = D
  s = sample(c(1:sum(D$N)), nss, replace = F)
  D$N = rep(round(nss/D$ntimepoints),D$ntimepoints)
  D$expr = D$expr[s,]
  D$timepoint = D$timepoint[s]
}

run_slingshot(D, act, transf,"tSNE")
run_slingshot(D, act, transf,"PCA")
run_slingshot(D, act, transf,"diffMaps")
run_scorpius(D, act, transf, "distPear")
run_scorpius(D, act, transf, "distSpear")
run_scorpius(D, act, transf, "distCos") # error in the eigenvalue decomposition
run_scorpius(D, act, transf, "distEucl")
run_scorpius(D, act, transf, "distManh")
run_monocle(D, act, transf, "ICA") # high sample (cell) dimensionality does not allow for ordering the cells
run_monocle(D, act, transf, "DDRTree")
run_prinCurves(D, act, transf, "tSNE")
run_prinCurves(D, act, transf, "diffMaps")
run_TSCAN(D, act, transf)
run_seriation(D, act, transf)
run_pseudogp(D, act, transf, "tSNE")
run_pseudogp(D, act, transf, "PCA")
run_slicer(D, act, transf, "lle")
run_DeLorean(D, act, transf)

resultFiles = list.files(file.path(wd, "results"), pattern=paste0(".*output.*"), full.names=TRUE)

plot_markers(act, transf, resultFiles,shift.start = T, circular = F) # these are not indicative figures, must show the density over pseudotime


fig1(D$expr,D$timepoint, fpath = wd)

fig2(act, transf, resultFiles,shift.start = T, circular = F) # these are not indicative figures, must show the density over pseudotime



plot_dr(D,paste0(file.path(wd, "results"),'/int_costim_10tp_asinh_tSNE.txt'))


#all_plots <- plotting(intervention = intervention, number_of_timepoints = number_of_timepoints, Normalization = Normalization, Algorithms = plot_algs)












