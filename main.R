# Script for running all algorithms in a single pipeline

library(flowCore)
library(tidyverse)
library(cowplot)
library(foreach)
library(doParallel)

wd <- "" # enter your working directory
setwd(wd)

### ~ ~ ~ Functions ~ ~ ~ ###
source("load_data.R")
source("plot_markers.R")
source("plot_dr.R")
source("shift_start.R")
source("reverse_pseudotime.R")
source("calc_spline.R")
source("cycle_pseudotime.R")
source("Rough.R")
source("Dist_cor.R")
source("Bio_con.R")
source("check_pairs.R")
source("Robustness.R")
source("time_metric.R")
source("fig1.R")
source("fig2.R")

source("run_slingshot.R")
source("run_scorpius.R")
source("run_monocle.R")
source("run_prinCurves.R")
source("run_TSCAN.R")
source("run_seriation.R")
source("run_pseudogp.R")
source("run_slicer.R")
source("run_DeLorean.R")

#** Bodenmiller 2012 Dataset **##

# Choose right files depending to interverntion and number of timepoints
# activators = {01:'pVO4', 02:'IL3',  03:'IL2',  04:'IL12', 05:'GCSF', 06:'GMCSF', 
#               07:'BCR',  08:'IFNg', 09:'IFNa', 10:'LPS',  11:'PMA-IONO'};
activation = '01'
dataFolder = "" # enter your data folder
dataFiles = list.files(dataFolder, pattern=paste0("cd4.*",activation,".fcs"), full.names=TRUE)

t = as.matrix(c(0,1,5,15,30,60,120,240))

param.functional = c(7, 8, 13, 15, 16, 17, 18, 19, 22, 24, 25, 26, 28, 30)
param.lineage = c(3, 4, 9, 11, 12, 14, 21, 29, 31, 33)


#############
#** BEGIN **#
#############
dir.create(file.path(wd, "results"), showWarnings = FALSE)

measurement.time = t
index = param.functional

# Load data
transformation = 'asinh'
proteins = "Phosph"

D = load.Data(dataFiles, index, measurement.time, transformation)

# draw the protein abundance against experimental time
fig1(D$expr,D$timepoint, fpath = wd)

run_slingshot(D, activation, transformation,"tSNE")
run_slingshot(D, activation, transformation,"PCA")
run_slingshot(D, activation, transformation,"diffMaps")
run_scorpius(D, activation, transformation, "distPear")
run_scorpius(D, activation, transformation, "distSpear")
run_scorpius(D, activation, transformation, "distCos") # error in the eigenvalue decomposition
run_scorpius(D, activation, transformation, "distEucl")
run_scorpius(D, activation, transformation, "distManh")
run_monocle(D, activation, transformation, "ICA") # high sample (cell) dimensionality does not allow for ordering the cells
run_monocle(D, activation, transformation, "DDRTree")
run_prinCurves(D, activation, transformation, "tSNE")
run_prinCurves(D, activation, transformation, "diffMaps")
run_TSCAN(D, activation, transformation)

# the next algorithms take long time to run so subsampling is recommended (around 2000 cells should be sufficient)
subsample = F
if (subsample) {
  set.seed(123)
  nss = sum(D$N)*0.2
  s = sample(c(1:sum(D$N)), nss, replace = F)
  
  Dsub = D
  Dsub$expr = D$expr[s,]
  Dsub$N = rep(round(nss/D$ntimepoints),D$ntimepoints)
  Dsub$timepoint = D$timepoint[s]
}
run_DeLorean(Dsub, activation, transformation)
run_slicer(Dsub, activation, transformation, "lle")

resultFiles = list.files(file.path(wd, "results"), pattern=paste0(".*output.*"), full.names=TRUE)

fig2(activation, transformation, resultFiles,shift.start = T, circular = T) # some algorithms (TSCAN) do not conform with this function

# Roughness metric
R = Rough(resultFiles)

# Distance correlation metric
dC = Dist_cor(resultFiles)

# Biological consistency metric
BC = Bio_con(D, "Pathway_Hierarchy.csv", resultFiles, nruns=28)

# Robustness metric
Rob = Robustness(resultFiles, nruns = 7, cell.subset = 0.2)

# Time metric
Tm = time_metric(resultFiles,nruns = 100)
