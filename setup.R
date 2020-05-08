#CPLEX <- FALSE
CPLEX <- TRUE
library(compiler)
library(kernlab)
library(qlcMatrix)
library(MASS)
library(R.matlab)
library(data.table)

source("./other_functions.R")
source("./data_preprocessing.R")
source("./evaluation_functions.R")
if(CPLEX==TRUE){
  library(Rcplex)
  source("./algorithms.R")
}else{
  library(lpSolveAPI)
  source("./algorithms_lpsolveapi.R")
}
#source("./visualization_beta_ver.R")