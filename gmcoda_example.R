#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------
#   This is a simple example of "gmcoda: Graphical model for multiple 
# compositional vectors."
#---------------------------------------
#   Please contact Huaying Fang (hyfang@cnu.edu.cn) for any questions about 
# gmcoda.
#-------------------------------------------------------------------------------
require(huge);
source("gmcoda.R");
#---------------------------------------
#-1 Basic example (no edges in the conditional dependence network)
#---------------------------------------
# 1.1 Generate logistic normal variables
set.seed(100);
n = 100;
p = 50;
p1 =25;
x = matrix(rnorm(n * p), nrow = n);
y1 = exp(x[, 1:p1]) / rowSums(exp(x[, 1:p1]));
y2 = exp(x[, (p1+1):p]) / rowSums(exp(x[, (p1+1):p]));
x.frac.list = list(y1, y2);
totCount1 = round(runif(n = n,  min = 1000, max = 2000));
totCount2 = round(runif(n = n,  min = 1000, max = 2000));
x.count.list = list(round(y1 * totCount1), round(y2 * totCount2));
#---------------------------------------
# 1.2 Run gCoda 
##using fraction
res_gmcoda_frac = gmcoda(data = x.frac.list, isCnt = F);
##using counts
res_gmcoda_count = gmcoda(data = x.count.list, isCnt = T);
#---------------------------------------
# 1.3 Get the estimation of the inverse covariance matrix
{
  cat("Gmcoda using fraction data:\n");
  print(round(res_gmcoda_frac$opt_icov, 2));
  cat("Gmcoda using count data:\n");
  print(round(res_gmcoda_count$opt_icov, 2));
}
#-------------------------------------------------------------------------------