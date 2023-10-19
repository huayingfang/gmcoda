#!/usr/bin/env Rscript
#-------------------------------------------------------------------------------
aa = c("glasso", "huge", "SpiecEasi");
bb = sapply(aa, require, char = T, warn = F, quiet = F);
source("gmcoda.R");
file1 = "OTU_counts_bacteria_example_gmcoda.csv";
file2 = "OTU_counts_fungi_example_gmcoda.csv";
file_out = "output_gmcoda_example.csv";
dat1 = data.matrix(read.table(file1, head = T, sep = ","));
dat2 = data.matrix(read.table(file2, head = T, sep = ","));
out = gmcoda_wrap(mat1 = dat1, mat2 = dat2, propFilt = 0.8, only.gmcoda = F,);
write.table(xxb, file = file_out, sep = ",", quote = F, row.names = F, col.names = T);
#-------------------------------------------------------------------------------
#-
#-------------------------------------------------------------------------------
#-