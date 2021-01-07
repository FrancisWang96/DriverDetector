library(data.table)
library(ggplot2)
library(tidyverse)
library(cowplot)
setwd('/Users/wangzeyuan/Desktop/DriverGenePathway/R')
source('preprocessing.R')
source('preSmallFunctions.R')
source('plots.R')
source('sub_BMR_allGenes.R')
source('BMRsmallFunctions.R')
source('sub_sigGenes20190409.R')
setwd('/Users/wangzeyuan/Desktop/学术/DriverGenePathway相关/SigGenes20180616/sigGenes/tests')
M <- fread('LUSC.mutations.txt')
C <- fread('exome_full192.coverage.txt')
V <- fread('gene.covariates.txt')
dict <- fread('mutation_type_dictionary_file.txt')
chr_files_directory="hg19"
categ_flag=NaN
output_filestem = "preprocessed"
pre_out <- preprocessing(M,C,dict,V,chr_files_directory="hg19",categ_flag=NaN,output_filestem)

M <- pre_out$M
C <- pre_out$C
V <- pre_out$V
bmr = 1.2e-6

bmr_out <- BMR(pre_out$M,pre_out$C,pre_out$V,bmr=1.2e-6)

BMR_out <- bmr_out
p_class = "allTest"
output_filestem = "sigout"
sigThreshold = 0.1
filter = TRUE
