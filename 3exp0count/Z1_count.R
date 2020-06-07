library(dplyr)
library(ggplot2)
library(tidyverse)
# run Z1_dup_state.R first to get ??? 
setwd("~/cloud/para/suz_res/")
load("./Z1_dup_state")
load("./Z2_dup_gene_posit.R")
load("./Z3_check_map.R") # get all class of 325 genes.

# a3=1, a4=2
a31_a42 <- a3a4_differ_a31_a42 %>% arrange(chrom2, posit)
