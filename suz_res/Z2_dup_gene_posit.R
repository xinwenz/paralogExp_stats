library(dplyr)
library(ggplot2)
library(tidyverse)
# run Z1_dup_state.R first to get ??? 
setwd("~/cloud/para/suz_res/")
load("./Z1_dup_state")

## all the 13914 genes annotated in ISO1
ref_strain_tf <- ref_strain_tf %>% rowwise() %>% mutate(sm = turnNato0(sum(c(dtf_a1,dtf_a2,dtf_a3,dtf_a4,dtf_a5,dtf_a6,dtf_a7,dtf_ab8,dtf_b1,dtf_b2,dtf_b3,dtf_b4,dtf_b6), na.rm = TRUE)))
ref_strain_tf$sm <- as.factor(ref_strain_tf$sm)
cols <- c("0" = "black", "1" = "red", "2" = "red", "3" = "red",
          "4" = "yellow", "5" = "yellow", "6" = "yellow", "7" = "orange",
          "8" = "orange", "9" = "orange", "10" = "cyan", "11" = "cyan",
          "12" = "cyan", "13" = "cyan", "NA" = "grey50")

ggplot(data = chromlen) + 
    geom_bar(aes(x=chrom2, y=len),stat = "identity", fill=NA, color="black") + 
    geom_point(data=ref_strain_tf, aes(x=chrom2, y= posit, color=sm),alpha=0.5 ,size=2, shape=20) +
    scale_color_manual(values=cols) + 
    theme(legend.title = element_blank())

##  check 531 genes in 2L 2R 3L 3R 4 X. 
ref_strain_tf_part <- ref_strain_tf %>% filter(dtf_a1 ==T|dtf_a2==T|dtf_a3==T|dtf_a4==T|dtf_a5==T|dtf_a6==T|dtf_a7==T|dtf_ab8==T|dtf_b1==T|dtf_b2==T|dtf_b3==T|dtf_b4==T|dtf_b6==T)

cols2 <- c( "1" = "red", "2" = "red", "3" = "red",
          "4" = "black", "5" = "black", "6" = "black", "7" = "blue",
          "8" = "blue", "9" = "blue", "10" = "cyan", "11" = "cyan",
          "12" = "cyan", "13" = "cyan", "NA" = "grey50")
ggplot(data = chromlen) + 
    geom_bar(aes(x=chrom2, y=len),stat = "identity", fill=NA, color="black") + 
    geom_point(data=ref_strain_tf_part, aes(x=chrom2, y= posit, color=sm), size=2, shape=20) + 
    scale_color_manual(values=cols2) + 
    theme(legend.title = element_blank())

ggplot() + geom_histogram(aes(x=as.numeric(ref_strain_tf_part$sm)))

##  check 325 gene postions in iso1 
ref_strain_tf_a3a4 <- ref_strain_tf %>% filter(dtf_a3 ==T | dtf_a4 == T)

ggplot(data = chromlen) + 
    geom_bar(aes(x=chrom2, y=len),stat = "identity", fill=NA, color="black") + 
    geom_point(data=ref_strain_tf_a3a4, aes(x=chrom2, y= posit, color=sm), size=2, shape=20) + 
    scale_color_manual(values=cols2) + 
    theme(legend.title = element_blank())

ggplot() + geom_histogram(aes(x=as.numeric(ref_strain_tf_a3a4$sm)))

