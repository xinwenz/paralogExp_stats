library(dplyr)
library(ggplot2)
library(tidyverse)
# run Z1_dup_state.R first to get ??? 
setwd("~/cloud/para/suz_res/")
load("./Z1_dup_state")
load("./Z2_dup_gene_posit.R")

#### number of alleles and its relationship with prevelance of duplication. 
dim(ref_strain_mp)
ref_strain_sm_numal <- ref_strain_mp %>% 
    mutate(dtf_a1=a1_mhits > 1 ,dtf_a2 = a2_mhits > 1 ,dtf_a3 =a3_mhits > 1 ,dtf_a4 =a4_mhits > 1 ,dtf_a5=a5_mhits > 1 ,dtf_a6 =a6_mhits > 1 ,dtf_a7 =a7_mhits > 1 ,dtf_ab8= ab8_mhits > 1 ,dtf_b1 = b1_mhits > 1 ,dtf_b2 = b2_mhits > 1 , dtf_b3 = b3_mhits > 1 , dtf_b4 = b4_mhits > 1 , dtf_b6 = b6_mhits > 1 ) %>% rowwise() %>%
    mutate(sm = turnNato0(sum(c(dtf_a1,dtf_a2,dtf_a3,dtf_a4,dtf_a5,dtf_a6,dtf_a7,dtf_ab8,dtf_b1,dtf_b2,dtf_b3,dtf_b4,dtf_b6), na.rm = TRUE))) %>% rowwise() %>%
    mutate(numOfalle = length(table(c(a1_mhits,a2_mhits,a3_mhits,a4_mhits,a5_mhits,a6_mhits,a7_mhits,ab8_mhits,b1_mhits,b2_mhits,b3_mhits,b4_mhits,b6_mhits)))) %>% 
    #filter(dtf_a1 ==T|dtf_a2==T|dtf_a3==T|dtf_a4==T|dtf_a5==T|dtf_a6==T|dtf_a7==T|dtf_ab8==T|dtf_b1==T|dtf_b2==T|dtf_b3==T|dtf_b4==T|dtf_b6==T) %>% 
    mutate(a3a4 = (dtf_a3 == T | dtf_a4 ==T) ) %>%
    select(fbna, chrom2,posit, sm, numOfalle, a3a4)

dim(ref_strain_sm_numal) #13928 6
names(ref_strain_sm_numal)

ggplot(data=ref_strain_sm_numal) + geom_jitter(aes(x=sm, y=numOfalle, color=a3a4), width=0.2,height=0.2, size=0.1) +
    scale_y_continuous(breaks=seq(0,13,1)) + 
    scale_x_continuous(breaks=seq(0,13,1))


################
ref_a3a4_mp <- ref_strain_mp %>% filter(a3_mhits > 1 | a4_mhits > 1) %>% rowwise() %>% 
    #mutate(a3_mcd = min(a3_hits_1,a3_hits_2,a3_hits_3,a3_hits_4,a3_hits_5, na.rm = T),  a4_mcd = min(a4_hits_1,a4_hits_2,a4_hits_3,a4_hits_4,a4_hits_5, na.rm = T)) %>% 
    select(fbna, chrom, chrom2,len, posit, starts_with("ref_"), starts_with("a3_"), starts_with("a4_"))

ggplot() + geom_histogram(aes(x=ref_a3a4_mp$a3_mhits))
ggplot() + geom_histogram(aes(x=ref_a3a4_mp$a4_mhits))

names(ref_a3a4_mp)
mhits_melt <- melt(ref_a3a4_mp, id.vars = "fbna", measure.vars = c("a3_mhits","a4_mhits"))
names(mhits_melt)

ggplot(data=mhits_melt) + geom_histogram(aes(x=value, fill=variable), alpha=0.5, binwidth = 1) + scale_x_continuous(breaks = seq(0,45,2))


###########
ref_a3a4_freq_cp <- merge.data.frame(ref_a3a4_mp, ref_strain_tf_a3a4, by=c("fbna","chrom2","posit"))
ggplot(data=ref_a3a4_freq_cp) + geom_jitter(aes(x=a3_mhits, y=a4_mhits, color=sm),shape=20, width=0.18,height=0.18) +
    scale_color_manual(values=cols) + 
    theme(legend.title = element_blank()) +
    coord_cartesian( xlim = c(0,45), ylim=c(0,45))

########### human readable data frame #### 
ref_a3a4_sm <- ref_a3a4_freq_cp %>% select(-starts_with("dtf")) %>% arrange(chrom2,posit)

## filter ref_a3a4_sm into several groups ### double check their positons  ## 
# 0,2
a3a4_differ_a30_a42 <- ref_a3a4_sm %>% filter(a3_mhits==0 & a4_mhits==2)
dim(a3a4_differ_a30_a42) #2 

# 0,>2
a3a4_differ_a30_a49 <- ref_a3a4_sm %>% filter(a3_mhits==0 & a4_mhits>2)
dim(a3a4_differ_a30_a49)  #0

# 1,2
a3a4_differ_a31_a42 <- ref_a3a4_sm %>% filter(a3_mhits==1 & a4_mhits==2)
dim(a3a4_differ_a31_a42) # 30

# 1,>2
a3a4_differ_a31_a49 <- ref_a3a4_sm %>% filter(a3_mhits==1 & a4_mhits>2)
dim(a3a4_differ_a31_a49) # 8 

# 2,0
a3a4_differ_a32_a40 <- ref_a3a4_sm %>% filter(a3_mhits==2 & a4_mhits==0)
dim(a3a4_differ_a32_a40) # 1

# 2,1
a3a4_differ_a32_a41 <- ref_a3a4_sm %>% filter(a3_mhits==2 & a4_mhits==1)
dim(a3a4_differ_a32_a41) #21

# 2,2
a3a4_differ_a32_a42 <- ref_a3a4_sm %>% filter(a3_mhits==2 & a4_mhits==2)
dim(a3a4_differ_a32_a42) # 73 

# 2, >2
a3a4_differ_a32_a49 <- ref_a3a4_sm %>% filter(a3_mhits==2 & a4_mhits > 2)
dim(a3a4_differ_a32_a49) # 11

# >2,0 
a3a4_differ_a39_a40 <- ref_a3a4_sm %>% filter(a3_mhits > 2 & a4_mhits == 0)
dim(a3a4_differ_a39_a40) # 0

# >2,1
a3a4_differ_a39_a41 <- ref_a3a4_sm %>% filter(a3_mhits > 2 & a4_mhits ==1)
dim(a3a4_differ_a39_a41) # 21

# >2,2
a3a4_differ_a39_a42 <- ref_a3a4_sm %>% filter(a3_mhits > 2 & a4_mhits == 2)
dim(a3a4_differ_a39_a42) #11

# >2,>2, differ
a3a4_differ_a39_a49 <- ref_a3a4_sm %>% filter( a3_mhits >2 & a4_mhits > 2 & a3_mhits != a4_mhits )
dim(a3a4_differ_a39_a49)  # 144

# >2,>2, same
a3a4_same_a39_a49 <- ref_a3a4_sm %>% filter(a3_mhits >2 & a4_mhits > 2 & a3_mhits== a4_mhits)
dim(a3a4_same_a39_a49)  # 19

2 + 30 + 8 + 1 + 21 + 73 + 11  + 5 + 11 + 19 + 144
dim(ref_a3a4_sm)

y <- unlist(ref_a3a4_sm$fbna)
write(y, file="~/cloud/para/2cdsPgene0posit/fbgn325.txt")
