library(dplyr)
setwd("~/cloud/para/suz_res/")

# cyp28d1 :FBgn0031689
# Ugt86Dh :	FBgn0040252

a3_cds <- read.table("suz_a3_cdsInGene.res",header=F,fill=T,stringsAsFactors = F,na.strings = "na")
a4_cds <- read.table("suz_a4_cdsInGene.res",header=F,fill=T,stringsAsFactors = F,na.strings = "na")
ref_cds <- read.table("suz_ref_cdsInGene.res",header=F,fill=T,stringsAsFactors = F,na.strings = "na")

names(a3_cds) <- c('fbna','chrom','posit',paste0(c("cds_"),1:5),paste0(c("hits_"),1:5) )
names(a4_cds) <- c('fbna','chrom','posit',paste0(c("cds_"),1:5),paste0(c("hits_"),1:5) )
names(ref_cds) <- c('fbna','chrom','posit',paste0(c("cds_"),1:5),paste0(c("hits_"),1:5) )

a4_check <- a4_cds %>% group_by(fbna) %>% mutate(mhits = mean(c(hits_1,hits_2,hits_3,hits_4,hits_5) , na.rm=T) ) 
a3_check <- a3_cds %>% group_by(fbna) %>% mutate(mhits = mean(c(hits_1,hits_2,hits_3,hits_4,hits_5) , na.rm=T) ) 
ref_check <- ref_cds %>% group_by(fbna) %>% mutate(mhits = mean(c(hits_1,hits_2,hits_3,hits_4,hits_5) , na.rm=T) ) 

tmp <- merge.data.frame(a4_check,a3_check,by=1,all=T)

a4_a3_ref_check <- merge.data.frame(tmp,ref_check,by=1,all=T)  %>% filter(mhits.y > 1 | mhits.x > 1 | mhits > 1  ) %>% select(-posit.x, -posit.y, -chrom.x , -chrom.y) %>% arrange(chrom,posit)

# 62 genes
a4_a3_1ref_check <- merge.data.frame(tmp,ref_check,by=1,all=T)  %>% filter((mhits.y >= 2 | mhits.x >= 2 ) & mhits < 2  ) %>% select(-posit.x, -posit.y, -chrom.x , -chrom.y) %>% arrange(chrom,posit)


write.csv(a4_a3_ref_check,file="a4_a3_ref_check.csv")

tmp2 <- merge.data.frame(tmp,ref_check,by=c(1,2,3),all=T)
