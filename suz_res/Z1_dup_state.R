library(dplyr)
library(ggplot2)
library(UpSetR)
setwd("~/cloud/para/suz_res/")
# cyp28d1 :FBgn0031689
# Ugt86Dh :	FBgn0040252

options()

suz_cds_file <- paste0("./raw_suz_13strains/","ref","_cds.suz")
st_cds <- read.table(file = suz_cds_file,header=F, fill=T, stringsAsFactors = F, na.strings = "na")
names(st_cds) <- c('fbna', 'chrom', 'posit', paste0("cds_",1:5), paste0("hits_",1:5) )

suz_pg_file <- paste0("./raw_suz_13strains/","ref","_pgene.suz")
st_pg <- read.table(file=suz_pg_file, header=F, fill=T, stringsAsFactors = F, na.strings = "na")
names(st_pg)  <- c('fbna',"phits")

res <- merge.data.frame(st_cds, st_pg, by="fbna", all=T) %>% rowwise() %>%mutate(mhits=min(c(hits_1,hits_2,hits_3,hits_4,hits_5,phits), na.rm=TRUE) )

names(res) <- c('fbna','chrom','posit',paste0("ref", "_cds_",1:5), paste0("ref", "_hits_",1:5),paste0("ref", "_phits"),paste0("ref","_mhits"))

chromlen <- data.frame(chrom=c("NT_033778.4"  ,  "NT_037436.4"   , "NT_033779.5"  ,  "NT_033777.3"   , "NC_004354.4"  ,  "NC_024511.2" ,   "NC_004353.4"  , "NW_001846712.1" , "NW_001845284.1", "NC_024512.1"), chrom2=c("2R","3L", "2L", "3R", "X","mito", "4", "Y_p", "X_p","Y"), len=c(25286936, 28110227 , 23513712, 32079331, 23542271, 19524 , 1348131, 1144, 3347, 3667352))

ref_mp <- merge.data.frame(chromlen, res, by="chrom")
###################

strains <- c("a1", "a2","a3","a4","a5","a6","a7","ab8","b1","b2","b3","b4","b6")
#strains <- c("a1","a2")
#strain <- "a1"

turnNato0 <- function(x) {
    if(is.na(x) | is.infinite(x)) {return(0)}
    else {return(x)}
}

#strain="a1"
ref_strain_mp <- ref_mp
for (strain in strains) {
    print(strain)
    
    suz_cds_file <- paste0("./raw_suz_13strains/",strain,"_cds.suz")
    st_cds <- read.table(file = suz_cds_file,header=F, fill=T, stringsAsFactors = F, na.strings = "na")
    names(st_cds) <- c('fbna', 'chrom', 'posit', paste0("cds_",1:5), paste0("hits_",1:5) )

    suz_pg_file <- paste0("./raw_suz_13strains/",strain,"_pgene.suz")
    st_pg <- read.table(file=suz_pg_file, header=F, fill=T, stringsAsFactors = F, na.strings = "na")
    names(st_pg)  <- c('fbna', "phits")
    
    tmp <- merge.data.frame(st_cds,st_pg,by='fbna',all=T) %>% rowwise() %>% mutate(mhits_cds=turnNato0(min(c(hits_1,hits_2,hits_3,hits_4,hits_5), na.rm=TRUE)) ) %>% mutate(mhits = turnNato0( min(c(mhits_cds , phits), na.rm=FALSE))) %>% select(-c(chrom,posit))
    
    names(tmp) <- c('fbna', paste0(strain,"_cds_",1:5), paste0(strain,"_hits_",1:5), paste0(strain, "_phits"), paste0(strain, "_mhits_cds"), paste0(strain,"_mhits"))
    
    ref_strain_mp <- merge.data.frame(ref_strain_mp,tmp, by='fbna', all=T)
    
}


############  UpSet Plot ################# 

#need <- grep("^[a-b]+\\d+$", names(ref_strain_mp), value = T, perl=T)

ref_strain_tf <- ref_strain_mp %>% 
    mutate(dtf_a1=a1_mhits > 1 ,dtf_a2 = a2_mhits > 1 ,dtf_a3 =a3_mhits > 1 ,dtf_a4 =a4_mhits > 1 ,dtf_a5=a5_mhits > 1 ,dtf_a6 =a6_mhits > 1 ,dtf_a7 =a7_mhits > 1 ,dtf_ab8= ab8_mhits > 1 ,dtf_b1 = b1_mhits > 1 ,dtf_b2 = b2_mhits > 1 , dtf_b3 = b3_mhits > 1 , dtf_b4 = b4_mhits > 1 , dtf_b6 = b6_mhits > 1 ) %>% select(fbna, chrom2,posit, starts_with("dtf_"))

dim(ref_strain_mp)


colSums(ref_strain_tf[,4:16], na.rm=T)


tf_slic <- ref_strain_tf[,c(4:16)]  * 1 
tf_slic <- apply(tf_slic, FUN = turnNato0, MARGIN = c(1,2))
tf_slic <- as.data.frame(tf_slic)
names(tf_slic) <- substr(names(tf_slic), start=5, stop=20)
upset(tf_slic, nsets = 13, order.by = "freq")

# there are in total 531 genes. 
tf_slic_dup <- tf_slic %>% filter(a1 ==T|a2==T|a3==T|a4==T|a5==T|a6==T|a7==T|ab8==T|b1==T|b2==T|b3==T|b4==T|b6==T)

dim(tf_slic_dup)
#[1] 531  13

# extract a3 a4 related genes. 

a3a4_slic <- tf_slic %>% filter(a3 == 1 | a4 == 1)

dim(a3a4_slic)
# [1] 325  13
upset(a3a4_slic,nset=13,order.by="freq")

## check the position of these 325 genes in iso1: 

## summary 
# ref_strain_mp  : 13928  186 contains ref starin all info cds, hits, phit, mhits 
# ref_strain_tf : 13928 16    true of false of dtf: duplication?
# tf_slic : 13928 13
# tf_slic_dup : 531 13
# a3a4_slic: 325 13

