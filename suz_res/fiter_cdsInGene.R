library(dplyr)
setwd("~/cloud/para/suz_res/")
# cyp28d1 :FBgn0031689
# Ugt86Dh :	FBgn0040252

suz_cds_file <- paste0("./raw_suz_13strains/","ref","_cds.suz")
st_cds <- read.table(file = suz_cds_file,header=F, fill=T, stringsAsFactors = F, na.strings = "na")
names(st_cds) <- c('fbna', 'chrom', 'posit', paste0("cds_",1:5), paste0("hits_",1:5) )
st_cds <- st_cds %>% group_by(fbna) %>% mutate(mhits = mean(c(hits_1,hits_2,hits_3,hits_4,hits_5) , na.rm=T)) 
names(st_cds) <- c('fbna','chrom','posit',paste0("ref", "_cds_",1:5), paste0("ref", "_hits_",1:5),paste0("ref","_mhits"))

suz_pg_file <- paste0("./raw_suz_13strains/","ref","_pgene.suz")
st_pg <- read.table(file=suz_pg_file, header=F, fill=T, stringsAsFactors = F, na.strings = "na")
names(st_pg)  <- c('fbna', paste0("ref","_phits"))

ref <- merge.data.frame(st_cds, st_pg, by="fbna", all=T) 
#ref_single0 <- ref %>% filter(ref_mhits == 1 & ref_phits == 1) #13212 genes
#ref_singleM1 <- ref %>% filter(ref_mhits == 1 & ref_phits > 1 ) # 252 genes, since the pgene threshould is loose, it sometimes get lots hit
#ref_singleM2 <- ref %>% filter(ref_mhits > 1 & ref_phits == 1) # 85 genes, since some cds may gets multiple hits.  
ref_good_sg <-  ref %>% filter( (ref_hits_1 ==1 | ref_hits_2 == 1 | ref_hits_3 == 1 | ref_hits_4 == 1 | ref_hits_5 ==1 ) & ref_phits == 1) %>% mutate(class="good_sg", numL = 1 : 13286)  # 13286 good singleton
ref_fuzzy_sg <-  ref %>% filter( (ref_hits_1 ==1 | ref_hits_2 == 1 | ref_hits_3 == 1 | ref_hits_4 == 1 | ref_hits_5 ==1 ) & ref_phits > 1) %>% mutate(class = "fuzzy_sg", numL= 1: 299) # 299 fuzzy singleton
ref_para <- ref %>% filter( !(fbna %in% ref_good_sg$fbna |  fbna %in% ref_fuzzy_sg$fbna)) # 344  # manual label ##### 
write.csv(ref_para,file="ref_para.csv")
ref_para_manual <- read.csv(file = "ref_para_manforce.csv", header = T, na.strings = "-")

ref_base <- rbind(ref_good_sg, ref_fuzzy_sg,ref_para_manual)
res <- ref_base
strains <- c("a1", "a2","a3","a4","a5","a6","a7","ab8","b1","b2","b3","b4","b6")
#strains <- c("a1","a2")
for (strain in strains) {
    print(strain)

    suz_cds_file <- paste0("./raw_suz_13strains/",strain,"_cds.suz")
    st_cds <- read.table(file = suz_cds_file,header=F, fill=T, stringsAsFactors = F, na.strings = "na")
    names(st_cds) <- c('fbna', 'chrom', 'posit', paste0("cds_",1:5), paste0("hits_",1:5) )
    st_cds <- st_cds %>% group_by(fbna) %>% mutate(mhits = mean(c(hits_1,hits_2,hits_3,hits_4,hits_5) , na.rm=T)) %>% select(-c(chrom, posit))
    names(st_cds) <- c('fbna', paste0(strain, "_cds_",1:5), paste0(strain, "_hits_",1:5),paste0(strain,"_mhits"))
    
    suz_pg_file <- paste0("./raw_suz_13strains/",strain,"_pgene.suz")
    st_pg <- read.table(file=suz_pg_file, header=F, fill=T, stringsAsFactors = F, na.strings = "na")
    names(st_pg)  <- c('fbna', paste0(strain,"_phits"))
    
    tmp <- merge.data.frame(st_cds,st_pg,by='fbna',all=T)
    res <- merge.data.frame(res,tmp, by='fbna', all=T)
}

a3a4_phits_B1 <- res %>% filter(a3_phits > 1 | a4_phits > 1 ) %>% select(fbna,chrom,posit,class,numL,grep(x=names(res), pattern = "*_[mp]hits_*", value=T)) #464 
tmp1 <- a3a4_mhits_B1 %>% group_by( class,numL) %>% summarise(n()) # 255 groups
write.csv(a3a4_mhits_B1,file="a3a4_mhits_B1.csv")

a3a4_mpmatch <- res %>% group_by(fbna) %>% filter(min(a3_hits_1,a3_hits_2,a3_hits_3,a3_hits_4,a3_hits_5,na.rm=T) == a3_phits & min(a4_hits_1,a4_hits_2,a4_hits_3,a4_hits_4,a4_hits_5,na.rm=T) == a4_phits & (a3_phits > 1 | a4_phits > 1)) %>% select(fbna,chrom,posit,class,numL,grep(x=names(res), pattern = "*_[mp]hits_*", value=T)) # 169 


a3a4_mpmatch_uniq <- a3a4_mpmatch %>% group_by(class, numL) %>% summarise(count=n())  #89 
write.csv(a3a4_mpmatch, file="to_google_drive.csv")
