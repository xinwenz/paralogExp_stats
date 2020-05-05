library(dplyr)
setwd("~/cloud/para/suz_res/")
# cyp28d1 :FBgn0031689
# Ugt86Dh :	FBgn0040252
#options(digits=2)
options()

suz_cds_file <- paste0("./raw_suz_13strains/","ref","_cds.suz")
st_cds <- read.table(file = suz_cds_file,header=F, fill=T, stringsAsFactors = F, na.strings = "na")
names(st_cds) <- c('fbna', 'chrom', 'posit', paste0("cds_",1:5), paste0("hits_",1:5) )

#st_cds <- st_cds %>% rowwise() %>% mutate(mhits = round(min(c(hits_1,hits_2,hits_3,hits_4,hits_5), na.rm=T),3)) 
#names(st_cds) <- c('fbna','chrom','posit',paste0("ref", "_cds_",1:5), paste0("ref", "_hits_",1:5),paste0("ref","_mhits"))

suz_pg_file <- paste0("./raw_suz_13strains/","ref","_pgene.suz")
st_pg <- read.table(file=suz_pg_file, header=F, fill=T, stringsAsFactors = F, na.strings = "na")
names(st_pg)  <- c('fbna',"phits")

res <- merge.data.frame(st_cds, st_pg, by="fbna", all=T) %>% mutate(mhits=min(c(hits_1,hits_2,hits_3,hits_4,hits_5,phits)) )

names(res) <- c('fbna','chrom','posit',paste0("ref", "_cds_",1:5), paste0("ref", "_hits_",1:5),paste0("ref", "_phits"),paste0("ref","_mhits"))

###################

strains <- c("a1", "a2","a3","a4","a5","a6","a7","ab8","b1","b2","b3","b4","b6")
#strains <- c("a1","a2")
strain <- "a1"
for (strain in strains) {
    print(strain)
    
    suz_cds_file <- paste0("./raw_suz_13strains/",strain,"_cds.suz")
    st_cds <- read.table(file = suz_cds_file,header=F, fill=T, stringsAsFactors = F, na.strings = "na")
    names(st_cds) <- c('fbna', 'chrom', 'posit', paste0("cds_",1:5), paste0("hits_",1:5) )
    #st_cds <- st_cds %>% rowwise() %>% mutate(mhits = min(c(hits_1,hits_2,hits_3,hits_4,hits_5) , na.rm=T)) %>% select(-c(chrom, posit))
    #names(st_cds) <- c('fbna', paste0(strain, "_cds_",1:5), paste0(strain, "_hits_",1:5),paste0(strain,"_mhits"))
    
    suz_pg_file <- paste0("./raw_suz_13strains/",strain,"_pgene.suz")
    st_pg <- read.table(file=suz_pg_file, header=F, fill=T, stringsAsFactors = F, na.strings = "na")
    #names(st_pg)  <- c('fbna', paste0(strain,"_phits"))

    tmp <- merge.data.frame(st_cds,st_pg,by='fbna',all=T) %>% mutate(mhits=min(c(hits_1,hits_2,hits_3,hits_4,hits_5,phits)) ) %>% select(-c(chrom, posit))
    names(tmp) <- c('fbna', paste0(strain, "_cds_",1:5), paste0(strain, "_hits_",1:5), paste0(strain,"_phits") ,paste0(strain,"_mhits"))
    
    res <- merge.data.frame(res,tmp, by='fbna', all=T)
}

################## 
# find any duplication in all the founder strains.... by gene length.... 
observ1 <- res %>% filter(a1_phits > 1 | a2_phits > 1 |a3_phits > 1 |a4_phits > 1 |a5_phits > 1 |a6_phits > 1 |a7_phits > 1 |ab8_phits > 1 | b1_phits > 1 | b2_phits > 1 | b3_phits > 1 | b4_phits > 1 | b6_phits > 1 ) %>% select(fbna,chrom,posit,grep(x=names(res), pattern = "*_[mp]hits_*", value=T)) # 964

## the duplications should be polymorphic, so variance exist : 
observ1_1 <- res %>% rowwise() %>% filter(a1_phits > 1 | a2_phits > 1 |a3_phits > 1 |a4_phits > 1 |a5_phits > 1 |a6_phits > 1 |a7_phits > 1 |ab8_phits > 1 | b1_phits > 1 | b2_phits > 1 | b3_phits > 1 | b4_phits > 1 | b6_phits > 1 ) %>% mutate(ckvar = var(c(a1_phits , a2_phits ,a3_phits ,a4_phits ,a5_phits ,a6_phits ,a7_phits,ab8_phits , b1_phits , b2_phits , b3_phits , b4_phits , b6_phits), na.rm=TRUE)  ) %>% filter(ckvar > 0.01) %>% select(fbna,chrom,posit,ckvar,grep(x=names(res), pattern = "*_[mp]hits_*", value=T))  # 758


################ 
observ2 <- res %>% filter(a1_mhits > 1 | a2_mhits > 1 |a3_mhits > 1 |a4_mhits > 1 |a5_mhits > 1 |a6_mhits > 1 |a7_mhits > 1 |ab8_mhits > 1 | b1_mhits > 1 | b2_mhits > 1 | b3_mhits > 1 | b4_mhits > 1 | b6_mhits > 1 ) %>% select(fbna,chrom,posit,grep(x=names(res), pattern = "*_[mp]hits_*", value=T)) #848 / 562 when mhits= "minimum"

observ2_1 <- res %>% rowwise() %>% filter(a1_mhits > 1 | a2_mhits > 1 |a3_mhits > 1 |a4_mhits > 1 |a5_mhits > 1 |a6_mhits > 1 |a7_mhits > 1 |ab8_mhits > 1 | b1_mhits > 1 | b2_mhits > 1 | b3_mhits > 1 | b4_mhits > 1 | b6_mhits > 1 ) %>% mutate(ckvar = var(c(a1_mhits , a2_mhits ,a3_mhits ,a4_mhits ,a5_mhits ,a6_mhits ,a7_mhits,ab8_mhits , b1_mhits , b2_mhits , b3_mhits , b4_mhits , b6_mhits), na.rm=TRUE)  ) %>% filter(ckvar > 0.01) %>% select(fbna,chrom,posit,ckvar,grep(x=names(res), pattern = "*_[mp]hits_*", value=T)) 











tmp1 <- a3a4_mhits_B1 %>% group_by( class,numL) %>% summarise(n()) # 255 groups
#write.csv(a3a4_mhits_B1,file="a3a4_mhits_B1.csv")

a3a4_mpmatch <- res %>% group_by(fbna) %>% filter(min(a3_hits_1,a3_hits_2,a3_hits_3,a3_hits_4,a3_hits_5,na.rm=T) == a3_phits & min(a4_hits_1,a4_hits_2,a4_hits_3,a4_hits_4,a4_hits_5,na.rm=T) == a4_phits & (a3_phits > 1 | a4_phits > 1)) %>% select(fbna,chrom,posit,class,numL,grep(x=names(res), pattern = "*_[mp]hits_*", value=T)) # 169 


a3a4_mpmatch_uniq <- a3a4_mpmatch %>% group_by(class, numL) %>% summarise(count=n())  #89 
write.csv(a3a4_mpmatch, file="to_google_drive.csv")
