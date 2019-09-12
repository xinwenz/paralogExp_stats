library(dplyr)
# cyp28d1 :FBgn0031689
# Ugt86Dh :	FBgn0040252
a3_gene <- read.table("suz_a3.res",header=F,fill=T)
a4_gene <- read.table("suz_a4.res",header=F,fill=T)

a3_gene[is.na(a3_gene)] <- 0
tmp <- paste0(c("bs",'idt','lth'),rep(1:5,each=3))
names(a3_gene) <- c('fbna','nhits',tmp)
a3_gene <- a3_gene %>% group_by(fbna) %>% 
    mutate(ndup= sum( c(lth1,lth2,lth3,lth4,lth5) > 0.9*lth1  & c(idt1,idt2,idt3,idt4,idt5) > 95 ) )

a4_gene[is.na(a4_gene)] <- 0
names(a4_gene) <- c('fbna','nhits',tmp)
a4_gene <- a4_gene %>% group_by(fbna) %>% 
    mutate(ndup= sum( c(lth1,lth2,lth3,lth4,lth5) > 0.9*lth1  & c(idt1,idt2,idt3,idt4,idt5) > 95 ) )


a3a4_gene <-  merge.data.frame(a3_gene,a4_gene,by=1) 

##################################################3
a3a4_nodp <- a3a4_gene %>% filter(ndup.x == 1 & ndup.y == 1)


useless <- c(grep(pattern=".*[2,3,4,5].x",names(a3a4_gene),value=T),grep(pattern=".*[3,4,5].y",names(a3a4_gene),value=T) )
a4_2copy <- a3a4_gene %>% group_by(fbna) %>%
    filter(ndup.x < 2 & ndup.y == 2 )   %>%
            select(-useless)
a4_2copy[a4_2copy$fbna=='FBgn0040252',]
a4_2copy[a4_2copy$fbna=='FBgn0031689',]
# rubics :  Identy all> 90, length can't have more than 1.5 fold difference. ####### 

useless <- c(grep(pattern=".*[3,4,5].x",names(a3a4_gene),value=T),grep(pattern=".*[2,3,4,5].y",names(a3a4_gene),value=T) )
a3_2copy <- a3a4_gene %>% group_by(fbna) %>%
    filter(ndup.x == 1 & ndup.y < 1 & 
               idt1.x > 90 & idt2.x > 90 & idt1.y > 90  )   %>%
    select(-useless)


useless <- c(grep(pattern=".*[3,4,5].x",names(a3a4_gene),value=T),grep(pattern=".*[3,4,5].y",names(a3a4_gene),value=T) )
a3a4_2copy <- a3a4_gene %>% group_by(fbna) %>%
    filter(ndup.x == 1 & ndup.y == 1 & 
               idt1.x > 90 & idt2.x > 90 & idt2.x > 90 & idt1.y > 90  )   %>%
    select(-useless)

######## a3=2,a4=1 
useless <- c(grep(pattern=".*[3,4,5].x",names(a3a4_gene),value=T),grep(pattern=".*[3,4,5].y",names(a3a4_gene),value=T) )
a3a4_2copy <- a3a4_gene %>% group_by(fbna) %>%
    filter(ndup.x == 1 & ndup.y == 1 & 
               idt1.x > 90 & idt2.x > 90 & idt2.x > 90 & idt1.y > 90  )   %>%
    select(-useless)

ans <- c()
for(k in 0:5) { 
print(k)
tmp0 <- a3a4_gene %>% filter(ndup.x == k & ndup.y == 0 ) 
ans <- c(ans,nrow(tmp0)) 
tmp1 <- a3a4_gene %>% filter(ndup.x == k & ndup.y == 1 )
ans <- c(ans,nrow(tmp1)) 
tmp2 <- a3a4_gene %>% filter(ndup.x == k & ndup.y == 2 ) 
ans <- c(ans,nrow(tmp2)) 
tmp3 <- a3a4_gene %>% filter(ndup.x == k & ndup.y == 3 ) 
ans <- c(ans,nrow(tmp3)) 
tmp4 <- a3a4_gene %>% filter(ndup.x == k & ndup.y == 4 ) 
ans <- c(ans,nrow(tmp4)) 
tmp5 <- a3a4_gene %>% filter(ndup.x == k & ndup.y == 5 ) 
ans <- c(ans,nrow(tmp5)) 
}
matrix(ans,nrow=6,byrow = T)


##### read in cds ####
a4_cds <- read.table("a4_suz_cds.res",header=F)
colnames(a4_cds) <- c("cds_nm","geneFB","lth","hits")
a4mhit <- a4_cds %>% group_by(geneFB) %>% summarise(cds_amt=n(),mhits=median(hits))
tmp <- as.data.frame(table(a4mhit$mhits))

a <- a4mhit[which(a4mhit$mhits == 1.5),]
merge.data.frame(a,a4_cds,by="geneFB",all=F) %>% arrange(geneFB,desc(lth))

a <- a4mhit[which(a4mhit$mhits == 2),]
b <- merge.data.frame(a,a4_cds,by="geneFB",all=F) %>% arrange(geneFB,desc(lth))
