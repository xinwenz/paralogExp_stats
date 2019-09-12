library(ggplot2)
setwd("~/cloud/para/plot_pos/")
a3_dfs=list.files(path='./a3_48gene',full.names = T)
a4_dfs=list.files(path='./a4_48gene',full.names = T)
a4_dfs_nocds =list.files(path='./a4_48gene_nocds/', full.name=T)
ref_dfs=list.files(path='./ref_48gene',full.names = T)

plot3 <- function(v) {
    genes_a3 <- data.frame()
    genes_a4 <- data.frame()
    genes_ref <- data.frame()
    for (i in v) {
        tmp <- read.table(file=a3_dfs[i], header=F)
        tmp1 <- cbind(tmp,i)
        genes_a3 <- rbind(genes_a3,tmp1)
        tmp <- read.table(file=a4_dfs[i], header=F)
        tmp1 <- cbind(tmp,i)
        genes_a4 <- rbind(genes_a4,tmp1)
        tmp <- read.table(file=ref_dfs[i], header=F)
        tmp1 <- cbind(tmp,i)
        genes_ref <- rbind(genes_ref,tmp1)
    }
    colnames(genes_a3) <- c("chrom",'start','end','idt','lpr','labg')
    colnames(genes_a4) <- c("chrom",'start','end','idt','lpr','labg')   
    colnames(genes_ref) <- c("chrom",'start','end','idt','lpr','labg') 
    
  p3 <-  ggplot(genes_a3) + 
        geom_segment(aes(x=start,xend=end,y=1:nrow(genes_a3),yend=1:nrow(genes_a3),color=abs(start-end)),size= 4,alpha=1) +  
        geom_text(aes(x=(start+end)/2,y=1:nrow(genes_a3),label=round(lpr,2)),vjust=1.3,color='red') +
        geom_text(aes(x=max(end,start), y = 1:nrow(genes_a3), label=labg), color='purple',hjust=-1) +
        geom_text(aes(x=(start+end)/2,y=1:nrow(genes_a3),label=round(idt,2)),color='green') 
        #ylim(c(0,20)) +

    p4 <-  ggplot(genes_a4) + 
        geom_segment(aes(x=start,xend=end,y=1:nrow(genes_a4),yend=1:nrow(genes_a4),color=abs(start-end)),size= 4,alpha=1) +  
        geom_text(aes(x=(start+end)/2,y=1:nrow(genes_a4),label=round(lpr,2)),vjust=1.3,color='red') +
        geom_text(aes(x=max(start,end), y = 1:nrow(genes_a4), label=labg), color='purple',hjust=-1) +
     geom_text(aes(x=(start+end)/2,y=1:nrow(genes_a4),label=round(idt,2)),color='green') 
       # xlim(c(21100000,21250000)) 
    
    pf <-  ggplot(genes_ref) + 
        geom_segment(aes(x=start,xend=end,y=1:nrow(genes_ref),yend=1:nrow(genes_ref),color=abs(start-end)),size= 4,alpha=1) +  
        geom_text(aes(x=(start+end)/2,y=1:nrow(genes_ref),label=round(lpr,2)),vjust=1.3,color='red') +
        geom_text(aes(x=max(start,end), y = 1:nrow(genes_ref), label=labg), color='purple',hjust=-1) +
    geom_text(aes(x=(start+end)/2,y=1:nrow(genes_ref),label=round(idt,2)),color='green') 
    
    grid.arrange(p3,p4,pf,nrow=3)
}

# x chrom together 
plot3(c(19,20,21,22,17)) # x chrome
plot3(c(21,22)) 


# mito 
plot3(c(11,7,6,5,8,12,14,45,13,15,9,10)) # inversion in a4,  
plot3(c(11,7,6,5,8)) # inversion in a4,  
plot3(c(11,7)) 
plot3(c(6,5,8,12,14,45,13,15,9,10))  #  mito is circle, circle duplicate. 


#3R 

plot3(c(30,43,38,31,41,32,33,34,16,35,36,37))
plot3(c(30,43)) # a4, 2 fold 
plot3(c(41,32)) #a4 ,3 fold
plot3(c(33,34)) # not together 
plot3(c(16,35)) #a3, 2 fold 


# 3L 
plot3(c(3,28,29,39))
plot3(c(3,28)) # a3 , 2/3 fold 
plot3(c(39)) # a4, very far away 2 copy. 
# 




# 2R
plot3(c(25,47,46,18,26,27))
plot3(c(25)) # a4, 2 fold, far away 
plot3(c(25,47,18,26,27))
plot3(c(25,26,27)
plot3(27)


# 2L 
plot3(c(1,23,42))

######################################### only A4 ## one part 
genes_a4 <- data.frame()
for (i in c(25,47,46,18,26,27)) {
    tmp <- read.table(file=a4_dfs_nocds[i], header=F)
    tmp1 <- cbind(tmp,i)
    genes_a4 <- rbind(genes_a4,tmp1)
}

colnames(genes_a4) <- c("chrom",'start','end','idt','lpr','labg')   

ggplot(genes_a4) + 
    geom_segment(aes(x=start,xend=end,y=1:nrow(genes_a4),yend=1:nrow(genes_a4),color=abs(start-end)),size= 4,alpha=1) +
    #geom_text(aes(x=(start+end)/2,y=1:nrow(genes_a4),label=round(lpr,2)),vjust=1.3,color='red') +
    geom_text(aes(x=max(start,end), y = 1:nrow(genes_a4), label=labg), color='purple',hjust=-1) 
    #geom_text(aes(x=(start+end)/2,y=1:nrow(genes_a4),label=round(idt,2)),color='green') 
# xlim(c(21100000,21250000)) 
