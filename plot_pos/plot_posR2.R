library(ggplot2)

setwd("~/cloud/para/plot_pos/")

a3_dfs=list.files(path='./a3_48gene',full.names = T)
a4_dfs=list.files(path='./a4_48gene',full.names = T)
ref_dfs=list.files(path='./ref_48gene',full.names = T)

# x chrom together


for(i in 40:40) {
    print(i)
    test_a3 <- read.table(file = a3_dfs[i], header = F)
    test_a4 <- read.table(file = a4_dfs[i], header = F)
    test_ref <- read.table(file = ref_dfs[i], header = F)
    
    colnames(test_a3) <- c("chrom",'start','end','idt','lpr')
    colnames(test_a4) <- c("chrom","start","end",'idt','lpr')
    colnames(test_ref) <- c("chrom","start","end","idt","lpr")
    
    
    #jit3 <- rnorm( nrow(test_a3), mean = -5, sd = 5)
    #jit4 <- rnorm(nrow(test_a4), mean = -5, sd = 5)
    #jitf <- rnorm(nrow(test_ref), mean = -5, sd = 5)
    
    p3 <-   ggplot(test_a3) + 
        geom_segment(aes(x=start,xend=end,y=1:nrow(test_a3),yend=1:nrow(test_a3),color=abs(start-end)),size= 4,alpha=1) +  
        geom_text(aes(x=(start+end)/2,y=1:nrow(test_a3),label=round(lpr,2)),vjust=1.3,color='red') +
        geom_text(aes(x=(start+end)/2,y=1:nrow(test_a3),label=round(idt,2)),color='green') +
        #ylim(c(0,20)) +
        ggtitle(a3_dfs[i])
    
    p4 <-   ggplot(test_a4) + 
        geom_segment(aes(x=start,xend=end,y=1:nrow(test_a4),yend=1:nrow(test_a4),color=abs(start-end)),size= 4,alpha=1) +  
        geom_text(aes(x=(start+end)/2,y=1:nrow(test_a4),label=round(lpr,2)),vjust=1.3,color='red') +
        geom_text(aes(x=(start+end)/2,y=1:nrow(test_a4),label=round(idt,2)),color='green') +
        #geom_text(aes(x=(start+end)/2,y=idt,label=abs(start-end)),hjust=5,color='yellow') + 
        #ylim(c(0,20)) +
        ggtitle(a4_dfs[i])
    
    pf <-   ggplot(test_ref) + 
        geom_segment(aes(x=start,xend=end,y=1:nrow(test_ref),yend=1:nrow(test_ref),color=abs(start-end)),size= 4,alpha=1) +  
        geom_text(aes(x=(start+end)/2,y=1:nrow(test_ref),label=round(lpr,2)),vjust=1.3,color='red') +
        geom_text(aes(x=(start+end)/2,y=1:nrow(test_ref),label=round(idt,2)),color='green') +
        #geom_text(aes(x=(start+end)/2,y=idt,label=abs(start-end)),hjust=5,color='yellow') + 
        #ylim(c(0,20)) +
        ggtitle(ref_dfs[i])  
    
    grid.arrange(p3,p4,pf,nrow=3)
    
    #  assign(x=paste0("plot",i),value = p)
    
}

