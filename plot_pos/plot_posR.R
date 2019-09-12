library(ggplot2)

setwd("~/cloud/para/plot_pos/")

a3_dfs=list.files(path='./a3_48gene',full.names = T)
a4_dfs=list.files(path='./a4_48gene',full.names = T)
ref_dfs=list.files(path='./ref_48gene',full.names = T)

for(i in 17:17) {
    print(i)
    test_a3 <- read.table(file = a3_dfs[i], header = F)
    test_a4 <- read.table(file = a4_dfs[i], header = F)
    test_ref <- read.table(file = ref_dfs[i], header = F)
    
    colnames(test_a3) <- c("chrom",'start','end','idt','lpr')
    colnames(test_a4) <- c("chrom","start","end",'idt','lpr')
    colnames(test_ref) <- c("chrom","start","end","idt","lpr")
    

    jit3 <- rnorm( nrow(test_a3), mean = -5, sd = 5)
    jit4 <- rnorm(nrow(test_a4), mean = -5, sd = 5)
    jitf <- rnorm(nrow(test_ref), mean = -5, sd = 5)
    
 p3 <-   ggplot(test_a3) + 
        geom_segment(aes(x=start,xend=end,y=idt+jit3,yend=idt+jit3,color=abs(start-end)),size= 4,alpha=0.5) +  
        geom_text(aes(x=(start+end)/2,y=idt+jit3,label=lpr),vjust=-0.8,color='red') +
        geom_text(aes(x=(start+end)/2,y=idt+jit3,label=idt),color='green') +
        #geom_text(aes(x=(start+end)/2,y=idt,label=abs(start-end)),hjust=5,color='yellow') + 
        ylim(c(50,110)) +
        ggtitle(a3_dfs[i])
    
 p4 <-   ggplot(test_a4) + 
     geom_segment(aes(x=start,xend=end,y=idt+jit4,yend=idt+jit4,color=abs(start-end)),size= 4,alpha=0.5) +  
     geom_text(aes(x=(start+end)/2,y=idt+jit4,label=lpr),vjust=-0.8,color='red') +
     geom_text(aes(x=(start+end)/2,y=idt+jit4,label=idt),color='green') +
     #geom_text(aes(x=(start+end)/2,y=idt,label=abs(start-end)),hjust=5,color='yellow') + 
     ylim(c(50,110)) +
     ggtitle(a4_dfs[i])
 
 pf <-   ggplot(test_ref) + 
     geom_segment(aes(x=start,xend=end,y=idt+jitf,yend=idt+jitf,color=abs(start-end)),size= 4,alpha=0.5) +  
     geom_text(aes(x=(start+end)/2,y=idt+jitf,label=lpr),vjust=-0.8,color='red') +
     geom_text(aes(x=(start+end)/2,y=idt+jitf,label=idt),color='green') +
     #geom_text(aes(x=(start+end)/2,y=idt,label=abs(start-end)),hjust=5,color='yellow') + 
     ylim(c(50,110)) +
     ggtitle(ref_dfs[i])  
 
 grid.arrange(p3,p4,pf,nrow=3)
    
#  assign(x=paste0("plot",i),value = p)

}

