library(ggplot2)
library(ggrepel)

setwd("~/cloud/para/gene_expression")
tp <- read.table('total_point.txt', header=T)
tp$copy_max <- as.character(tp$copy_max)

fp5 <- ggplot(data = tp) + 
     geom_point(aes(x=log2( 1 + (count_a3 + count_a4)/2 - abs(count_a3-count_a4)/2), y=log2( 1 + (count_a3+  count_a4)/2 + abs(count_a3-count_a4)/2 ),shape=condition), size = 1.8) + 
    scale_shape_manual(values=c(16,17)) +  
    #scale_fill_manual(values = c('1'='gray60','2'='darkcyan','3'='hotpink2','5'='green','7'='black','9'='orange')) +
    #geom_point(aes(x=log2(count_a3), y=log2(count_a4), color=Gene, shape=condition), size=3) + 
    geom_text_repel(aes(log2( 1 + (count_a3 + count_a4)/2 - abs(count_a3-count_a4)/2), y =log2( 1 + (count_a3+  count_a4)/2 + abs(count_a3-count_a4)/2 ) , label = geneName,color=geneName), size= 2.5) + 
    #geom_text_repel(aes(x=log2(count_a3), y= log2(count_a4), label = Gene, color= Gene)) + 
    coord_fixed(ratio=1, xlim=c(0,20), ylim=c(0,20)) + 
    geom_abline(slope=1, intercept = 0, linetype='dashed') +
    geom_abline(aes(slope = slope, intercept = int)) + 
    #geom_abline(slope=1, intercept = log2(5)) + 
    #geom_abline(slope=1, intercept = log2(7)) + 
    labs(x="log2(read counts) of lower expression" , y="log2(read counts) of higher expression", title ="total expression read counts in two conditions") +
    guides(color=FALSE, fill=FALSE) + 
    theme_bw() +
    theme(
        #legend.position = c(0.8,0.1),
        legend.title = element_text(size=15),
        legend.text=element_text(size=13),
        panel.grid = element_line(color="grey85"),  
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        #axis.title = element_blank(),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold")) + 
    facet_wrap(~copy_max)

