library(ggplot2)
library(ggrepel)

setwd("~/cloud/para/for_slides/gene_expression")
sig_df <- read.table('s_para_log_sig.txt', header=T)


ggplot(data = sig_df) + 
   # geom_point(aes(x=log2( 1 + round((count_n1 + count_n2)/2 - abs(count_n1-count_n2)/2)), y=log2( 1 + round((count_n1+  count_n2)/2 + abs(count_n1-count_n2)/2 )), color=Gene, shape=condition), size = 3) + 
    geom_point(aes(x=log2(count_n1), y=log2(count_n2), color=Gene), size=4) + 
    #geom_text_repel(aes(log2( 1 + round((count_n1 + count_n2)/2 - abs(count_n1-count_n2)/2)), y =log2( 1 + round((count_n1+  count_n2)/2 + abs(count_n1-count_n2)/2 )) , label = Gene,color=Gene)) + 
    geom_text_repel(aes(x=log2(count_n1), y= log2(count_n2), label = Gene,color=Gene), size=5) + 
    #scale_color_brewer() +
    #scale_shape_manual(values=c(16,2)) + 
    coord_fixed(ratio=1, xlim=c(0,16), ylim=c(0,16)) + 
    geom_abline(slope = 1, intercept = 0) + 
    labs(x="log2(read counts) of copy 1" , y="log2(read counts) of copy 2", title ="paralog-specific read counts") +
    guides(color=FALSE) + 
    theme_bw() +
    theme(
        legend.position = c(0.8,0.1),
        legend.title = element_text(size=15),
        legend.text=element_text(size=13),
        panel.grid = element_line(color="grey85"),  
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        #axis.title = element_blank(),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))

ggsave(file="~/cloud/para/gene_expression/fig4_paralog_point.pdf", figure4, device = "pdf", width = 15, height = 15, units = "cm",dpi=200)

