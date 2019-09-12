library(ggplot2)
library(ggrepel)
library(dplyr)
library(plyr)
library(ggpubr)
library(reshape2)

setwd("~/cloud/para/for_slides/gene_expression")
tp <- read.table('s_35_exprssion.txt', header=T)
copyMax <- read.table('s_total_point.txt', header=T)
copyMax <- copyMax %>% select(geneName, copy_max, slope, int) 


#tp$copy_max <- as.character(tp$copy_max)
tp1 <- tp %>% select(geneName,copy,depth_n_total)
tp2 <- dcast(tp1,value.var = "depth_n_total",id.vars="geneName", geneName~copy)
colnames(tp2) <- c("geneName",'copy1','copy2','copy3','copy5','copy7')
tp3 <- tp2 %>% group_by(geneName) %>% mutate(copyn = sum(copy2,copy3,copy5,copy7,na.rm = T) ) 

#join(copyMax,tp3,by='geneName')

tp4 <- join(copyMax,tp3,by='geneName') %>% mutate(lb = as.character(1:nrow(tp3))) %>% filter(geneName != 'Ugt303B3')  

ggplot(data = tp4) + 
     geom_point(aes(x=log2(copy1), y=log2(copyn), color=geneName), size = 1.9) + 
    #scale_shape_manual(values=c(16,17)) +  
    #scale_fill_manual(values = c('1'='gray60','2'='darkcyan','3'='hotpink2','5'='green','7'='black','9'='orange')) +
    #geom_point(aes(x=log2(count_a3), y=log2(count_a4), color=Gene, shape=condition), size=3) + 
    #geom_text_repel(aes(x=copy1, y = copyn , label = geneName,color=geneName), size= 4) + 
    geom_text_repel(aes(x=log2(copy1), y= log2(copyn), label = lb, color= geneName)) + 
    coord_fixed(ratio=1, xlim=c(0,20), ylim=c(0,20)) + 
    geom_abline(slope=1, intercept = 0, linetype='dashed') +
    geom_abline(aes(slope = slope, intercept = int)) + 
    #geom_abline(slope=1, intercept = log2(5)) + 
    #geom_abline(slope=1, intercept = log2(7)) + 
    labs(x="log2(read counts of singleton)" , y="log2(read counts of n copy)", title ="total expression read counts ") +
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


tp5 <- tp4 %>% mutate(ollr = log2(copyn/copy1) - log2(copy_max/1 )) 
tp5$copy_max <- as.character(tp5$copy_max)
## copy and pasted into excel, human force, add labels ### 
tp6 <- read.table('tp5_manu_yheight.txt', header=T)
tp6$copy_max <- as.character(tp6$copy_max)

ans <- ggplot(tp6) + 
    geom_histogram( aes(x= ollr, fill=copy_max, color=geneName), bins = 20) + 
    scale_y_continuous(breaks =0:15, labels=0:15, limits=c(0,8) )+ 
    scale_x_continuous(breaks =seq(-2.5,5.0,by=0.5), labels=seq(-2.5,5.0,by=0.5))+
    scale_color_manual(values=rep("white",34),guide='none') +  
    scale_fill_manual(values = c('1'='gray60','2'='darkcyan','3'='hotpink2','5'='green','7'='black','9'='orange'))  + 
    geom_text_repel(aes(x=ollr, y= yheight, label=geneName), color="blue")+
    labs(x="log2(expression:copyn/copy1) - log2(genome:copyn)" , y="count", title ="Distribution of log odds ratio", fill="copy number") +
    theme_bw() +
    theme(
        legend.position = c(0.8,0.8),
        legend.title = element_text(size=15),
        legend.text=element_text(size=13),
        panel.grid = element_line(color="grey85"),  
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        #axis.title = element_blank(),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold")) 
    #geom_text(aes(x=ollr, y= yhrandom, label=geneName))

ggsave(filename = '~/cloud/para/gene_expression/fig_hist.jpg',ans, device = 'jpg', width = 30, height = 20, units = "cm",dpi=200)






