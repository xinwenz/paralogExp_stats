library(ggplot2)
library(ggrepel)

library(dplyr)
library(gridExtra)
library(reshape2)
library(grid)
library(gridBase)

library(ggpubr)

setwd("~/cloud/para/for_slides/gene_expression")
pexp16 <- read.table('s_16_paralog.txt', header=T)
#exp35$copy <- as.character(exp35$copy)

plotpara <- function(x) {
    ggplot(data=pexp16) + 
        geom_bar(aes(x=paralog, y = ncount, fill=strain), stat='identity', color="black",position = "dodge",data=pexp16[pexp16$Gene == x,]) +
        scale_fill_manual(values = c('A3'='grey','A4'='grey'))+
        #scale_alpha_manual(values = c('nic'=0.8, 'ctl'=0.2)) +
        labs(title=x, y="FPKM", x="relative position %") + 
        theme_bw() +
        #theme(axis.text.x=element_text(angle=0, hjust=1)) +
        #coord_cartesian(xlim=c(1,20), ylim=c(0, 0.8)) +
        #scale_x_continuous(breaks =seq(0,100,by=10), limits = c(0,100)) + 
        #scale_y_continuous(breaks = c(1,2,3),labels = c('A3','iso1','A4')) + 
        theme(
            #panel.grid = element_line(color="grey85"),  
            legend.position = 'none',
            #legend.background =  element_rect( size=0.5, linetype="solid",colour ="darkblue"),
            #legend.title = element_blank(),
            #legend.text = element_blank(),
            #plot.title = element_text(size = 13, face = "bold",hjust=.05),
            plot.title = element_text(size = 18, face = "bold",hjust=.5),
            axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
            axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"),
            axis.title.x = element_blank()
        )
}

#plotpara('Gs1')
p1 <-plotpara('cyp28d1')
#plotpara('Prosβ5R2')

p2 <- plotpara('Haspin')
p3 <- plotpara('IntS3')

p4 <- plotpara('CG6472')

p5 <- plotpara('Lapsyn')
p6 <- plotpara('Glycogenin')
p7 <- plotpara('CG44245')

p8 <- plotpara('spook')

p9 <- plotpara('Snakeskin')
p10 <- plotpara('mRpS5')

p11 <- plotpara('Ugt86Dh')

p12 <- plotpara('CG6912')
p13 <- plotpara('CG6125')


p14 <- plotpara('CG2233')
p15 <- plotpara('fiz')
#plotpara('CG15044')
p16 <- plotpara('CG15043')
#plotpara('Mgstl')

tx1=textGrob("FPKM of the genes", gp=gpar(fontface="bold", fontsize=15))

pa <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,nrow=4, top=tx1)

ggsave(file="~/cloud/para/gene_expression/fig4_paralog_fpkm.pdf",pa, device = "pdf", width = 20, height = 22.5, units = "cm",dpi=200)
