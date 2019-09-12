library(ggplot2)
library(ggrepel)

library(dplyr)
library(gridExtra)
library(reshape2)
library(grid)
library(gridBase)

library(ggpubr)

setwd("~/cloud/para/for_slides/gene_pos")
#twoL <- read.table('slide_2L.txt', header=T)
#twoL$fill <- as.character(twoL$fill)
twoR <- read.table('slide_2R.txt', header=T)
twoR$fill <- as.character(twoR$fill)
#threeL <- read.table('3L.txt', header=T)
#threeL$fill <- as.character(threeL$fill)
threeR <- read.table('slide_3R.txt', header=T)
threeR$fill <- as.character(threeR$fill)
#X <- read.table('X.txt', header=T)
#X$fill <- as.character(X$fill)


plotonechrome <- function(xx,tx,tg=FALSE) {
ggplot(data=xx) + 
    geom_point(aes(x=position,y=strain,color=fill),position=position_jitter(width=0,height=0.001),size=3) +
    geom_abline(slope = 0,intercept = 1, alpha=0.4) +
    geom_abline(slope = 0,intercept = 2, alpha=0.5) +
    geom_abline(slope = 0,intercept = 3, alpha=0.4) +
    scale_color_manual(values = c('1'='gray60','2'='darkcyan','3'='hotpink2','5'='green','7'='black','to other chromosome'='orange'))+ 
    #geom_text(aes(x=position,y=strain,label=geneName, color=fill), vjust=-1.8,position=position_jitter(width=0,height=0.3)) + 
    geom_text_repel(aes(x = position, y = strain, label = geneName,color=fill), data=xx[xx$label==TRUE,],size=4) + 
    labs(title=paste("chromosome arm",tx) , y="strain", x="relative position %", color="copy number") +
    theme_bw() +
    #theme(axis.text.x=element_text(angle=0, hjust=1)) +
    #coord_cartesian(xlim=c(1,20), ylim=c(0, 0.8)) +
    #scale_x_discrete(labels=1:20) + 
    scale_x_continuous(breaks =seq(0,100,by=10), limits = c(0,100)) + 
    scale_y_continuous(breaks = c(1,2,3),labels = c('A3','iso1','A4')) + 
    theme(
        panel.grid = element_line(color="grey85"),  
        legend.position = tg,
        legend.background =  element_rect( size=0.5, linetype="solid",colour ="darkblue"),
        legend.title = element_text(face='bold', size = 15),
        legend.text = element_text(face='plain', size = 15),
        plot.title = element_text(size = 16, face = "bold", hjust=0.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.2,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"),
        #axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
        #axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))
}


plotonechrome_withx <- function(xx,tx,tg=FALSE) {
    ggplot(data=xx) + 
        geom_point(aes(x=position,y=strain,color=fill),position=position_jitter(width=0,height=0.001),size=3) +
        geom_abline(slope = 0,intercept = 1, alpha=0.4) +
        geom_abline(slope = 0,intercept = 2, alpha=0.5) +
        geom_abline(slope = 0,intercept = 3, alpha=0.4) +
        scale_color_manual(values = c('1'='gray60','2'='darkcyan','3'='hotpink2','5'='green','7'='black','to other chromosome'='orange'))+ 
        #geom_text(aes(x=position,y=strain,label=geneName, color=fill), vjust=-1.8,position=position_jitter(width=0,height=0.3)) + 
        geom_text_repel(aes(x = position, y = strain, label = geneName,color=fill), data=xx[xx$label==TRUE,],size=4) + 
        labs(title=paste("chromosome arm",tx) , y="strain", x="relative position %") +
        theme_bw() +
        #theme(axis.text.x=element_text(angle=0, hjust=1)) +
        #coord_cartesian(xlim=c(1,20), ylim=c(0, 0.8)) +
        #scale_x_discrete(labels=1:20) + 
        scale_x_continuous(breaks =seq(0,100,by=10), limits = c(0,100)) + 
        scale_y_continuous(breaks = c(1,2,3),labels = c('A3','iso1','A4')) + 
        theme(
            panel.grid = element_line(color="grey85"),  
            legend.position = tg,
            legend.background =  element_rect( size=0.5, linetype="solid",colour ="darkblue"),
            legend.title = element_blank(),
            legend.text = element_text(face='bold', size = 18),
            plot.title = element_text(size = 13, face = "bold",hjust=.05),
            axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
            axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
            axis.title.x = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"), 
            #axis.title.x = element_blank(),
            axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=.5,face="bold"))
}
#1200 * 400 size # lock for 8 in slides
#p1 <- plotonechrome(twoL,"2L",tg="none")
plotonechrome(twoR,"2R",tg="none")
#p3 <- plotonechrome(threeL,"3L",tg="none")
plotonechrome(threeR,"3R",tg="none")
#p5 <- plotonechrome_withx(X,"X",tg="none")
