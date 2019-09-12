library(ggplot2)
library(ggrepel)

setwd("~/cloud/para/for_slides/gene_expression")
exp35 <- read.table('s_35_exprssion.txt', header=T)
exp35$copy <- as.character(exp35$copy)

plotonegene <- function(x) {
ggplot(data=exp35) + 
    geom_bar(aes(x=strain, y = fpkm, fill=copy), stat='identity', color="black",position = "dodge",data=exp35[exp35$geneName == x,]) +
    scale_fill_manual(values = c('1'='gray60','2'='darkcyan','3'='hotpink2','5'='green','7'='black','to other chromosome'='orange'))+
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
        plot.title = element_text(size = 18, face = "bold",hjust=.5),
        axis.text.x = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=13,hjust=.5,vjust=.5,face="plain"),
        axis.title.y = element_text(colour="grey20",size=15,hjust=.5,vjust=0,face="bold"),
        axis.title.x = element_blank()
        )
}
plotonegene('Gs1')
plotonegene('cyp28d1')
plotonegene('Prosβ5R2')

plotonegene('Haspin')
plotonegene('IntS3')
plotonegene('photorepair')
plotonegene('θTrypsin')
plotonegene('CG6472')
plotonegene('Lapsyn')
plotonegene('Glycogenin')
plotonegene('CG44245')s
plotonegene('Ugt49B1')
plotonegene('Ugt49C1')
plotonegene('CG17658')

plotonegene('spook')
plotonegene('Alp9')
plotonegene('Snakeskin')
plotonegene('mRpS5')

plotonegene('Or85f')
plotonegene('CG34117')

plotonegene('Ugt86Dh')
plotonegene('DNAlig3')

plotonegene('CG31157')
plotonegene('CG7966')

plotonegene('CG6912')
plotonegene('CG6125')

plotonegene('Dlc90F')
plotonegene('CG18600')

plotonegene('Ugt303B3')
plotonegene('CG1894')

plotonegene('CG2233')
plotonegene('fiz')
plotonegene('CG15044')
plotonegene('CG15043')
plotonegene('Mgstl')
