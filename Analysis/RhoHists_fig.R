setwd("~/Documents/ZosteraGenomics/DivEffects/")

library(tidyverse)
library(patchwork)
library(viridis)
library(fastman)

i.rho <- readRDS("Abbott Results/Abbott.Init.Combined.Rho.rds")
f.rho <- readRDS("Abbott Results/Abbott.Final.Combined.Rho.rds")

rho <- left_join(i.rho,f.rho,by=c("Chr","Pos"))
saveRDS(rho,"Abbott Results/Abbott.All.Combined.Rho.rds")

ir <- ggplot(rho,aes(x=Rho.x)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.3,alpha=0.2,option="B")) + 
  geom_vline(aes(xintercept=0),linetype="dashed",size=1) +
  theme_bw() + xlab("Rho") + ylab("Count") + ggtitle("Initial") +
  theme(plot.title = element_text(hjust = 0.5))
fr <- ggplot(rho,aes(x=Rho.y)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.5,alpha=0.2,option="B")) + 
  geom_vline(aes(xintercept=0),linetype="dashed",size=1) +
  theme_bw() + xlab("Rho") + ylab("Count") + ggtitle("Final") +
  theme(plot.title = element_text(hjust = 0.5))
ip <- ggplot(rho,aes(x=pval.x)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.3,alpha=0.2,option="B")) + 
  theme_bw() + xlab("p-value") + ylab("Count")
fp <- ggplot(rho,aes(x=pval.y)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.5,alpha=0.2,option="B")) + 
  theme_bw() + xlab("p-value") + ylab("Count")



plot_grid(ir,fr,ip,fp,
          nrow=2)

#######
##Manhattan plot
#######

rho.man <- rho
rho.man$Chr <- as.numeric(rho.man$Chr)
rho.man$snp <- paste("snp",1:nrow(rho),sep="")
rho.man$padj <- p.adjust(rho.man$pval.y,method="fdr")

manhattan(rho.man,chr="Chr",bp="Pos",p="pval.y",snp="snp")

write.table(rho[,c("Chr","Pos")],"All_SNPs.txt",
            row.names=F,col.names=F,quote=F)

cols <- viridis(2,begin=0.3,end=0.5,option="B",alpha=0.6)

m <- fastman_gg(rho.man,chr="Chr",bp="Pos",p="pval.y",speedup=T,col=cols,
           suggestiveline = NA,genomewideline=-log10(1e-5),ylim=c(0,15),xlab="",cex=2) +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + theme_bw() + ylim(0,7) +
  theme(legend.position="none") + xlab("Chromosome")

layout <- "
AB
CD
EE
"

##Figure 2
pdf("Figures/Rho_dists.pdf",height=7,width=6)
ir+fr+ip+fp +m +plot_layout(guides="collect",design=layout)+plot_annotation(tag_levels='A')
dev.off()


###Supplementary figure for six plot only
i.rho <- readRDS("Abbott Results/Abbott.Init.Six.Rho.rds")
f.rho <- readRDS("Abbott Results/Abbott.Final.Six.Rho.rds")

rho <- left_join(i.rho,f.rho,by=c("Chr","Pos"))

ir <- ggplot(rho,aes(x=Rho.x)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.3,alpha=0.2,option="B")) + 
  geom_vline(aes(xintercept=0),linetype="dashed",size=1) +
  theme_bw() + xlab("Rho") + ylab("Count") + ggtitle("Initial") +
  theme(plot.title = element_text(hjust = 0.5))
fr <- ggplot(rho,aes(x=Rho.y)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.5,alpha=0.2,option="B")) + 
  geom_vline(aes(xintercept=0),linetype="dashed",size=1) +
  theme_bw() + xlab("Rho") + ylab("Count") + ggtitle("Final") +
  theme(plot.title = element_text(hjust = 0.5))
ip <- ggplot(rho,aes(x=pval.x)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.3,alpha=0.2,option="B")) + 
  theme_bw() + xlab("p-value") + ylab("Count")
fp <- ggplot(rho,aes(x=pval.y)) + 
  geom_histogram(color="grey30",fill=viridis(1,begin=0.5,alpha=0.2,option="B")) + 
  theme_bw() + xlab("p-value") + ylab("Count")

pdf("Figures/Rho_dists_Six.pdf",height=6,width=6)
ir+fr+ip+fp +plot_annotation(tag_levels='A')
dev.off()


