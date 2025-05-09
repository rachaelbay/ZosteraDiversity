setwd("~/Documents/ZosteraGenomics/DivEffects/")

set.seed(12345)

library(tidyverse)
library(patchwork)
library(viridis)
library(randomForest)

sum <- readRDS("Abbott Results/Abbott.Summary.rds")

init.mean <- readRDS("Abbott Results/Abbott.Init.mean.het.var.rds")
fin.mean <- readRDS("Abbott Results/Abbott.Final.mean.het.var.rds")
rich <- read.csv("StachData/Abbott et al. 2017 Ecology/survival_bio_all.csv")
trait <- read.csv("StachData/Abbott et al. 2017 Ecology/Chapter3_data.csv")


merge.a <- sum %>% 
          left_join(init.mean,by="Plot") %>%
          left_join(fin.mean,by="Plot") %>%
          left_join(rich[,c("Plot","number")],by="Plot")

merge <- left_join(merge.a,trait[,c("Plot","Rao_intitial",
                                    "Rao_final_weighted","Ave_R_Intitial","R_weighted_final",
                                    "genotypic_diversity","genotypic_evenness")],by="Plot")


## Relatedness
i.rel <- ggplot(merge, aes(x=Ave_R_Intitial,y=total,color=as.factor(tot))) + geom_point(size=2) +
  theme_bw() + ylab("Total Biomass") + xlab("Relatedness") + 
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets"))
f.rel <- ggplot(merge, aes(x=R_weighted_final,y=total,color=as.factor(tot))) + geom_point(size=2) +
  theme_bw() + ylab("Total Biomass") + xlab("Relatedness") + 
  stat_smooth(method="lm",color="grey40") +
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets"))

## Heterozygosity
i.het <- ggplot(merge, aes(x=het.x,y=total,color=as.factor(tot))) + geom_point(size=2) +
  theme_bw() + ylab("Total Biomass") + xlab("Heterozygosity") + 
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets"))
f.het <- ggplot(merge, aes(x=het.y,y=total,color=as.factor(tot))) + geom_point(size=2) +
  theme_bw() + ylab("Total Biomass") + xlab("Heterozygosity") + 
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets")) +
  theme(plot.title = element_text(hjust = 0.5))

i.var <- ggplot(merge, aes(x=var.x,y=total,color=as.factor(tot))) + geom_point(size=2) +
  theme_bw() + ylab("Total Biomass") + xlab("Allelic variance") + 
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets")) +
  ggtitle("Initial") +  theme(plot.title = element_text(hjust = 0.5))
f.var <- ggplot(merge, aes(x=var.y,y=total)) + geom_point(aes(color=as.factor(tot)),size=2) +
  theme_bw() + ylab("Total Biomass") + xlab("Allelic variance") + 
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets")) +
  stat_smooth(method="lm",color="grey40") +
  ggtitle("Final") + theme(plot.title = element_text(hjust = 0.5))

## Genotype Richness
i.gen <- ggplot(merge, aes(x=as.factor(tot),y=total,color=as.factor(tot))) + 
  geom_boxplot(show.legend=F) + geom_jitter(size=2,width=0.1) +
  theme_bw() + ylab("Total Biomass") + xlab("Genet richness") +
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets"))
f.gen <- ggplot(merge, aes(x=number,y=total)) +
  theme_bw() + ylab("Total Biomass") + xlab("Genet richness") +
  geom_point(size=2,aes(col=as.factor(tot))) + stat_smooth(method="lm",color="grey40") +
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets"))


plot <- i.var + f.var + 
  i.rel + f.rel +
  i.gen + f.gen +
  i.het + f.het + 
  plot_annotation(tag_levels='A')  +plot_layout(guides="collect",nrow=4) 


##Supplementary figure
pdf("Figures/GenomeWide_variation_biomass_04.14.25.pdf",height=10,width=6)
plot &   theme(legend.position="bottom") &
  theme(legend.title=element_blank())
dev.off()

#Stats
summary(lm(total~var.x,data=merge))
summary(lm(total~var.y,data=merge))
summary(lm(total~Ave_R_Intitial,data=merge))
summary(lm(total~R_weighted_final,data=merge))
summary(lm(total~het.y,data=merge))
summary(lm(total~het.y,data=merge))
summary(lm(total~number,data=merge))
t.test(total~as.factor(tot),data=merge)



##Random Forest
library(randomForest)
sub.com <- merge[,c("total","number","het.y","var.y",
                    "Rao_final_weighted","R_weighted_final","genotypic_evenness")]
names(sub.com) <- c("Final_Biomass","Genet_Richness","Heterozygosity","Allelic_Variance","Trait_Diversity","Relatedness","Genet_Evenness")
dim(na.omit(sub.com))
com.rf <- randomForest(Final_Biomass~., data=na.omit(sub.com),importance=T,ntree=1000,nPerm=100)
com.rf
varImpPlot(com.rf,type=1,pch=16,pt.cex=1.5,main="",n.var=6)
sub.com.imp <- data.frame(importance(com.rf))

rf.fig <- ggplot(sub.com.imp,aes(x=X.IncMSE,y=reorder(rownames(sub.com.imp),X.IncMSE))) + 
  geom_dotplot(dotsize=1.5,stackdir="center") + theme_bw() + 
  ylab("") + xlab("Mean Decrease in Accuracy")



### Change in allelic variance
merge$delta.var <- merge$var.y-merge$var.x

ggplot(merge, aes(x=delta.var,y=var.x)) + geom_point(aes(color=as.factor(tot)),size=2)

## Delta variance v. biomass
d.var <- ggplot(merge, aes(x=delta.var,y=total)) + geom_point(aes(color=as.factor(tot)),size=2) +
  theme_bw() + ylab("Total Biomass") + xlab(expression(Delta~ "Allelic variance")) + 
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets")) +
  stat_smooth(method="lm",color="grey40") +
  theme(plot.title = element_text(hjust = 0.5))

## Final allelic variance vs. biomass with no title
f.var.notitle <- ggplot(merge, aes(x=var.y,y=total)) + geom_point(aes(color=as.factor(tot)),size=2) +
  theme_bw() + ylab("Total Biomass") + xlab("Allelic variance") + 
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets")) +
  stat_smooth(method="lm",color="grey40") +
  theme(plot.title = element_text(hjust = 0.5))


summary(lm(delta.var~total,data=merge))

pdf("Figures/Figure1_04.14.25.pdf",height=3,width=11)
rf.fig + f.var.notitle + d.var +  + plot_annotation(tag_levels='A') & 
  theme(legend.title=element_blank())  &
  theme(plot.margin = unit(c(0.1,0.5,0.1,0.5),"cm"))
dev.off()
