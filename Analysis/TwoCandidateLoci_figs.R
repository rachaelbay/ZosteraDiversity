setwd("~/Documents/ZosteraGenomics/DivEffects/")

library(tidyverse)
library(nlme)
library(lmtest)
library(viridis)
library(patchwork)
library(ggeffects)

###Abbott 2017
Ab18 <- readRDS("Abbott Results/Abbott2017_plot_2loc.rds")
Ab18$fail <- is.na(Ab18$above)
Ab18$value=1

b1 <- ggplot(Ab18, aes(x=inv1.tf,y=total)) + geom_boxplot() + 
  ylab("Total Biomass") + xlab("Loc1 polymorphic") + 
  geom_point(aes(color=as.factor(Ab18$tot))) + theme_bw() +
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets"))

b2 <- ggplot(Ab18, aes(x=inv2.tf,y=total)) + geom_boxplot() +
  ylab("Total Biomass") + xlab("Loc2 polymorphic") + 
  geom_point(aes(color=as.factor(Ab18$tot))) + theme_bw() +
  scale_color_viridis(discrete=T,end=0.7,begin=0.3,labels=c("2 genets","6 genets"))

f1 <- ggplot(Ab18, aes(fill=fail,x=inv1.tf,y=value)) + 
  geom_bar(position="fill",sta="identity") + 
  ylab("Proportion of plots") + xlab("Loc1 polymorphic") +
  scale_fill_viridis(discrete=T,begin=0.2,end=0.8,option="F",labels=c("Succeed","Fail")) +
  theme_bw()
f2 <- ggplot(Ab18, aes(fill=fail,x=inv2.tf,y=value)) + 
  geom_bar(position="fill",sta="identity") + 
  ylab("Proportion of plots") + xlab("Loc2 polymorphic") +
  scale_fill_viridis(discrete=T,begin=0.2,end=0.8,option="F",labels=c("Succeed","Fail")) +
  theme_bw()

###Plot

row1 <- (f1 + f2) + plot_layout(guides = "collect") & 
  theme(legend.position = "right") &
  theme(legend.title=element_blank())

row2 <- (b1 + b2) + plot_layout(guides = "collect") & 
  theme(legend.position = "right") &
  theme(legend.title=element_blank())

plot <- row1/row2

pdf("Figures/2loc_all.pdf",width=6,height=8)
plot +plot_annotation(tag_levels='A') 
dev.off()
