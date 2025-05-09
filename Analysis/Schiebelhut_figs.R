setwd("~/Documents/ZosteraGenomics/DivEffects/")

library(tidyverse)
library(SeqArray)
library(SNPRelate)
library(sf)
library(raster)
library(viridis)
library(pheatmap)
library(ggplotify)
library(patchwork)

###################
### Map
###################

#Read in data and shapefiles
dat <- read.delim("../BodTom12.21/construct/Zm_TomBod_LatLon.txt")
Count.shp <- read_sf("../BodTom12.21/Map/California_County_Boundaries/California_County_Boundaries.shp")

#Crop and format data
BB <- extent(min(dat$Lon)-0.1,max(dat$Lon)+0.1, min(dat$Lat)-0.1,max(dat$Lat)+0.1)
Count.shp2 <- st_crop(Count.shp,BB)
Count.df <- fortify(Count.shp2)

#Plot
map <- ggplot() +
  theme_bw() +
  geom_sf(data = Count.df, fill = 'grey90')+
  geom_point(data = dat, aes(x=Lon, y=Lat,color=Bay), size=3.0)+
  labs(title = '', xlab = ' ', ylab = ' ') +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Longitude") + ylab("Latitude") +
  scale_color_viridis(discrete="T",option="F",begin=0.3,end=0.7)

map


######################
## PCA of candidates
######################
pca <- readRDS("Cands.pca.2loci.rds")

pcplot <- ggplot(pca, aes(x=PC1,y=PC2,color=Bay)) + geom_point(size=2) +
  theme_bw() +
  scale_color_viridis(discrete="T",option="F",begin=0.3,end=0.7) +  
  theme(legend.position="none") +
  theme(legend.title=element_blank())
pcplot


#####################
## Freq v. Temp
#####################
agg <- pca %>%
  group_by(Site,Bay,depth) %>%
  summarize(Temp=mean(Temp),
            Longitude=mean(Longitude),
            Latitude=mean(Latitude),
            Loc1=mean(inv1)/2,
            Loc2=mean(inv2)/2)

loc1 <- ggplot(agg, aes(x=Temp,y=Loc1,col=Bay)) +geom_point(size=2) +
  theme_bw() +
  scale_color_viridis(discrete="T",option="F",begin=0.3,end=0.7) +  
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  ylab("Allele Frequency Loc1") + xlab("Mean Temperature")
lm1 <- lm(Loc1~Temp, data=agg)
summary(lm1)
lm1.2 <- lm(Loc1~Temp,data=agg %>% filter(Site!="PB"))
summary(lm1.2)

loc2 <- ggplot(agg, aes(x=Temp,y=Loc2)) +geom_point(size=2,aes(color=Bay)) +
  theme_bw() +
  scale_color_viridis(discrete="T",option="F",begin=0.3,end=0.7) +  
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  ylab("Allele Frequency Loc2") + xlab("Mean Temperature") + stat_smooth(method="lm",color="grey30")
lm2 <- lm(Loc2~Temp, data=agg)
summary(lm2)


#######################
## LD of candidates
#######################
##Get metadata
meta <- read.csv("../BodTom12.21/ZM_meta.csv")

##Get genomic data 
genos <- meta$INDV
genofile <- seqOpen("Zm_AllBodTom.gds")

##Get candidate SNPs
cands <- read.delim("annot/Cands.pos",header=F,sep=" ")
cands$snp <- paste(cands$V1,cands$V2,sep=".")

##Pull candidates from genofile
samp.list <- seqGetData(genofile,"sample.id")
chr <- seqGetData(genofile,"chromosome")
pos <- seqGetData(genofile,"position")
snp.pos <- data.frame(Chr=paste("Chr",chr,sep=""),Pos=pos)
snp.pos$snp <- paste(snp.pos$Chr,snp.pos$Pos,sep=".")
cand.snp <-which(snp.pos$snp%in%cands$snp==T) 

ld <- snpgdsLDMat(genofile,snp.id=cand.snp,sample.id=genos,slide=-1,num.thread=2) 
ld.mat <- ld$LD^2

ld.mat[upper.tri(ld.mat)] <- NA
rownames(ld.mat) <- cands$snp
colnames(ld.mat) <- cands$snp
chr.groups <- data.frame(Chr=cands$V1)
rownames(chr.groups) <- cands$snp

table(chr.groups)
gaps <- c(6,11,20,27,33 )

heat <- as.ggplot(pheatmap(ld.mat,cluster_cols = F,cluster_rows = F,
                           annotation_row=chr.groups,annotation_col=chr.groups,
                           show_rownames=F,show_colnames=F,
                           gaps_row=gaps,gaps_col = gaps, border_color=NA,na_col="white",
                           scale="none",
                           annotation_names_row=F,annotation_names_col=F))


pdf(file="Figures/Cand_heatmap_raw.pdf",width=5,height=4.5)
heat
dev.off()

layout <- "
AB#
ACD
"

#plot - will have to put in the heatmap manually - wah wah
pdf(file="Figures/Schieblhut_figs_raw.pdf",width=10,height=8)
map + pcplot + loc1 + loc2 + plot_layout(design=layout) + plot_annotation(tag_levels='A')
dev.off()
  