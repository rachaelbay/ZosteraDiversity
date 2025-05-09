setwd("~/Documents/ZosteraGenomics/DivEffects/")

library(tidyverse)
library(SeqArray)
library(SNPRelate)

##Get metadata
meta <- read.csv("../BodTom12.21/ZM_meta.csv")

##Get genomic data 
genos <- meta$INDV
genofile <- seqOpen("Zm_AllBodTom.gds")
#samp.list <- seqGetData(genofile,"sample.id")
good.snp <- snpgdsSelectSNP(genofile,sample.id=genos,autosome.only=F,maf=0.1,missing.rate=0.1)
seqClose(genofile)

##Get candidate SNPs
cands <- read.delim("annot/Cands.pos",header=F,sep=" ")
cands$snp <- paste(cands$V1,cands$V2,sep=".")

#Find SNPs with at least 2 minor alleles
genofile <- seqOpen("Zm_AllBodTom.gds")
samp.list <- seqGetData(genofile,"sample.id")
chr <- seqGetData(genofile,"chromosome")
pos <- seqGetData(genofile,"position")
snp.pos <- data.frame(Chr=paste("Chr",chr,sep=""),Pos=pos)
snp.pos$snp <- paste(snp.pos$Chr,snp.pos$Pos,sep=".")
cand.snp <- which(snp.pos$snp%in%cands$snp==T) 

##PCA
cand.pca <- snpgdsPCA(genofile,sample.id=meta$INDV,snp.id=cand.snp)
samp.ids <- cand.pca$sample.id
cand.pca.frame <- data.frame(INDV=samp.ids,
                             PC1=cand.pca$eigenvect[,1],
                             PC2=cand.pca$eigenvect[,2],
                             PC3=cand.pca$eigenvect[,3],
                             PC4=cand.pca$eigenvect[,4])

plot(cand.pca$eigenvect,pch=".",cex=4)
abline(v=c(0,0.12),col="red")
abline(h=c(0.011,-0.085),col="red")

cand.pca.frame$inv1 <- 1
cand.pca.frame$inv1[cand.pca.frame$PC1<0] <- 0
cand.pca.frame$inv1[cand.pca.frame$PC1>0.12] <- 2

cand.pca.frame$inv2 <- 1
cand.pca.frame$inv2[cand.pca.frame$PC2>0.011] <- 0
cand.pca.frame$inv2[cand.pca.frame$PC2<(-0.085)] <- 2

#saveRDS(cand.pca.frame,"Cands.pca.2loci.rds")
cand.pca.frame <- readRDS("Cands.pca.2loci.rds")

