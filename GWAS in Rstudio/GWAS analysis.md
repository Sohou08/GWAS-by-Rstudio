
# Purpose: Genetic determinism of one trait through GWAS approach  #


```{r}
# In case ofseveral variables, possibility to choice one variable as folling:
i<-4
Y<-PG[,i]
liste<- which(!(is.na(Y)))
Y<-Y[liste]
hist(Y, main=names(PG)[i])
# the markers start at column 18
G<-PG[liste, 18:dim(PG)[2]]
length(Y)
```

#### Compute Kinship ####

```{r}
# genMatrix function of QTLRel make the compute of matrix
# 1000 markers
liste <- sample(colnames(G),1000,replace=FALSE)
K<-genMatrix(G[,liste])
I<-diag(length(Y))
hist(K$AA)
```

#### Compute the variance of each components ####

This model estimate the additive and environnmental variance.
```{r}
mod1<-estVC(y=Y,v=list(AA=K$AA,DD=NULL,HH=NULL,AD=NULL,MH=NULL,EE=I))
mod1

#Make loop for all markers selected
print(date())
GWAS.mod1<-scanOne(y=Y,gdat=G,vc=mod1,test="Chisq")
print(date())

#Store p-value
pval.mod1<-data.frame(marqueurs=colnames(G),pvalue=GWAS.mod1$p)
varnom<-names(PG)[i]
file<-paste("GRAINS_GWAS.mod1_",varnom,".Rdata", sep="") 
save(GWAS.mod1, file=file, compress= TRUE)
file<-paste("GRAINS.pval.mod1_",varnom,".Rdata", sep="")
save(pval.mod1, file=file, compress= TRUE)
```

#### Aftertreatment ####

```{r}
library("data.table")
library(ggplot2)
#http://www.bioconductor.org/packages/release/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
#https://github.com/jdstorey/qvalue
#install.packages("devtools")
library("devtools")
#install_github("jdstorey/qvalue")
library(qvalue)
```

#### P-value Data ####

```{r}
varnom <- "moy_prot"

file<-paste("GRAINS_GWAS.mod1_",varnom,".Rdata", sep="")
load(file)
pvalues<-GWAS.mod1$p
length(GWAS.mod1$p)
pval.mod1<-data.frame(Marker_ID=names(GWAS.mod1$p),pvalue=pvalues)
min(pval.mod1$pvalue)
#Distribution of p-values
hist(pval.mod1$pvalue,nclass=20, main=paste("Pvalues on ", varnom, sep=""))
```

#### qvalues ####
```{r}
pvalues<-pval.mod1$pvalue
qobj<-qvalue(p=pvalues)

summary(qobj)
# Controle of Fdr rate
#qobj<-qvalue(p=pvalues, fdr.level=0.1)
# name of marker having the controle fdr rate
#as.character(pval.mod1[which(qobj$significant==TRUE),1])
#  qvalues and pvalues
plot(qobj)
# hist of q values
hist(qobj)
# qq plot
qq(pval.mod1$pvalue,main=paste("QQplot on", varnom, sep=""))
```

####  Exploration of interested p-value ####

```{r}
# threshold à 4
seuil.mod1<--log10(0.0001)
assoc.mod1<-pval.mod1[-log10(pval.mod1$pvalue)>=seuil.mod1,]
dim(assoc.mod1)
# hist of positive value
hist(-log10(assoc.mod1[,2]), 
     main = paste("Distribution du -log de la pvalue pour ",varnom,sep=""),
     xlab= " logarithme de la pvalue" ,
     ylab= " Effectif" ,
     col="blue", border="white")
```

#### Position of marker according the study specie ####

```{r}
file_Phys_DRW<-"BREEDWHEAT_on_durum_physic.Rdata"

# Name of the file --> BLAST
load(file_Phys_DRW)

names(BLAST)
SNP_PHYS<-BLAST
SNP_PHYS<-as.data.frame(SNP_PHYS)
names(SNP_PHYS)

Positif_280K<-merge(SNP_PHYS,assoc.mod1, by.x=1, by.y=1 )
names(Positif_280K)

Positif_280K

```





