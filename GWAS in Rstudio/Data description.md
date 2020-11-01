
# Genetic determinism of one trait through GWAS approach (Section 1) #

```{r}
# Clean environment
rm(list=ls())
# Call library 
library(lattice)
library(QTLRel)
library(qqman)
library(igraph)
```
## Prepare phenotypic data ##
```{r}
GRAINS<-read.table(file, header=TRUE, sep=";", dec=".")
dim(GRAINS)
names(GRAINS)
head(GRAINS,1)
names(GRAINS)[1]<-"Code.unique.2010"
# distributions of all variables 
for ( i in 2:ncol(DATA)) { hist(DATA[,i], main=names(DATA)[i]) , freq= FALSE}

#hist in pairs
pairs(DATA[,c(2:ncol(DATA) )])

#Correlation traits, Outlier , replace outlier by NA and delete NA values ####
```{r}
cor(DATA[,-1], DATA[,"Zn"], use="pairwise.complete.obs")
cor(DATA[,-1], DATA[,"Fe"], use="pairwise.complete.obs")
DATA[which(DATA$Zn.Fe > 4 ),]
DATA$Fe[which(DATA$Zn.Fe > 4 )]<-NA
names(DATA)
for (i in 2:ncol(DATA_sub)) { 
   DATA_sub<-DATA_sub[which(!is.na(DATA_sub[,i])),]
  }
```

### Search synthetics variables ###

```{r}
m1<-lm(DATA[,"Cu"] ~ DATA[,"Fe"] )
A<-DATA[,"Fe"]^2
m2<-lm(DATA[,"Cu"] ~ DATA[,"Fe"] + A )
anova(m1,m2)
# Principal Components Analysis
# entering raw data and extracting PCs
# from the correlation matrix
fit <- princomp(DATA[,2:ncol(DATA_sub)], cor=TRUE)
summary(fit) # print variance accounted 
loadings(fit) # pc loadings
plot(fit,type="lines") 
fit$scores # the principal components
biplot(fit)

```
 
 ## Prepare Genotype data ##

```{r}
load("G_EPO.Rdata")
dim(G_EPO)
##Compute minor allelic frequency 
Freq_EPO<-freqall(G_EPO)
hist(Freq_EPO, main="Distribution of minor  allelic")

# Delete Minor allelic frequency for 5% threshold
MAF<-which(Freq_EPO<0.05|Freq_EPO>0.95)
G_EPO<-G_EPO[,-MAF]
dim(G_EPO)

# QTL Rel request AA, BB et AB format
GA<-G_EPO
GA<-as.matrix(GA)

GA[which(GA==0)]<-"AA"
GA[which(GA==1)]<-"AB"
GA[which(GA==2)]<-"BB"
```

## Merge phenotype and genotype data ##

```{r}
names(DATA)
head(rownames(G_EPO),1)
PG<-merge(DATA, GA, by.x=1, by.y=0)
dim(PG)
PG[1:3, 1:10]
```

