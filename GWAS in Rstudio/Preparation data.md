
__Purpose: Genetic determinism of one trait through GWAS approach : Data description__


rm(list=ls())

```{r}
source("fonctions_EPO.r")##create at the begining
library(lattice)
library(QTLRel)
library(qqman)
library(igraph)
```

## Importation data


```{r}
# Importation data

GRAINS<-read.table(file, header=TRUE, sep=";", dec=".")

dim(GRAINS)
names(GRAINS)
head(GRAINS,1)

```

###### phenotypic data ######

```{r}
# Rename the first column which having the same name of the first column of genotypes files 
names(GRAINS)[1]<-"Code.unique.2010"
```
## Distribution trait

```{r}
# look the distributions of all variables 
for ( i in 2:ncol(DATA)) { hist(DATA[,i], main=names(DATA)[i]) , freq= FALSE}
#hist in pairs
pairs(DATA[,c(2:ncol(DATA) )])

```

## Correlation traits, Outlier , replace outlier by NA and delete NA values 
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

## Search synthetics variables

```{r}

m1<-lm(DATA[,"Cu"] ~ DATA[,"Fe"] )
A<-DATA[,"Fe"]^2
m2<-lm(DATA[,"Cu"] ~ DATA[,"Fe"] + A )
anova(m1,m2)
# Pricipal Components Analysis
# entering raw data and extracting PCs
# from the correlation matrix

fit <- princomp(DATA[,2:ncol(DATA_sub)], cor=TRUE)

summary(fit) # print variance accounted for
loadings(fit) # pc loadings
plot(fit,type="lines") # scree plot
fit$scores # the principal components
biplot(fit)

```
###### Genotype data ######

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

# QTL Rel demande des formats AA, BB et AB
GA<-G_EPO
GA<-as.matrix(GA)

GA[which(GA==0)]<-"AA"
GA[which(GA==1)]<-"AB"
GA[which(GA==2)]<-"BB"
```

###### Merger data phenotype and genotype ######

```{r}
names(DATA)
head(rownames(G_EPO),1)
PG<-merge(DATA, GA, by.x=1, by.y=0)
dim(PG)
PG[1:3, 1:10]
```

# GWAS & le calcul des pvalues

## Choix de la variable à analyser et préparation du fichier.

```{r}
head(names(PG), 20)

```

On voit que le premier marqueur commence en 18.
Revenir à cette étape pour changer de variable. Changer la valeur de i.


```{r}
# mettre ici le numero de la variable
i<-4

Y<-PG[,i]

liste<- which(!(is.na(Y)))
Y<-Y[liste]
  
hist(Y, main=names(PG)[i])

# les marqueurs commencent en 18
G<-PG[liste, 18:dim(PG)[2]]

length(Y)


```

## Calcul de la matrice de Kinship

```{r}

# avec 10000 marqueurs tirés au hasard ca passe en 2-3 minutes 
# la foncton genMatrix de QTLRel fait le calcul des matrices 
# ici il n'y en a que 1000
liste <- sample(colnames(G),1000,replace=FALSE)
K<-genMatrix(G[,liste])
I<-diag(length(Y))

hist(K$AA)

```


## Calcul des composantes de la variance et estimation des covariances dues à la kinship

Le modèle estime une variance génétique additive et une variance environnementale.

```{r}

# le premier modele sert a calculer la variance
mod1<-estVC(y=Y,v=list(AA=K$AA,DD=NULL,HH=NULL,AD=NULL,MH=NULL,EE=I))
mod1


```

## Boucle sur tous les marqueurs retenus

La fonction scanOne va faire une analyse par marqueur en utilisant les composantes de la variance calculées précédemment.

```{r}

print(date())
GWAS.mod1<-scanOne(y=Y,gdat=G,vc=mod1,test="Chisq")
print(date())

# il faut 25 min pour les 93000 marqueurs

```

## Stockage des pvalues pour une étude ultérieure
```{r}

pval.mod1<-data.frame(marqueurs=colnames(G),pvalue=GWAS.mod1$p)

# stocker les valeurs 
varnom<-names(PG)[i]

file<-paste("GRAINS_GWAS.mod1_",varnom,".Rdata", sep="") 

save(GWAS.mod1, file=file, compress= TRUE)

file<-paste("GRAINS.pval.mod1_",varnom,".Rdata", sep="")
save(pval.mod1, file=file, compress= TRUE)

```


# Post traitement des données

Cette étape peut repartir des fichiers de pvalues stockées.

```{r}
library("data.table")
library(ggplot2)

```


## Charger les packages pour l'étude des qvalues 

Pas indispensables ici mais utile de réfléchir au nombre de faux positifs
```{r}
#http://www.bioconductor.org/packages/release/bioc/vignettes/qvalue/inst/doc/qvalue.pdf
#https://github.com/jdstorey/qvalue

#install.packages("devtools")

library("devtools")

#install_github("jdstorey/qvalue")
library(qvalue)

```


## Données de pvalues
Il s'agit là de récuperer pour l instant simplement les noms des variables étudiées


```{r}
varnom <- "moy_prot"

file<-paste("GRAINS_GWAS.mod1_",varnom,".Rdata", sep="")

load(file)

pvalues<-GWAS.mod1$p
length(GWAS.mod1$p)

pval.mod1<-data.frame(Marker_ID=names(GWAS.mod1$p),pvalue=pvalues)


```

## La valeur de la pvalue du SNP le plus significatif
```{r}

min(pval.mod1$pvalue)


```

## La distribution des pvalues
```{r}
hist(pval.mod1$pvalue,nclass=20, main=paste("Pvalues on ", varnom, sep=""))



```

## les qvalues
```{r}
pvalues<-pval.mod1$pvalue

qobj<-qvalue(p=pvalues)

summary(qobj)

# controle du taux de fdr
#qobj<-qvalue(p=pvalues, fdr.level=0.1)
# nom des marqueurs a taux de fdr controle
#as.character(pval.mod1[which(qobj$significant==TRUE),1])

#  qvalues and pvalues
plot(qobj)

# histogrammes des q values
hist(qobj)


# qq plot
qq(pval.mod1$pvalue,main=paste("QQplot on", varnom, sep=""))
```

Conclusions 
En ce qui concerne le caractère moy_prot il ne semble pas y avoir de SNP lié au caractère...


## Exploration des valeurs intéressantes
```{r}

# seuil à 4
seuil.mod1<--log10(0.0001)

assoc.mod1<-pval.mod1[-log10(pval.mod1$pvalue)>=seuil.mod1,]
dim(assoc.mod1)

# histrogramme des pvalues positives
hist(-log10(assoc.mod1[,2]), 
     main = paste("Distribution du -log de la pvalue pour ",varnom,sep=""),
     xlab= " logarithme de la pvalue" ,
     ylab= " Effectif" ,
     col="blue", border="white")


```


## Les données de positions des marqueurs sur le blé dur

```{r}
file_Phys_DRW<-"BREEDWHEAT_on_durum_physic.Rdata"

# le fichier s'appelle BLAST
load(file_Phys_DRW)

names(BLAST)
SNP_PHYS<-BLAST
SNP_PHYS<-as.data.frame(SNP_PHYS)
names(SNP_PHYS)

```


## Les positions des marqueurs intéressants sont elle connues ?

```{r}

Positif_280K<-merge(SNP_PHYS,assoc.mod1, by.x=1, by.y=1 )
names(Positif_280K)

Positif_280K
```

Chr le chromosome ou se trouve le SNP candidat
Start sa position physique sur le génome de référence
ensuite suivent les paramètres de BLAST
et enfin la pvalue est redonnée




