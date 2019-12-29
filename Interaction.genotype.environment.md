## Library 
#library(dplyr)
#library(FactoMineR)

## Description data, manage NA ##

```{r}
rm(list=ls())
data<-read.table("rapeseed_data_hybrides_techno.csv", sep=";", dec=".", h=T)
summary(data)
for (i in c(1:7)) { data[,i]<-as.factor(data[,i])}
data<-data%>%filter(!is.na(trait))%>%droplevels
data<-data%>%filter(!is.na(CU))%>%droplevels ## CU is the environment
hist(data$trait)
plot(trait~ YEAR, data)
year<-unique(data$YEAR  )

# Characterization of environment in each year 
list=c("2018")## Choice one year
data_2018<-data[data$YEAR%in%list,]
# Consider Environment as combination between site and year  
IQDFF<-as.data.frame(tapply(data$trait, list(data$MATERIAL, data$CU), mean,row.names=T))
IQDFF[1:2,1:2]

#Filter NA
x<-as.matrix(apply(IQDFF,1,function(x) sum(is.na(x))/length(x)))
colnames(x)<-"%NAM"
sokhna<-cbind(IQDFF ,x)
sokhna[1:3,1:3]
dax<-subset(sokhna, sokhna$`%NAM`<0.50)
y<-as.matrix(apply(dax,2,function(x) sum(is.na(x))/length(x)))
Y<-t(y)
rownames(Y)<-"%NAG"
dax_G<-rbind(dax, Y)
dax_F<-dax_G[-19,-c(2,17,18,19,22,26)]
Soh=dax_F

# Compute % NA after filter 
Nb<-as.data.frame(apply(Soh,1,function(x) sum(is.na(x))))
colnames(Nb)<-"NbNA"
dim(Soh)
dady<-(sum(Nb)/(nrow(Soh)*ncol(Soh)))*100###
write.table(Soh, file="traitIGxE.csv", sep=";")

# Imputation NA 
library(missForest)
data<-missForest(Soh)
data$OOBerror #criterion NRMSE is used to evaluate the precision of the imputation
test<-data$ximp ## extract imputed data
one<-t(test)
```

#### Typology of environments
```{r}
cs<-scale(one, center = TRUE, scale = TRUE) 
d1<- dist(cs, method = "euclidean") 
fit<- hclust(d1) # fait le dendrogramme 
plot(fit) # 
plot<-heatmap(as.matrix(test), Colv = NA, Rowv = NA, scale="row", xlab="environnement", ylab="mat?riel", main="heatmap")
```

#### ECOVALENCE
```{r}
IQFF<-as.matrix(Soh) ##Soh is the matrix having the mean of each genotype in each environment
nblieu<-ncol(IQFF)
nbgeno<-nrow(IQFF)
Ige<-matrix(0, nrow=nbgeno, ncol=nblieu)## Matrix name of interaction GxE

for(i in 1:nbgeno)
  for(j in 1:nblieu)
    Ige[i,j] <- IQFF[i,j] - mean(IQFF[,j]) - mean(IQFF[i,]) + mean(IQFF[,])
Ige[1:3,1:3]

Eco <- matrix(0, nrow=nbgeno, ncol=nblieu,dimnames=list(rownames(Ige), colnames(Ige)))# Stock the ecovalence value

for(i in 1:nbgeno)
  for(j in 1:nblieu)
    Eco[i,j] <- (Ige[i,j])^2

Ecovalence_Env <- rep(0, ncol(Eco))
names(Ecovalence_Env) <- colnames(Eco)
Ecovalence_geno <- rep(0, nrow(Eco))
names(Ecovalence_geno) <- rownames(Eco)


## compute ecovalence
for(i in 1:nbgeno)
  Ecovalence_geno[i] <- sum(Eco[i,])
for(j in 1:nblieu)
  Ecovalence_Env[j] <- sum(Eco[,j])

Eco_Env<-as.data.frame(Ecovalence_Env)
Eco_Env$Env<-colnames(Soh)
Eco_geno<-as.data.frame(Ecovalence_geno)
Eco_geno$geno<-paste("G", 1:nbgeno, sep="")#

## Compare the value of each ecovalence and their mean environment or genotype 

meanEnv<-as.data.frame(colMeans(soh))
colnames(meanEnv)<-"mean"
row.names(meanEnv)<-colnames(Soh)
meanEn$Env<-row.names(meanEnv)

meanGen<-as.data.frame(rowMeans(soh))
colnames(meanGen)<-"mean"
row.names(meanGen)<-row.names(Soh)
meanGen$geno<-row.names(meanGen)

##Ecovalence_environment 
P<-Eco_Env%>%
  mutate(Env = fct_reorder(Env, Ecovalence_Env)) %>%
  ggplot( aes(x=Env, y=Ecovalence_Env)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  ggtitle("Ecovalence Env")
theme_bw()

D<-meanEnv%>%
  mutate(Env = fct_reorder(Env, mean)) %>%
  ggplot( aes(x=Env, y=mean)) +
  geom_segment( aes(xend=Env, yend=0)) +
  geom_point( size=4, color="orange") +
  coord_flip() +
  theme_bw() +
  ggtitle("Mean Env")
xlab("")
library(gridExtra)
grid.arrange(P,D,ncol=2)

# Ecovalence_genotype

K<-Eco_geno%>%
  mutate(geno= fct_reorder(geno, Ecovalence_geno)) %>%
  ggplot( aes(x=geno, y=Ecovalence_geno)) +
  geom_bar(stat="identity", fill="#E69F00", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  ggtitle("Ecovalence Geno")
theme_bw()

Q<-meanGen%>%
  mutate(geno= fct_reorder(geno,mean)) %>%
  ggplot( aes(x=geno, y=mean)) +
  geom_segment( aes(xend=geno, yend=0)) +
  geom_point( size=2, color="orange") +
  coord_flip() +
  theme_bw() +
  ggtitle("Mean Geno")
xlab("")
library(gridExtra)
grid.arrange(K,Q,ncol=2)
```

## AMMI model

```{r}
colnames(Ige)<-colnames(data)## Ige is the matrix interaction GxE
row.names(Ige)<-row.names(data)
res.pca <- PCA(Ige, scale.unit=FALSE, graph = FALSE) 
res.pca$eig
barplot(res.pca$eig[,2], names=paste("Dim", 1:nrow(res.pca$eig)), las=1,
        ylab="% of the total variance", main="ACP of interaction matrix")

coord_var <- res.pca$var$coord
coord_var[,1]

coord_ind <- res.pca$ind$coord
coord_ind[,1]

moy_ind <- rowMeans(data)
moy_env <- colMeans(data)

plot(moy_ind, coord_ind[,1], pch=16, xlab=" Mean value of the genotypes",
     ylab="AMMI-CP1", ylim=c(-90,150))
text(moy_ind, coord_ind[,1], labels=names(moy_ind), cex=0.8, pos=3, offset=0.3)
abline(h=0, lty=2)
abline(v=mean(moy_ind), lty=2)

## 
plot(moy_env, coord_var[,1], pch=16, xlab="Valeur moyenne des env",
     ylab="AMMI-CP1", ylim=c(-70,70))
text(moy_env, coord_var[,1], labels=names(moy_env), cex=0.8, pos=3, offset=0.3)
abline(h=0, lty=2)
abline(v=mean(moy_env), lty=2)
```

## GGEBiplot
```{r}
library(GGEBiplots)
GGE<-GGEModel(Soh)### Soh is the matrix having the mean of each genotype in each environment
MeanStability(GGE)
WhichWon(GGE)
DiscRep(GGE)
RankEnv(GGE)
```


