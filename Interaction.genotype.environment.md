##Library 
library(dplyr)

##Description data 

rm(list=ls())
data<-read.table("rapeseed_data_hybrides_techno.csv", sep=";", dec=".", h=T)
summary(data)
for (i in c(1:7)) { data[,i]<-as.factor(data[,i])}
data<-data%>%filter(!is.na(trait))%>%droplevels
data<-data%>%filter(!is.na(CU))%>%droplevels ## CU is the environment
hist(data$trait)

plot(trait~ YEAR, data)
year<-unique(data$YEAR  )
####dans le cas de la caract?risation des lieux dans chaque ann?e
list=c("2018")###faire pour chaque ann?e
data_year<-data[data$YEAR%in%list,]
summary(data_year)

####Prendre en compte ? la fois lieux et ann?e 
###Matrice renfermant la moyenne de chaque g?notype dans chaque environnement
IQDFF<-as.data.frame(tapply(data$trait, list(data$MATERIAL, data$CU), mean,row.names=T))
IQDFF[1:2,1:2]
###Filtrer les donn?es manquantes( g?notype et environnement)
x<-as.matrix(apply(IQDFF,1,function(x) sum(is.na(x))/length(x)))
colnames(x)<-"%NAM"
sokhna<-cbind(IQDFF ,x)
sokhna[1:3,1:3]
dax<-subset(sokhna, sokhna$`%NAM`<0.50)
y<-as.matrix(apply(dax,2,function(x) sum(is.na(x))/length(x)))
Y<-t(y)
rownames(Y)<-"%NAG"
dax_G<-rbind(dax, Y)
dax_F<-dax_G[-19,-c(2,17,18,19,22,26)]#####suppirmer les lieux ayant plus de 50% de NA
Soh=dax_F
##Calculer le pourcentage de NA apres le filtre
Nb<-as.data.frame(apply(Soh,1,function(x) sum(is.na(x))))
colnames(Nb)<-"NbNA"
dim(Soh)
dady<-(sum(Nb)/(nrow(Soh)*ncol(Soh)))*100###
write.table(Soh, file="traitIGxE.csv", sep=";")

################Imputer les donn?es manquantes apres le filter
library(missForest)
data<-missForest(Soh)
data$OOBerror###crit?re NRMSE pour pr?cision l'imputation des donn?es
test<-data$ximp##extraire les donn?es imput?es
one<-t(test)
cs<-scale(one, center = TRUE, scale = TRUE) # cr?ation du fichier de donn?es centr?-r?duit

### Typologie des environnements

d1<- dist(cs, method = "euclidean") # Cr?ation de la matrice de distance
fit<- hclust(d1) # fait le dendrogramme 
plot(fit) # 
plot<-heatmap(as.matrix(test), Colv = NA, Rowv = NA, scale="row", xlab="environnement", ylab="mat?riel", main="heatmap")

###calcul ECOLVALENCE
###Calcul? d'abords la moyenne de la teneur en huile de chaque g?notype et environnement

##Calcul d'?covalence
##Reprendre la matrice de compl?te imput?e

IQFF<-as.matrix(Soh)
nblieu<-ncol(IQFF)
nbgeno<-nrow(IQFF)
Ige<-matrix(0, nrow=nbgeno, ncol=nblieu)##Matrice d'interaction GxE

for(i in 1:nbgeno)
  for(j in 1:nblieu)
    Ige[i,j] <- IQFF[i,j] - mean(IQFF[,j]) - mean(IQFF[i,]) + mean(IQFF[,])
Ige[1:3,1:3]

Eco <- matrix(0, nrow=nbgeno, ncol=nblieu,dimnames=list(rownames(Ige), colnames(Ige)))##matrix pour le calcul ?covalence pour stocker la veleur d'?covalance

for(i in 1:nbgeno)
  for(j in 1:nblieu)
    Eco[i,j] <- (Ige[i,j])^2

## initialisation des objets et noms
Ecovalence_Env <- rep(0, ncol(Eco))
names(Ecovalence_Env) <- colnames(Eco)
Ecovalence_geno <- rep(0, nrow(Eco))
names(Ecovalence_geno) <- rownames(Eco)


## calcul de l'?covalence
for(i in 1:nbgeno)
  Ecovalence_geno[i] <- sum(Eco[i,])
for(j in 1:nblieu)
  Ecovalence_Env[j] <- sum(Eco[,j])

Eco_Env<-as.data.frame(Ecovalence_Env)
Eco_Env$Env<-colnames(Soh)
Eco_geno<-as.data.frame(Ecovalence_geno)
Eco_geno$geno<-paste("G", 1:nbgeno, sep="")###cod? le nom des g?notypes de G1-GX

###Moyenne de la teneur en huile de chaque environnement
meanEnv<-as.data.frame(colMeans(soh))
colnames(meanEnv)<-"mean"
row.names(meanEnv)<-colnames(Soh)
meanEn$Env<-row.names(meanEnv)

###Moyenne de la teneur en huile de chaque g?notype
meanGen<-as.data.frame(rowMeans(soh))
colnames(meanGen)<-"mean"
row.names(meanGen)<-row.names(Soh)
meanGen$geno<-row.names(meanGen)


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

###Analyse Eco GEno


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

###########AMMI model
library(FactoMineR)
###Pour le model AMMI il faut reprendre la matrice d'interaction Ige en changeant le nom des colonnes et le noms des lignes
colnames(Ige)<-colnames(data)
row.names(Ige)<-row.names(data)
res.pca <- PCA(Ige, scale.unit=FALSE, graph = FALSE) # dans la librarie FactoMineR

## ----Valeurs propres---------------------------------------------------------
res.pca$eig
barplot(res.pca$eig[,2], names=paste("Dim", 1:nrow(res.pca$eig)), las=1,
        ylab="Pourcentage de la variance totale", main="ACP sur matrice des interactions")

## ----Coordonn?es des environnements -----------------------------------------------------------
coord_var <- res.pca$var$coord
coord_var[,1]

## ----Coordonn?es des individus-----------------------------------------------------------
coord_ind <- res.pca$ind$coord
coord_ind[,1]

## ------------------------------------------------------------------------
moy_ind <- rowMeans(data)
moy_annee <- colMeans(data)

plot(moy_ind, coord_ind[,1], pch=16, xlab="Valeur moyenne des genotypes",
     ylab="AMMI-CP1", ylim=c(-90,150))
text(moy_ind, coord_ind[,1], labels=names(moy_ind), cex=0.8, pos=3, offset=0.3)
abline(h=0, lty=2)
abline(v=mean(moy_ind), lty=2)




## Graphique des annees
plot(moy_annee, coord_var[,1], pch=16, xlab="Valeur moyenne des env",
     ylab="AMMI-CP1", ylim=c(-70,70))
text(moy_annee, coord_var[,1], labels=names(moy_annee), cex=0.8, pos=3, offset=0.3)
abline(h=0, lty=2)
abline(v=mean(moy_annee), lty=2)

####GGEBiplot

library(GGEBiplots)

GGE<-GGEModel(Soh)###Matrice filtrer et imputÃ©e
MeanStability(GGE)
WhichWon(GGE)
DiscRep(GGE)
RankEnv(GGE)


