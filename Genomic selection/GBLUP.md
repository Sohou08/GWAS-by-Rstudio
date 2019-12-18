
## GBLUP by mixed.model
```{r}
rm(list = ls());gc()

Y <- read.csv2('OIL_OLE_TNS_persee -comp.csv',dec = '.')# data pheno

K<-read.csv("kinship_okTNS.csv", sep=";", h=T,row.names = 1)# Kinship

Y<-Y%>%filter(!is.na(CU))
Y<-Y%>%filter(!is.na( trait))

Blup <- lmer(trait~ (1|CU) + (1|Material)+(1|REP), data =Y, REML=T)
rm(Y)

blup <- ranef(Blup)$Material
blup$geno<-row.names(blup)

rm(Blup); gc()## release space

Yblup<-blup[which(blup$geno %in% row.names(K)),]

rm(blup);gc()
summary(Yblup)

name.genotype <- Yblup$geno # 

Yblup <- Yblup[,1]

names(Yblup) <- name.genotype

Kblup<-K[row.names(K) %in% name.genotype,row.names(K) %in% name.genotype ]
Kblup[1:3,1:3]

K<-Kblup
Y<-Yblup
rm(Kblup,Yblup),gc()

K=as.matrix(K)
heatmap(K)

# Order geno and kinship
ordre <- numeric()
for(i in 1:length(Y))
  ordre <- c(ordre, which(name.genotype==row.names(K)[i]))
names(Y)
Y<- Y[ordre]

library(rrBLUP)

Gblup<-mixed.solve(y=Y, K=K, SE=T)

vu<-Gblup$Vu
Ve<-Gblup$Ve
H2 <- vu/(vu+Ve)

u<-Gblup$u
write.table(u, file="GEBV.trait.csv", sep=";")
```

## GBLUP by Synbreed

```{r}
#library(synbreedData)
library(synbreed)
X<-read.csv("Matrix_012_TNS.csv", sep=";", h=T,row.names = 1) ## geno matrix
X<-X[row.names(X) %in% names(Y), ]
X[1:3,1:3]
setequal(row.names(Y), row.names(X))

##subdivide the set  2: 3/4 train and 1/4 train
train.sample <- sample(1:109,81, replace = F)
test.sample <- c(1:109)[-train.sample]

train.K <- K[train.sample, train.sample]
test.K <- K[test.sample, test.sample]

train.Y <- Y[train.sample]
test.Y <- Y[test.sample]

Gblup_train<-mixed.solve(y=train.Y, K=train.K)

a_G_tran<-Gblup_train$u## additive value 
Gblup_test<-mixed.solve(y=test.Y, K=test.K )#
a_G_test<-Gblup_test$u

train.X <- X[train.sample,]

test.X <- X[test.sample, ]

fit_train <- mixed.solve(y=train.Y , Z=train.X)

beta.hat <- fit_train$u
a_X_tran<-as.matrix(train.X)%*%beta.hat

a.X.test <- as.matrix(test.X ) %*% beta.hat

tmp1 <- cor(a_G_tran, a_X_tran)

tmp2 <-cor(a_G_test,a.X.test)####resultat du cross validation
```

## Glup Cross model
```{r}
rm(list=ls()); gc()
K<-read.csv("kinship_okTNS.csv", sep=";", h=T,row.names = 1)# kinship matrix
table<-read.csv("QDFF.reel.csv", sep=";", dec=".")## data pheno
table$Material<-as.factor(table$Material)
table<-OLE
summary(table)
table$FEMALE
table$MALE

K1<-K[row.names(K)%in%table$FEMALE, row.names(K)%in%table$FEMALE];dim(K1)
K2<-K[row.names(K)%in%table$MALE,row.names(K)%in%table$MALE];dim(K2)

K1<-as.matrix(K1)
K2<-as.matrix(K2)
y<-table$QDFF
S<-kronecker(K1, K2)
library(sommer)
anss <- mmer(OLE~1,
             random=~vs(FEMALE,Gu=K1) + vs(MALE,Gu=K2),
             rcov=~units,
             data=table)

pin(anss , h2 ~ (V1+V2) / ( V1+V2+V3) )###0.60

summary(anss)$varcomp
##h2 femelle
vt<-summary(anss2)$varcomp
F<-as.data.frame(anss$U$`u:FEMALE`)
M<-as.data.frame(anss$U$`u:MALE`)
#write.table(F, file="F.U.csv", sep=";")
#write.table(M, file="M.U.csv", sep=";")

##Cross validation
y.trn<-table
vv1 <- which(!is.na(table$OLE))
vv2 <- sample(vv1, 28)
y.trn[vv2,"OLE"] <- NA

anss2 <- mmer(OLE~1,
              random=~vs(FEMALE,Gu=K1) + vs(MALE,Gu=K2),
              rcov=~units,
              data=y.trn)


zu1 <- model.matrix(~FEMALE-1,y.trn) %*% anss2$U$`u:FEMALE`$trait
zu2 <- model.matrix(~MALE-1,y.trn) %*% anss2$U$`u:MALE`$trait
u <- zu1+zu2+anss2$Beta[1,"Estimate"]
cor(u[vv2,], table$trait[vv2])

```
