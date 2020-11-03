
# Glup Cross model
```{r}
rm(list=ls()); gc()
K<-read.csv("kinship_okTNS.csv", sep=";", h=T,row.names = 1) # kinship matrix
table<-read.csv("QDFF.reel.csv", sep=";", dec=".") # data pheno
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

pin(anss , h2 ~ (V1+V2) / ( V1+V2+V3) ) # 0.60

summary(anss)$varcomp
# h2 femelle
vt<-summary(anss2)$varcomp
F<-as.data.frame(anss$U$`u:FEMALE`)
M<-as.data.frame(anss$U$`u:MALE`)
# write.table(F, file="F.U.csv", sep=";")
# write.table(M, file="M.U.csv", sep=";")

# Cross validation
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
