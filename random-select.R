tableau=tableau[!is.na(tableau$trait ),]
tableau=tableau[,-c(11,12)]

for (i in c(1:9)) { tableau[,i]<-as.factor(tableau[,i])}
for (i in c(10)){ tableau[,i]<-as.numeric(tableau[,i])}

summary(tableau)
tableau$choice <- paste(tableau$CU , tableau$MATERIAL, sep='_')

Nb.line <- length(unique(paste(tableau$CU , tableau$MATERIAL, sep='_'))) *2
Nb.col <- ncol(tableau)
list.choice <- unique(paste(tableau$CU , tableau$MATERIAL, sep='_'))

new.tableau <- data.frame(matrix(nrow = Nb.line, ncol = Nb.col))
colnames(new.tableau) <- colnames(tableau)
new.tableau$choice <- rep(unique(tableau$choice), 2)


for(i in list.choice){
  
  mini.tab <- tableau[tableau$choice==i,]
  if(nrow(mini.tab)!=1){
    
    new.tableau[new.tableau$choice==i,][2,1] <-as.character(mini.tab[2,1])
    new.tableau[new.tableau$choice==i,][2,2] <-as.character(mini.tab[2,2])
    new.tableau[new.tableau$choice==i,][2,3] <-as.character(mini.tab[2,3])
    new.tableau[new.tableau$choice==i,][2,4] <-as.character(mini.tab[2,4])
    new.tableau[new.tableau$choice==i,][2,5] <-as.character(mini.tab[2,5])
    new.tableau[new.tableau$choice==i,][2,6] <-as.character(mini.tab[2,6])
    new.tableau[new.tableau$choice==i,][2,7] <-as.character(mini.tab[2,7])
    new.tableau[new.tableau$choice==i,][2,8] <-as.character(mini.tab[2,8])
    new.tableau[new.tableau$choice==i,][2,9] <-as.character(mini.tab[2,9])
    new.tableau[new.tableau$choice==i,][2,10] <-as.numeric(mini.tab[2,10])
    new.tableau[new.tableau$choice==i,][2,11] <-as.character(mini.tab[2,11])
    
    
    new.tableau[new.tableau$choice==i,][1,1] <-as.character(mini.tab[1,1])
    new.tableau[new.tableau$choice==i,][1,2] <-as.character(mini.tab[1,2])
    new.tableau[new.tableau$choice==i,][1,3] <-as.character(mini.tab[1,3])
    new.tableau[new.tableau$choice==i,][1,4] <-as.character(mini.tab[1,4])
    new.tableau[new.tableau$choice==i,][1,5] <-as.character(mini.tab[1,5])
    new.tableau[new.tableau$choice==i,][1,6] <-as.character(mini.tab[1,6])
    new.tableau[new.tableau$choice==i,][1,7] <-as.character(mini.tab[1,7])
    new.tableau[new.tableau$choice==i,][1,8] <-as.character(mini.tab[1,8])
    new.tableau[new.tableau$choice==i,][1,9] <-as.character(mini.tab[1,9])
    new.tableau[new.tableau$choice==i,][1,10] <-as.numeric(mini.tab[1,10])
    new.tableau[new.tableau$choice==i,][1,11] <-as.character(mini.tab[1,11])
    
  }
  
} 

new.tableau <- new.tableau[order(new.tableau$choice),]
new.tableau$REP <- rep(c(1,2), Nb.line/2)
new.tableau=new.tableau[!is.na(new.tableau$trait ),]
write.table(new.tableau, file="data_trait_equilibred.csv", sep=";")
