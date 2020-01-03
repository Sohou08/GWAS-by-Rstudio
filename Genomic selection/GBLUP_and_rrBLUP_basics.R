
#  ```{r required libraries}

library(MASS)
library(breedR)
library(car)

source("functions_GBLUP_rrBLUP_basics.R")

# ```


# ```{r read pedigree + phenotypes + genotypes}

pedigree_pm <- read.table(file = "pedigree_GS3_pm.txt", header = F)
sel_var_pedigree <- c("D", "dad", "mum")
names(pedigree_pm) <- sel_var_pedigree

record_pedigree <- nrow(pedigree_pm)
print(c("number records pedigree: ", record_pedigree), quote = F)

phenotype_pm <- read.table(file = "phenotype_GS3_pm.txt", header = F)
sel_var_phenotype <- c("D", "gener", "ebv_dbh", "ebv_ht", "ebv_sweep","r_dbh", "r_ht", "r_sweep")
names(phenotype_pm) <- sel_var_phenotype

record_phenotype <- nrow(phenotype_pm)
print(c("number records phenotype: ", record_phenotype), quote = F)

ped_phenotype_pm <- merge(pedigree_pm, phenotype_pm, by = "D")
record_ped_phenotype <- nrow(ped_phenotype_pm)
print(c("number records ped_phenotype: ", record_ped_phenotype), quote = F)

# read genotypes in fixed format
# note that this requires looking beforehand how many markers
# read fixed width formatted data

genotype_pm <- read.fwf("genotypes_GS3_pm_sorted.txt",skip=1,widths=c(33,3,rep(1,2257)))
record_genotype <- nrow(genotype_pm)
columns_genotype <- ncol(genotype_pm)
print(c("number records genotype: ", record_genotype), quote = F)

# ```


# ```{r setup}

# convert into matrix
# some systems cannot handle that many loci for direct solving
# so next line is commented to avoid taking all markers

M <- as.matrix(genotype_pm[,3:columns_genotype])

# in following line we take just sample markers out of total number
# try more if your system has enough memory
# sample_snp <- sample(3:2257,1000,replace = F)
# M=as.matrix(geno[,sample_snp])

# number of individuals, number of snps

num_ind <- nrow(M); num_snp <- ncol(M)

# setup storage matrix for individual BLUPs

GBV <- matrix(NA,nrow=num_ind,ncol=4)
colnames(GBV)=c("rrBLUP","gBLUP","pBLUP","pheno")

# ```


# ```{r Calculate allele frequencies}

# compute allele freq *in the founders* 
# founder individuals are all those with mum == 0 and dad == 0

founders <- which(ped_phenotype_pm$dad==0 & ped_phenotype_pm$mum==0)

freq=c()
for (i in 1:num_snp){
  # number of copies divided by number of 2*individuals
  freq[i]=mean(M[founders,i])/2
}

summary(freq)
hist(freq)
sumtwopq=sum(2*freq*(1-freq))

# ```

# ```{r reads phenotypes & adds noise}

# read phenotypes
# choose column from 5 to 10, where traits are in file

ped_phenotype_pm$trait <- NA
print(colnames(ped_phenotype_pm[5:10]))
ped_phenotype_pm$trait <- ped_phenotype_pm[,5]
# alternative ways to access the trait of interest
# ped_phenotype_pm$trait <- ped_phenotype_pm[c("ebv_dbh")]

hist(ped_phenotype_pm$trait)

# add random noise to y so that things are not so good !!!
ped_phenotype_pm$trait=ped_phenotype_pm$trait+rnorm(ped_phenotype_pm$trait)
hist(ped_phenotype_pm$trait)
GBV[,4] <- ped_phenotype_pm$trait

# ```

# ```{r sets variance components}

# because we use quite accurate BLUP's as pseudo-phenotypes, we will consider that h2=0.5, thus

vary=var(ped_phenotype_pm$trait)
h2=0.5
varg=h2*vary
vare=vary-varg
vara=varg/sumtwopq

# ```


# ```{r RRBLUP explicit}

# ------------------------------
# BLUPSNP, also known as RRBLUP
# ------------------------------

# VanRaden's BLUPSNP 
#center genotypic matrix with allele frequencies

Z=matrix(NA,num_ind,num_snp)
for (i in 1:num_snp){
  Z[,i]=M[,i]-2*freq[i]
}

# because we have "only" 2600 markers, we can construct RRBLUP explicitly

X=matrix(1,ncol=1,nrow=num_ind)

# setup MME equations

LHS_rr_XX <- t(X) %*% X
# LHS_rr_XX <- crossprod(X)
LHS_rr_XZ <- t(X) %*% Z
# LHS_rr_XZ <- crossprod(X, Z)
LHS_rr_ZX <- t(Z) %*% X
# LHS_rr_ZX <- crossprod(Z, X)
LHS_rr_ZZ <- t(Z) %*% Z + diag(1, ncol(Z))*(vare/vara)
# LHS_rr_ZZ <- crossprod(Z) + diag(1, ncol(Z))*(vare/vara)

LHS_rr_upper <- cbind(LHS_rr_XX,LHS_rr_XZ)
LHS_rr_lower <- cbind(LHS_rr_ZX,LHS_rr_ZZ)
LHS_rr_explicit <- rbind(LHS_rr_upper,LHS_rr_lower)

num_eq <- ncol(LHS_rr_explicit)

print(c("number of equations: ", ncol(LHS_rr_explicit)), quote = F)

RHS_rr_Xy <- t(X) %*% ped_phenotype_pm$trait
RHS_rr_Zy <- t(Z) %*% ped_phenotype_pm$trait

RHS_rr_explicit <- rbind(RHS_rr_Xy,RHS_rr_Zy)

sol_rr_explicit <- solve(LHS_rr_explicit, RHS_rr_explicit)

blup_snp_explicit <- sol_rr_explicit[2:num_eq]

plot(abs(blup_snp_explicit))

hist(blup_snp_explicit)

GBV[,1] <- Z %*% blup_snp_explicit

hist(GBV[,1], main = "distribution EBV")

plot(ped_phenotype_pm$trait, GBV[,1], main = "phenotype versus EBV")

print(c("correlation observation vs prediction: ", cor(ped_phenotype_pm$trait,GBV[,1])), quote = F)

# ```


# ```{r prepare matrix for GBLUP}

# ----------------------------------
# GBLUPs with relationship matrices
# -----------------------------------
# In the GBLUPs, we will NOT use the inverse of G to avoid the problem of G being not invertible
# in this case, we can use the traditional method equivalent to BLUP 
# we will use Harville's MME that allow for singular G. They are *the same* as the equations commonly
# used for RKHS models:
# (X R-1 X    X'R-1 WG      ) beta = (X'R-1 y)
# (G W'R-1 X  GW'R-1 WG+G  )  alpha= (GW'R-1 y)
# GEBV=G%*%alpha
# where G is the matrix of covariances (G*varg in our case)

# W contains 1's if the individual has phenotype and 0 otherwise
# all individuals in genotype have phenotype
W <- diag(1,num_ind,num_ind)

# make heatmap of G and A

# ```


# ```{r GBLUP VanRaden}

# ------------------
# VanRaden's GBLUP 1
# ------------------

GVR=tcrossprod(Z)/sumtwopq
sol=nonsymMME(Xex=X,W=W,yex=ped_phenotype_pm$trait,G=GVR,varg=varg,vare=vare)

GBV[,2]=sol[2:length(sol)]

hist(GBV[,2], main = "distribution EBV")

plot(ped_phenotype_pm$trait, GBV[,2], main = "phenotype versus EBV")

print(c("correlation observation vs prediction: ", cor(ped_phenotype_pm$trait, GBV[,2])), quote = F)

heatmap(GVR)

# ```


# ```{r BLUP pedigree}

res.pedBlup <- remlf90(fixed = trait ~ 1,
                       genetic = list(model = 'add_animal',
                                      pedigree = ped_phenotype_pm[c("D", "dad", "mum")],
                                      id = "D"),
                       data = ped_phenotype_pm,
                       method = "ai")
summary(res.pedBlup)

amat <- as.matrix(res.pedBlup$effects$genetic$effects$direct$structure.matrix)

heatmap(amat)

# ```

# ```{r ABLUP}

# ------------------
# pedigree BLUP
# ------------------

sol=nonsymMME(Xex=X,W=W,yex=ped_phenotype_pm$trait,G=amat,varg=varg,vare=vare)

GBV[,3]=sol[2:length(sol)]

hist(GBV[,3], main = "distribution EBV")

plot(ped_phenotype_pm$trait, GBV[,3], main = "phenotype versus EBV")

print(c("correlation observation vs prediction: ", cor(ped_phenotype_pm$trait, GBV[,3])), quote = F)

# ```

# ```{r scatterplot EBVs}

scatterplotMatrix(GBV[,c("rrBLUP","gBLUP","pBLUP","pheno")],
                  diagonal="histogram", smooth=FALSE,
                  main = "EBV comparisons")

# ```


