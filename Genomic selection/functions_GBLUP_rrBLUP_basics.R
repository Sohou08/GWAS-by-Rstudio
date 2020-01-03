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

nonsymMME<-function(Xex,W,yex,G,varg,vare){
  #yex, W, Xex expanded with 0's for individuals with no phenotype
  nfix=ncol(Xex)
  WG=W%*%G*varg
  LHS=crossprod(cbind(Xex,WG))/vare
  neq=ncol(LHS)
  # add covariance structure 
  LHS[(nfix+1):neq,(nfix+1):neq]=LHS[(nfix+1):neq,(nfix+1):neq]+G*varg
  RHS=t(cbind(Xex,WG))%*%yex/vare
  # LHS is not guaranteed to be full rank and we use
  # a generalized inverse from MASS (this is *not* an approximation; see Harville 1976 for details)
  LHSi=ginv(LHS)
  sol=LHSi%*%RHS
  # put back in the right scale, from aplha to u
  sol[(nfix+1):neq]=(varg*G)%*%sol[(nfix+1):neq]
  sol
}
