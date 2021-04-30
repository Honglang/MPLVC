###########################################
## Library loading and functions  #########
###########################################
library("MASS")
library("mvtnorm")
library("Matrix")

library(tidyverse)
library(stringr)
library(huge)

myspline = function(xx, K, d){
  Jn = d+K+1
  
  S = matrix(0, nrow=length(xx), ncol=Jn)	 
  if (K == 0){ 
    for(i in 1:(d+1)){
      S[,i]=xx^(i-1)
    }
  }  else {
    for(i in 1:(d+1)){
      S[,i]=xx^(i-1)
    }
    
    knots=quantile(xx, probs = (1:K)/(K+1))
    
    for(j in 1:K){
      
      S[,d+1+j]=((xx-knots[j])^d)*(xx>=knots[j]) 
    } 
  }
  S
}

golden = function(f, lower, upper) { 
  
  golden.ratio = 2/(sqrt(5) + 1)
  x1 = upper - golden.ratio*(upper - lower)
  x2 = lower + golden.ratio*(upper - lower)
  
  f1 = f(x1)
  f2 = f(x2)
  
  iteration = 0
  
  while (abs(upper - lower) > 0.05)
  {
    iteration = iteration + 1
    
    if (f2 > f1)
    {
      ### Set the new upper bound
      upper = x2
      ### Set the new upper test point
      ### Use the special result of the golden ratio
      x2 = x1
      f2 = f1
      
      ### Set the new lower test point
      x1 = upper - golden.ratio*(upper - lower)
      f1 = f(x1)
    } 
    else 
    {
      # the minimum is to the right of x1
      # let x1 be the new lower bound
      # let x2 be the new lower test point
      
      ### Set the new lower bound
      lower = x1
      
      ### Set the new lower test point
      x1 = x2
      
      f1 = f2
      
      ### Set the new upper test point
      x2 = lower + golden.ratio*(upper - lower)
      f2 = f(x2)
    }
  }
  
  ### Use the mid-point of the final interval as the estimate of the optimzer
  estimated.minimizer = (lower + upper)/2
  estimated.minimizer
}



##################################################
### Functions for Model Selection using BIC ######
##################################################

BIC_marginal = function(t, Z, y, G,mvec, mvec0,mvec1, mlist,K0, d0, K1, d1){
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  D1 = c(rep(0,d0+1), rep(1,K0))
  D2 = c(rep(0,d1+1), rep(1,K1))
  D = diag(c(D1, rep(0,dz), D2))
  Nobs=length(t)
  
  ##step 0: initial values
  
  B0 = myspline(t, K0, d0)
  B1 = myspline(t, K1, d1)
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B_temp = cbind(B0, Z, B1*G_ext)
  
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value
  
  run = 0
  while(run <= 50){
    
    run = run+1
    
    ngamma = (K0+d0+1+K1+d1+1+dz)
    arsumg = matrix(rep(0,J*ngamma),nrow=J*ngamma)
    arsumc = matrix(rep(0,J*ngamma*J*ngamma),nrow=J*ngamma)
    arsumgfirstdev1 = matrix(rep(0,J*ngamma*ngamma),nrow=J*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi = Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)
      B1i = myspline(ti, K1, d1)
      
      G_exti = G_ext[seq,]
      
      Bi = cbind(B0i, Zi, B1i*G_exti)
      
      mui = Bi%*%gamma.old
      mudoti = Bi
      fmui = mui  #link function
      fmui_dev = diag(ni)  # first dev of link function
      vmui = diag(ni) #variance of mui
      
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% fmui_dev %*% vmui %*% mlist[[j]][[i]] %*% vmui
        gi=c(gi,(1/N)*wi %*% (yi-mui))
        firstdev=rbind(firstdev,- (1/N)* wi %*% fmui_dev %*% mudoti)#first dev of gij
      }
      
      arsumc = arsumc + gi %*% t(gi)
      arsumg = arsumg + gi
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc)
    
    Q = t(arsumg) %*% arcinv1 %*% arsumg
    
    arqif1dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumg)/N
    arqif2dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1)/N
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  
  # based on the estimation, calculate Q, and BIC
  
  QN = Q
  r=length(arsumg)
  k=length(gamma.old)
  BIC = QN+(r-k)*log(N)
  
  
  return(BIC)
  
}


BIC_marginal_null = function(t, Z, y,mvec, mvec0,mvec1, mlist, K0, d0){
  
  D1 = c(rep(0,d0+1), rep(1,K0))
  dz=ncol(Z)
  D = diag(c(D1, rep(0,dz)))
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  ##step 0: initial values
  
  B0 = myspline(t, K0, d0)
  B_temp = cbind(B0, Z)
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value
  
  run = 0
  while(run <= 50){
    
    run = run+1
    
    ngamma = (K0+d0+1+dz)
    arsumg = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumc = matrix(rep(0,2*ngamma*2*ngamma),nrow=2*ngamma)
    #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumgfirstdev1 = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi = Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)
      
      Bi = cbind(B0i, Zi)
      
      mui = Bi%*%gamma.old
      mudoti = Bi
      fmui = mui  #link function
      fmui_dev = diag(ni)  # first dev of link function
      vmui = diag(ni) #variance of mui
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% fmui_dev %*% vmui %*% mlist[[j]][[i]] %*% vmui
        gi=c(gi,(1/N)*wi %*% (yi-mui))
        firstdev=rbind(firstdev,- (1/N)* wi %*% fmui_dev %*% mudoti)#first dev of gij
      }
      
      
      arsumc = arsumc + gi %*% t(gi)
      arsumg = arsumg + gi
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc)
    
    Q = t(arsumg) %*% arcinv1 %*% arsumg
    
    arqif1dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumg)/N
    arqif2dev1 = (2*t(arsumgfirstdev1) %*% arcinv1 %*% arsumgfirstdev1)/N
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  
  
  QN = Q
  
  k=length(gamma.old)
  BIC = QN+(J-1)*k*log(N)
  
  return(BIC)
  
}



#############################################################
#### Estimation and Inference under full joint model  #######
#############################################################

##################################################
### Estimation ###################################
##################################################
## Qstat: with input the whole data (t, Z, G, ylist, mvec, mvec0, mvec1), model specification (Klist, dlist),
## as well as the working correlation structure (mlist) and the parameter for estimation.
## the output of the function: Qstat itself without penalization, the first derivative of Qstat
## (which will be used for the gradient decent method for the estimation of the parameters) and
## the second derivative of Qstat (which will be used for the GCV function for the selection of smoothing tuning parameter)

## Here the order of the parameter gamma is intercept, partial linear parameter, and the slope
## And in the GEE we assume the marginal variances are the same so we omit the A_i part since that will
## be canceled in the final formulation of Qstat
Qstat=function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma){
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  ngamma=length(gamma)
  arsumg = matrix(rep(0,J*ngamma),nrow=J*ngamma)
  arsumc = matrix(rep(0,J*ngamma*J*ngamma),nrow=J*ngamma)
  arsumgfirstdev1 = matrix(rep(0,J*ngamma*ngamma),nrow=J*ngamma)
  
  
  for(i in 1:N){
    
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    Zi = Z[seq,,drop=FALSE]
    ti = t[seq]
    Gi=G[seq]
    yi=NULL
    B0il=vector("list",L)
    Zil=vector("list",L)
    B1Gil=vector("list",L)
    for(l in 1:L){
      yi=c(yi,ylist[[l]][seq])
      K0=Klist[[l]][1]
      d0=dlist[[l]][1]  
      K1=Klist[[l]][2]
      d1=dlist[[l]][2]
      B0il[[l]]=myspline(ti,K0,d0)
      G_exti = matrix(rep(Gi,K1+d1+1), nrow=ni, ncol=K1+d1+1)
      B1Gil[[l]]=myspline(ti,K1,d1)*G_exti
      Zil[[l]]=Zi
    }
    
    
    B0i_2 = as.matrix(bdiag(B0il))
    Zi_2 = as.matrix(bdiag(Zil))
    B1i_2 = as.matrix(bdiag(B1Gil))
    
    Bi_2 = cbind(B0i_2, Zi_2, B1i_2)
    
    mui = Bi_2%*%gamma
    mudoti = Bi_2
    
    gi=NULL
    firstdev=NULL
    for(j in 1:J){
      wi = t(mudoti) %*% mlist[[j]][[i]] 
      gi=c(gi,wi %*% (yi-mui))
      firstdev=rbind(firstdev,-wi %*% mudoti)#first dev of gij
    }
    
    arsumc = arsumc + gi %*% t(gi) 
    arsumg = arsumg + gi 
    arsumgfirstdev1 = arsumgfirstdev1 + firstdev
  }
  
  arcinv1 = ginv(arsumc/N) 
  
  Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N)
  
  Qp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N)
  Qpp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  list(Q=Q,Qp=Qp,Qpp=Qpp,arcinv1=arcinv1,arsumgfirstdev1=arsumgfirstdev1)
}


estimation_full = function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,lambda){
  #t: a long vector which stacks all of the measurement time points for all individuals
  #Z: a long (tall) matrix which stacks all of the other covariates measured at t
  #ylist: a list of multiple responses, each response vector is a long vector which stacks the measurement at t
  #mvec: a vector of n_i
  #mvec0: a vector of start measuring points of individuals
  #mvec1: a vector of end measuring points of individuals
  #G: a long vector of genetic information
  #Klist: list of number of knots for coefficient functions
  #dlist: list of degrees for coefficient functions
  #mlist: list of basis matrices
  #lambda: smoothing tunning parameter
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  D0=NULL
  D1=NULL
  gamma.old0=NULL
  gamma.old1=NULL
  alpha.old=NULL
  for(l in 1:L){
    K0=Klist[[l]][1]
    d0=dlist[[l]][1]  
    K1=Klist[[l]][2]
    d1=dlist[[l]][2]  
    temp0=c(rep(0,d0+1),rep(1,K0))
    temp1=c(rep(0,d1+1),rep(1,K1))
    D0=c(D0,temp0)
    D1=c(D1,temp1)
    B0=myspline(t,K0,d0)
    G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
    B1G=myspline(t,K1,d1)*G_ext
    B_temp = cbind(B0, Z, B1G)
    gamma.l=ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%ylist[[l]]
    gamma.old0=c(gamma.old0,gamma.l[1:(K0+d0+1)])
    alpha.old=c(alpha.old,gamma.l[K0+d0+1+(1:dz)])
    gamma.old1=c(gamma.old1,gamma.l[(K0+d0+2+dz):length(gamma.l)])
  }
  
  D = diag(c(D0, rep(0,L*dz), D1))    
  gamma.old=c(gamma.old0,alpha.old,gamma.old1)##simple LS, initial value 
  ngamma=length(gamma.old)
  
  run = 0  
  while(run <= 50){
    run = run+1
    res=Qstat(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma.old)
    arqif2dev1=res$Qpp/N+ 2*lambda*D	
    invarqif2dev1 = ginv(arqif2dev1)
    arqif1dev1=res$Qp/N+2*lambda*D%*%gamma.old
    gamma.new = gamma.old -  invarqif2dev1%*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    gamma.old = gamma.new
    #print(gamma.old)
  } #loop for gamma
  res=Qstat(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma.old) 
  QNdotdot = res$Qpp
  QN = res$Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  
  list(GCV=GCV, gamma=gamma.old, QN=QN)
  
} 


estimation_cov = function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,lambda,gamma){
  #t: a long vector which stacks all of the measurement time points for all individuals
  #Z: a long (tall) matrix which stacks all of the other covariates measured at t
  #ylist: a list of multiple responses, each response vector is a long vector which stacks the measurement at t
  #mvec: a vector of n_i
  #mvec0: a vector of start measuring points of individuals
  #mvec1: a vector of end measuring points of individuals
  #G: a long vector of genetic information
  #Klist: list of number of knots for coefficient functions
  #dlist: list of degrees for coefficient functions
  #mlist: list of basis matrices
  #lambda: smoothing tunning parameter
  #alpha_ture: true coefficients
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  D0=NULL
  D1=NULL
  for(l in 1:L){
    K0=Klist[[l]][1]
    d0=dlist[[l]][1]  
    K1=Klist[[l]][2]
    d1=dlist[[l]][2]  
    temp0=c(rep(0,d0+1),rep(1,K0))
    temp1=c(rep(0,d1+1),rep(1,K1))
    D0=c(D0,temp0)
    D1=c(D1,temp1)
  }
  
  D = diag(c(D0, rep(0,L*dz), D1))    
  ngamma=length(gamma)
  
  
  res=Qstat(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma) 
  
  arqif2dev1=res$Qpp/N+ 2*lambda*D	
  arcinv1=res$arcinv1
  arsumgfirstdev1=res$arsumgfirstdev1
  lam_gam = lambda*D%*%gamma
  arsums = matrix(rep(0,ngamma*ngamma),nrow=ngamma)
  
  for (i in 1:N){
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    Zi = Z[seq,,drop=FALSE]
    ti = t[seq]
    Gi=G[seq]
    yi=NULL
    B0il=vector("list",L)
    Zil=vector("list",L)
    B1Gil=vector("list",L)
    for(l in 1:L){
      yi=c(yi,ylist[[l]][seq])
      K0=Klist[[l]][1]
      d0=dlist[[l]][1]  
      K1=Klist[[l]][2]
      d1=dlist[[l]][2]
      B0il[[l]]=myspline(ti,K0,d0)
      G_exti = matrix(rep(Gi,K1+d1+1), nrow=ni, ncol=K1+d1+1)
      B1Gil[[l]]=myspline(ti,K1,d1)*G_exti
      Zil[[l]]=Zi
    }
    
    
    B0i_2 = as.matrix(bdiag(B0il))
    Zi_2 = as.matrix(bdiag(Zil))
    B1i_2 = as.matrix(bdiag(B1Gil))
    
    Bi_2 = cbind(B0i_2, Zi_2, B1i_2)
    
    mui = Bi_2%*%gamma
    mudoti = Bi_2
    
    gi=NULL
    for(j in 1:J){
      wi = t(mudoti) %*% mlist[[j]][[i]]
      gi=c(gi,wi %*% (yi-mui))
    }
    si = t(arsumgfirstdev1/N) %*% arcinv1 %*% gi + lam_gam #changed!!!!!!
    
    arsums = arsums + si%*%t(si)
    
  }
  
  
  cov_part = ginv(arqif2dev1*N/2) # this is the original code and it's correct
  cov_gamma = cov_part%*%arsums%*%cov_part # sandwich formula
  # cov_gamma = ginv(t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N))/N #=ginv(res$Qpp/2)
  # this is the direct one not using sandwich formula
  
  return(cov_gamma)
  
} 



############################################################
### 2 Version of Joint Effect test            ##############
### And interaction test ###################################
############################################################

### Testing

########################################################
########## Joint Testing  ##############################
########################################################

Qstat.null=function(t, Z, ylist,mvec, mvec0,mvec1, Klist,dlist, mlist,gamma){
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  ngamma=length(gamma)
  arsumg = matrix(rep(0,J*ngamma),nrow=J*ngamma)
  arsumc = matrix(rep(0,J*ngamma*J*ngamma),nrow=J*ngamma)
  arsumgfirstdev1 = matrix(rep(0,J*ngamma*ngamma),nrow=J*ngamma)
  
  
  for(i in 1:N){
    
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    Zi = Z[seq,,drop=FALSE]
    ti = t[seq]
    yi=NULL
    B0il=vector("list",L)
    Zil=vector("list",L)
    
    for(l in 1:L){
      yi=c(yi,ylist[[l]][seq])
      K0=Klist[[l]][1]
      d0=dlist[[l]][1]  
      
      B0il[[l]]=myspline(ti,K0,d0)
      
      Zil[[l]]=Zi
    }
    
    
    B0i_2 = as.matrix(bdiag(B0il))
    Zi_2 = as.matrix(bdiag(Zil))
    
    
    Bi_2 = cbind(B0i_2, Zi_2)
    
    mui = Bi_2%*%gamma
    mudoti = Bi_2
    
    
    gi=NULL
    firstdev=NULL
    for(j in 1:J){
      wi = t(mudoti) %*% mlist[[j]][[i]] 
      gi=c(gi,wi %*% (yi-mui))
      firstdev=rbind(firstdev,-wi %*% mudoti)#first dev of gij
    }
    
    arsumc = arsumc + gi %*% t(gi) 
    arsumg = arsumg + gi 
    arsumgfirstdev1 = arsumgfirstdev1 + firstdev
  }
  
  arcinv1 = ginv(arsumc/N) 
  
  Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N)
  
  Qp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N)
  Qpp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  
  
  list(Q=Q,Qp=Qp,Qpp=Qpp,arcinv1=arcinv1,arsumgfirstdev1=arsumgfirstdev1)
}



estimation_null = function(t, Z, ylist,mvec, mvec0,mvec1, Klist,dlist, mlist,lambda){
  #t: a long vector which stacks all of the measurement time points for all individuals
  #Z: a long (tall) matrix which stacks all of the other covariates measured at t
  #ylist: a list of multiple responses, each response vector is a long vector which stacks the measurement at t
  #mvec: a vector of n_i
  #mvec0: a vector of start measuring points of individuals
  #mvec1: a vector of end measuring points of individuals
  #G: a long vector of genetic information
  #Klist: list of number of knots for coefficient functions
  #dlist: list of degrees for coefficient functions
  #mlist: list of basis matrices
  #lambda: smoothing tunning parameter
  
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  D0=NULL
  D1=NULL
  gamma.old0=NULL
  gamma.old1=NULL
  alpha.old=NULL
  for(l in 1:L){
    K0=Klist[[l]][1]
    d0=dlist[[l]][1]  
    temp0=c(rep(0,d0+1),rep(1,K0))
    D0=c(D0,temp0)
    B0=myspline(t,K0,d0)
    B_temp = cbind(B0, Z)
    gamma.l=ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%ylist[[l]]
    gamma.old0=c(gamma.old0,gamma.l[1:(K0+d0+1)])
    alpha.old=c(alpha.old,gamma.l[K0+d0+1+(1:dz)])
  }
  
  D = diag(c(D0, rep(0,L*dz)))    
  gamma.old=c(gamma.old0,alpha.old)##simple LS, initial value 
  ngamma=length(gamma.old)
  
  run = 0  
  while(run <= 50){
    run = run+1
    res=Qstat.null(t, Z, ylist,mvec, mvec0,mvec1,Klist,dlist, mlist,gamma.old)
    arqif2dev1=res$Qpp/N+ 2*lambda*D	
    invarqif2dev1 = ginv(arqif2dev1)
    arqif1dev1=res$Qp/N+2*lambda*D%*%gamma.old
    gamma.new = gamma.old -  invarqif2dev1%*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    gamma.old = gamma.new
    #print(gamma.old)
  } #loop for gamma
  res=Qstat.null(t, Z, ylist,mvec, mvec0,mvec1, Klist,dlist, mlist,gamma.old) 
  QNdotdot = res$Qpp
  QN = res$Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  
  list(GCV=GCV, gamma=gamma.old, QN=QN)
  
} 


joint_test <- function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist){
  result1 = estimation_full(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,0)
  #gamma1 = result1$gamma
  Q1 = result1$QN#Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma1)
  
  Klist0=list(Klist[[1]][1],Klist[[2]][1])
  dlist0=list(dlist[[1]][1],dlist[[2]][1])
  result0 = estimation_null(t, Z, ylist,mvec, mvec0,mvec1, Klist0,dlist0, mlist, 0)
  Q0=result0$QN
  # gamma_temp = result0[[2]]
  # gamma0 = c(gamma_temp, rep(0,2*(K1+d1+1)))
  # Q0 = Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma0)
  df_test =Klist[[1]][2]+Klist[[2]][2]+dlist[[1]][2]+dlist[[2]][2]+2 #2*(K1+d1+1)
  
  gamma0=c(result0$gamma,rep(0,df_test))
  result00=Qstat(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma0)
  #print(result00)
  Q0=result00$Q
  
  test = Q0 - Q1
  
  
  #cutoff = qchisq(0.95, df_test)
  
  #rej = (test > cutoff)*1
  pvalue=pchisq(test,df_test,lower.tail = FALSE)
  return(list(test=test,pvalue=pvalue,df_test=df_test))
}


########################################################
########## Marginal Testing  ###########################
########################################################
Qstat_marginal_full = function(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1, mlist,lambda){
  dz=ncol(Z)
  D1 = c(rep(0,d0+1), rep(1,K0))
  D2 = c(rep(0,d1+1), rep(1,K1))
  D = diag(c(D1, rep(0,dz), D2))    
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  
  ##step 0: initial values
  
  B0 = myspline(t, K0, d0)
  B1 = myspline(t, K1, d1)               
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B_temp = cbind(B0, Z, B1*G_ext)
  
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value 
  
  run = 0  
  while(run <= 50){
    
    run = run+1
    
    ngamma = (K0+d0+1+K1+d1+1+dz)
    arsumg = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumc = matrix(rep(0,2*ngamma*2*ngamma),nrow=2*ngamma)
    #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumgfirstdev1 = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi =Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)
      B1i = myspline(ti, K1, d1) 
      
      G_exti = G_ext[seq,]   		 
      
      Bi = cbind(B0i, Zi, B1i*G_exti)
      
      mui = Bi%*%gamma.old
      mudoti = Bi
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% mlist[[j]][[i]]
        gi=c(gi,wi %*% (yi-mui))
        firstdev=rbind(firstdev,- wi %*% mudoti)#first dev of gij
      }
      
      arsumc = arsumc + gi %*% t(gi) 
      arsumg = arsumg + gi 
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc/N)
    
    Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N) 
    
    arqif1dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N))/N + 2*lambda*D%*%gamma.old
    arqif2dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1)/N)/N + 2*lambda*D	
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  QN = Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  k=length(gamma.old)
  BIC = QN+(J-1)*k*log(N)
  
  list(GCV=GCV, QN=QN,gamma=gamma.old,BIC=BIC)
  
} 




Qstat_marginal_null = function(t, Z, y,mvec, mvec0,mvec1, K0, d0, mlist,lambda){
  
  D1 = c(rep(0,d0+1), rep(1,K0))
  dz=ncol(Z)
  D = diag(c(D1, rep(0,dz)))
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  
  ##step 0: initial values
  
  B0 = myspline(t, K0, d0)              
  B_temp = cbind(B0, Z)
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value 
  
  run = 0  
  while(run <= 50){
    
    run = run+1
    
    ngamma = (K0+d0+1+dz)
    arsumg = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumc = matrix(rep(0,2*ngamma*2*ngamma),nrow=2*ngamma)
    #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumgfirstdev1 = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi = Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)                                       
      
      Bi = cbind(B0i, Zi)
      
      mui = Bi%*%gamma.old
      mudoti = Bi
      
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% mlist[[j]][[i]]
        gi=c(gi,wi %*% (yi-mui))
        firstdev=rbind(firstdev,-wi %*% mudoti)#first dev of gij
      }
      
      arsumc = arsumc + gi %*% t(gi) 
      arsumg = arsumg + gi 
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc/N)
    
    Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N) 
    
    arqif1dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg)/N)/N + 2*lambda*D%*%gamma.old
    arqif2dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1)/N)/N + 2*lambda*D 	
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  cat("run:")
  print(run)
  cat("gammadiff:")
  print(gammadiff)
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  QN = Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  k=length(gamma.old)
  BIC = QN+(J-1)*k*log(N)
  
  list(GCV=GCV, QN=QN, gamma=gamma.old,BIC=BIC)
  
} 

# res1=Qstat_marginal_null(t, Z, ylist[[1]], mvec, mvec0,mvec1,K01, d01,mlist1, 0)
# res2=Qstat.marginal(t, Z, ylist[[1]], mvec, mvec0,mvec1, G, K01, d01, K11, d11, mlist1, c(res1$gamma,rep(0,K11+d11+1)))
# res1=Qstat_marginal_full(t, Z, ylist[[1]], mvec, mvec0,mvec1,G, K01, d01, K11, d11,mlist1, 0)
# res2=Qstat.marginal(t, Z, ylist[[1]], mvec, mvec0,mvec1, G, K01, d01, K11, d11, mlist1, res1$gamma)

Qstat.marginal = function(t, Z, y, mvec, mvec0,mvec1, G, K0, d0, K1, d1, mlist, gamma){      
  # based on the estimation, calculate Q
  
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  ngamma=length(gamma)
  
  #ngamma = (K0+d0+1+K1+d1+1+1)
  arsumg = matrix(rep(0,2*ngamma),nrow=2*ngamma)
  arsumc = matrix(rep(0,2*ngamma*2*ngamma),nrow=2*ngamma)
  #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
  arsumgfirstdev1 = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
  #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
  
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  
  
  for(i in 1:N){
    
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    yi = y[seq]              
    Zi = Z[seq,,drop=FALSE]
    ti=t[seq]
    B0i = myspline(ti, K0, d0)
    B1i = myspline(ti, K1, d1) 
    
    G_exti = G_ext[seq,]   		    		  
    Bi = cbind(B0i, Zi, B1i*G_exti)
    
    mui = Bi%*%gamma
    mudoti = Bi
    
    gi=NULL
    firstdev=NULL
    for(j in 1:J){
      wi = t(mudoti) %*% mlist[[j]][[i]]
      gi=c(gi,wi %*% (yi-mui))
      firstdev=rbind(firstdev,- wi %*% mudoti)#first dev of gij
    }
    arsumc = arsumc + gi %*% t(gi) 
    arsumg = arsumg + gi 
    arsumgfirstdev1 = arsumgfirstdev1 + firstdev
    
    
  }
  
  arcinv1 = ginv(arsumc/N)
  
  Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N) 
  
  Q
  
} 



marginal_test <- function(t, Z, y,mvec, mvec0,mvec1, G, K0, d0, K1, d1,mlist){
  result1 = Qstat_marginal_full(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1,mlist, 0)
  #gamma1 = result1$gamma
  Q1 = result1$QN#Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma1)
  
  result0 = Qstat_marginal_null(t, Z, y, mvec, mvec0,mvec1,K0, d0,mlist, 0)
  Q00=result0$QN
  # gamma_temp = result0[[2]]
  # gamma0 = c(gamma_temp, rep(0,2*(K1+d1+1)))
  # Q0 = Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma0)
  cat("1st version:")
  print(Q00)
  
  df_test = (K1+d1+1)
  gamma0=c(result0$gamma,rep(0,df_test))
  Q0=Qstat.marginal(t, Z, y, mvec, mvec0,mvec1, G, K0, d0, K1, d1, mlist, gamma0)
  cat("2nd version:")
  print(Q0)
  test = Q0 - Q1
  
  # cutoff = qchisq(0.95, df_test)
  # 
  # rej = (test > cutoff)*1
  # return(rej)
  pvalue=pchisq(test,df_test,lower.tail = FALSE)
  return(list(test=test,pvalue=pvalue,df_test=df_test))
  
}



### Testing for genetic effect --- 2nd version



########################################################
########## Joint Testing  ##############################
########################################################


## Qstat: with input the whole data (t, Z, G, ylist, mvec, mvec0, mvec1), model specification (Klist, dlist),
## as well as the working correlation structure (mlist) and the parameter for estimation.
## the output of the function: Qstat itself without penalization, the first derivative of Qstat
## (which will be used for the gradient decent method for the estimation of the parameters) and
## the second derivative of Qstat (which will be used for the GCV function for the selection of smoothing tuning parameter)

## Here the order of the parameter gamma is intercept, partial linear parameter, and the slope
## And in the GEE we assume the marginal variances are the same so we omit the A_i part since that will
## be canceled in the final formulation of Qstat

# Here the gamma is the free parameters for the null model, i.e., reduced from the full model gamma by deleting the positions with 0 elements
# Here we just need to call the Qstat function for the full model and plug in the expanded gamma by filling the deleted positions with 0's.
# Note that the order of the full gamma vector is: the intercept parameter for model 1. the intercept parameter for model 2, the non-varying parameter for model 1
# the non-varying parameter for model 2, and the slope parameter for model 1 and the slope parameter for model 2

Qstat.null2=function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma){
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  ngamma_ex=dz*L
  for(l in 1:L){
    K0=Klist[[l]][1]
    d0=dlist[[l]][1]  
    K1=Klist[[l]][2]
    d1=dlist[[l]][2]
    ngamma_ex=ngamma_ex+K0+d0+K1+d1+2
  }
  
  
  ngamma=length(gamma)
  arsumg = matrix(rep(0,J*ngamma_ex),nrow=J*ngamma_ex)# no matter for full model or null model, the dimension for g is the same
  arsumc = matrix(rep(0,J*ngamma_ex*J*ngamma_ex),nrow=J*ngamma_ex) # no matter for the full model or null model, the dimension for the cov(g) is the same
  arsumgfirstdev1 = matrix(rep(0,J*ngamma_ex*ngamma),nrow=J*ngamma_ex) # here for the null model, the derivative is taken wrt the reduced parameter so the dimension of the gradient matrix is different
  
  
  for(i in 1:N){
    
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    Zi = Z[seq,,drop=FALSE]
    ti = t[seq]
    Gi=G[seq]
    yi=NULL
    B0il=vector("list",L)
    Zil=vector("list",L)
    B1Gil=vector("list",L)
    for(l in 1:L){
      yi=c(yi,ylist[[l]][seq])
      K0=Klist[[l]][1]
      d0=dlist[[l]][1]  
      K1=Klist[[l]][2]
      d1=dlist[[l]][2]
      B0il[[l]]=myspline(ti,K0,d0)
      G_exti = matrix(rep(Gi,K1+d1+1), nrow=ni, ncol=K1+d1+1)
      B1Gil[[l]]=myspline(ti,K1,d1)*G_exti
      Zil[[l]]=Zi
    }
    
    
    B0i_2 = as.matrix(bdiag(B0il))
    Zi_2 = as.matrix(bdiag(Zil))
    B1i_2 = as.matrix(bdiag(B1Gil))
    
    Bi_2 = cbind(B0i_2, Zi_2, B1i_2)
    Bi_2_null=cbind(B0i_2, Zi_2)
    
    mui = Bi_2_null%*%matrix(gamma,ncol=1) # this is the mu under the null model
    mudoti = Bi_2
    mudoti_null = Bi_2_null # mui should also be the same as Bi_2_null%*%gamma
    
    
    
    gi=NULL
    firstdev=NULL
    for(j in 1:J){
      wi = t(mudoti) %*% mlist[[j]][[i]] 
      gi=c(gi,wi %*% (yi-mui))
      firstdev=rbind(firstdev,-wi %*% mudoti_null)#first dev of gij
    }
    
    arsumc = arsumc + gi %*% t(gi) 
    arsumg = arsumg + gi 
    arsumgfirstdev1 = arsumgfirstdev1 + firstdev
  }
  
  arcinv1 = ginv(arsumc/N) 
  
  Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N)
  
  Qp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N)
  Qpp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  list(Q=Q,Qp=Qp,Qpp=Qpp,arcinv1=arcinv1,arsumgfirstdev1=arsumgfirstdev1)
}


estimation_null2 = function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,lambda){
  #t: a long vector which stacks all of the measurement time points for all individuals
  #Z: a long (tall) matrix which stacks all of the other covariates measured at t
  #ylist: a list of multiple responses, each response vector is a long vector which stacks the measurement at t
  #mvec: a vector of n_i
  #mvec0: a vector of start measuring points of individuals
  #mvec1: a vector of end measuring points of individuals
  #G: a long vector of genetic information
  #Klist: list of number of knots for coefficient functions
  #dlist: list of degrees for coefficient functions
  #mlist: list of basis matrices
  #lambda: smoothing tunning parameter
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  D0=NULL
  gamma.old0=NULL
  alpha.old=NULL
  for(l in 1:L){
    K0=Klist[[l]][1]
    d0=dlist[[l]][1]  
    K1=Klist[[l]][2]
    d1=dlist[[l]][2]  
    temp0=c(rep(0,d0+1),rep(1,K0))
    D0=c(D0,temp0)
    B0=myspline(t,K0,d0)
    B_temp = cbind(B0, Z)
    gamma.l=ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%ylist[[l]]
    gamma.old0=c(gamma.old0,gamma.l[1:(K0+d0+1)])
    alpha.old=c(alpha.old,gamma.l[K0+d0+1+(1:dz)])
  }
  
  D = diag(c(D0, rep(0,L*dz)))  
  gamma.old=c(gamma.old0,alpha.old)##simple LS, initial value 
  ngamma=length(gamma.old)
  
  run = 0  
  while(run <= 50){
    run = run+1
    res=Qstat.null2(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma.old)
    arqif2dev1=res$Qpp/N+ 2*lambda*D	
    invarqif2dev1 = ginv(arqif2dev1)
    arqif1dev1=res$Qp/N+2*lambda*D%*%gamma.old
    gamma.new = gamma.old -  invarqif2dev1%*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    gamma.old = gamma.new
    #print(gamma.old)
  } #loop for gamma
  res=Qstat.null2(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma.old) 
  QNdotdot = res$Qpp
  QN = res$Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  
  list(GCV=GCV, gamma=gamma.old, QN=QN)
  
} 

joint_test2 <- function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist){
  result1 = estimation_full(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,0)
  #gamma1 = result1$gamma
  Q1 = result1$QN#Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma1)
  
  # Klist0=list(Klist[[1]][1],Klist[[2]][1])
  # dlist0=list(dlist[[1]][1],dlist[[2]][1])
  result0 = estimation_null2(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist, 0)
  Q0=result0$QN
  # gamma_temp = result0[[2]]
  # gamma0 = c(gamma_temp, rep(0,2*(K1+d1+1)))
  # Q0 = Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma0)
  df_test =Klist[[1]][2]+Klist[[2]][2]+dlist[[1]][2]+dlist[[2]][2]+2 #2*(K1+d1+1)-2
  
  # gamma0=c(result0$gamma,rep(0,df_test))
  # result00=Qstat(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma0)
  # #print(result00)
  # Q0=result00$Q
  
  test = Q0 - Q1
  
  
  #cutoff = qchisq(0.95, df_test)
  
  #rej = (test > cutoff)*1
  pvalue=pchisq(test,df_test,lower.tail = FALSE)
  return(list(test=test,pvalue=pvalue,df_test=df_test))
}



########################################################
########## Marginal Testing  ###########################
########################################################
Qstat_marginal_full = function(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1, mlist,lambda){
  dz=ncol(Z)
  D1 = c(rep(0,d0+1), rep(1,K0))
  D2 = c(rep(0,d1+1), rep(1,K1))
  D = diag(c(D1, rep(0,dz), D2))    
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  
  ##step 0: initial values
  
  B0 = myspline(t, K0, d0)
  B1 = myspline(t, K1, d1)               
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B_temp = cbind(B0, Z, B1*G_ext)
  
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value 
  
  run = 0  
  while(run <= 50){
    
    run = run+1
    
    ngamma = (K0+d0+1+K1+d1+1+dz)
    arsumg = matrix(rep(0,J*ngamma),nrow=J*ngamma)
    arsumc = matrix(rep(0,J*ngamma*2*ngamma),nrow=J*ngamma)
    #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumgfirstdev1 = matrix(rep(0,J*ngamma*ngamma),nrow=J*ngamma)
    #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi =Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)
      B1i = myspline(ti, K1, d1) 
      
      G_exti = G_ext[seq,]   		 
      
      Bi = cbind(B0i, Zi, B1i*G_exti)
      
      mui = Bi%*%gamma.old
      mudoti = Bi
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% mlist[[j]][[i]]
        gi=c(gi,wi %*% (yi-mui))
        firstdev=rbind(firstdev,- wi %*% mudoti)#first dev of gij
      }
      
      arsumc = arsumc + gi %*% t(gi) 
      arsumg = arsumg + gi 
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc/N)
    
    Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N) 
    
    arqif1dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N))/N + 2*lambda*D%*%gamma.old
    arqif2dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1)/N)/N + 2*lambda*D	
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  QN = Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  k=length(gamma.old)
  BIC = QN+(J-1)*k*log(N)
  
  list(GCV=GCV, QN=QN,gamma=gamma.old,BIC=BIC)
  
} 

Qstat_marginal_null2 = function(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1, mlist,lambda){
  dz=ncol(Z)
  D1 = c(rep(0,d0+1), rep(1,K0))
  D = diag(c(D1, rep(0,dz)))
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  
  ##step 0: initial values
  B0 = myspline(t, K0, d0)
  B_temp = cbind(B0, Z)
  
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value 
  ngamma = length(gamma.old)
  ngamma_ex=(K0+d0+1+K1+d1+1+dz)
  run = 0  
  while(run <= 50){
    run = run+1
    arsumg = matrix(rep(0,J*ngamma_ex),nrow=J*ngamma_ex)
    arsumc = matrix(rep(0,J*ngamma_ex*J*ngamma_ex),nrow=J*ngamma_ex)
    #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumgfirstdev1 = matrix(rep(0,J*ngamma_ex*ngamma),nrow=J*ngamma_ex)
    #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi =Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)
      B1i = myspline(ti, K1, d1) 
      
      G_exti = G_ext[seq,]   		 
      
      Bi = cbind(B0i, Zi, B1i*G_exti)
      
      #mui = Bi%*%gamma.old
      mudoti = Bi
      mudoti_null = cbind(B0i, Zi)
      mui = mudoti_null%*%gamma.old
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% mlist[[j]][[i]]
        gi=c(gi,wi %*% (yi-mui))
        firstdev=rbind(firstdev,- wi %*% mudoti_null)#first dev of gij
      }
      
      arsumc = arsumc + gi %*% t(gi) 
      arsumg = arsumg + gi 
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc/N)
    
    Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N) 
    
    arqif1dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N))/N + 2*lambda*D%*%gamma.old
    arqif2dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1)/N)/N + 2*lambda*D	
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  QN = Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  k=length(gamma.old)
  BIC = QN+(J-1)*k*log(N)
  
  list(GCV=GCV, QN=QN,gamma=gamma.old,BIC=BIC)
  
} 






marginal_test2 <- function(t, Z, y,mvec, mvec0,mvec1, G, K0, d0, K1, d1,mlist){
  result1 = Qstat_marginal_full(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1,mlist, 0)
  #gamma1 = result1$gamma
  Q1 = result1$QN#Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma1)
  
  result0 = Qstat_marginal_null2(t, Z, y, mvec, mvec0,mvec1,G,K0, d0,K1,d1,mlist, 0)
  Q0=result0$QN
  # gamma_temp = result0[[2]]
  # gamma0 = c(gamma_temp, rep(0,2*(K1+d1+1)))
  # Q0 = Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma0)
  cat("1st version:")
  print(Q0)
  
  df_test = (K1+d1+1)
  # gamma0=c(result0$gamma,rep(0,df_test))
  # Q0=Qstat.marginal(t, Z, y, mvec, mvec0,mvec1, G, K0, d0, K1, d1, mlist, gamma0)
  # cat("2nd version:")
  # print(Q0)
  test = Q0 - Q1
  
  # cutoff = qchisq(0.95, df_test)
  # 
  # rej = (test > cutoff)*1
  # return(rej)
  pvalue=pchisq(test,df_test,lower.tail = FALSE)
  return(list(test=test,pvalue=pvalue,df_test=df_test))
  
}


### Testing for Interaction


########################################################
########## Joint Testing  ##############################
########################################################


## Qstat: with input the whole data (t, Z, G, ylist, mvec, mvec0, mvec1), model specification (Klist, dlist),
## as well as the working correlation structure (mlist) and the parameter for estimation.
## the output of the function: Qstat itself without penalization, the first derivative of Qstat
## (which will be used for the gradient decent method for the estimation of the parameters) and
## the second derivative of Qstat (which will be used for the GCV function for the selection of smoothing tuning parameter)

## Here the order of the parameter gamma is intercept, partial linear parameter, and the slope
## And in the GEE we assume the marginal variances are the same so we omit the A_i part since that will
## be canceled in the final formulation of Qstat

# Here the gamma is the free parameters for the null model, i.e., reduced from the full model gamma by deleting the positions with 0 elements
# Here we just need to call the Qstat function for the full model and plug in the expanded gamma by filling the deleted positions with 0's.
# Note that the order of the full gamma vector is: the intercept parameter for model 1. the intercept parameter for model 2, the non-varying parameter for model 1
# the non-varying parameter for model 2, and the slope parameter for model 1 and the slope parameter for model 2

Qstat.null.interaction=function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma){
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  pos=0
  for(l in 1:L){
    K0=Klist[[l]][1]
    d0=dlist[[l]][1]  
    pos=pos+K0+d0+1+dz
  }
  gamma_ex=gamma[1:pos]
  for(l in 1:L){
    K1=Klist[[l]][2]
    d1=dlist[[l]][2]  
    gamma_ex=c(gamma_ex,gamma[pos+l],rep(0,K1+d1))
  }
  
  
  
  ngamma=length(gamma)
  ngamma_ex=length(gamma_ex)
  arsumg = matrix(rep(0,J*ngamma_ex),nrow=J*ngamma_ex)# no matter for full model or null model, the dimension for g is the same
  arsumc = matrix(rep(0,J*ngamma_ex*J*ngamma_ex),nrow=J*ngamma_ex) # no matter for the full model or null model, the dimension for the cov(g) is the same
  arsumgfirstdev1 = matrix(rep(0,J*ngamma_ex*ngamma),nrow=J*ngamma_ex) # here for the null model, the derivative is taken wrt the reduced parameter so the dimension of the gradient matrix is different
  
  
  for(i in 1:N){
    
    seq = c(mvec0[i]:mvec1[i])
    ni = mvec[i]
    Zi = Z[seq,,drop=FALSE]
    ti = t[seq]
    Gi=G[seq]
    yi=NULL
    B0il=vector("list",L)
    Zil=vector("list",L)
    B1Gil=vector("list",L)
    B1Gil_null=vector("list",L)
    for(l in 1:L){
      yi=c(yi,ylist[[l]][seq])
      K0=Klist[[l]][1]
      d0=dlist[[l]][1]  
      K1=Klist[[l]][2]
      d1=dlist[[l]][2]
      null_factor_matrix=matrix(c(1,rep(0,K1+d1)),ncol=1)
      B0il[[l]]=myspline(ti,K0,d0)
      G_exti = matrix(rep(Gi,K1+d1+1), nrow=ni, ncol=K1+d1+1)
      B1Gil[[l]]=myspline(ti,K1,d1)*G_exti
      Zil[[l]]=Zi
      B1Gil_null[[l]]=B1Gil[[l]]%*%null_factor_matrix
    }
    
    
    B0i_2 = as.matrix(bdiag(B0il))
    Zi_2 = as.matrix(bdiag(Zil))
    B1i_2 = as.matrix(bdiag(B1Gil))
    B1i_2_null=as.matrix(bdiag(B1Gil_null))
    
    Bi_2 = cbind(B0i_2, Zi_2, B1i_2)
    Bi_2_null=cbind(B0i_2, Zi_2, B1i_2_null)
    
    mui = Bi_2%*%gamma_ex # this is the mu under the null model
    mudoti = Bi_2
    mudoti_null = Bi_2_null # mui should also be the same as Bi_2_null%*%gamma
    
    
    
    gi=NULL
    firstdev=NULL
    for(j in 1:J){
      wi = t(mudoti) %*% mlist[[j]][[i]] 
      gi=c(gi,wi %*% (yi-mui))
      firstdev=rbind(firstdev,-wi %*% mudoti_null)#first dev of gij
    }
    
    arsumc = arsumc + gi %*% t(gi) 
    arsumg = arsumg + gi 
    arsumgfirstdev1 = arsumgfirstdev1 + firstdev
  }
  
  arcinv1 = ginv(arsumc/N) 
  
  Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N)
  
  Qp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N)
  Qpp = 2*N*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  list(Q=Q,Qp=Qp,Qpp=Qpp,arcinv1=arcinv1,arsumgfirstdev1=arsumgfirstdev1)
}


estimation_null_interaction = function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,lambda){
  #t: a long vector which stacks all of the measurement time points for all individuals
  #Z: a long (tall) matrix which stacks all of the other covariates measured at t
  #ylist: a list of multiple responses, each response vector is a long vector which stacks the measurement at t
  #mvec: a vector of n_i
  #mvec0: a vector of start measuring points of individuals
  #mvec1: a vector of end measuring points of individuals
  #G: a long vector of genetic information
  #Klist: list of number of knots for coefficient functions
  #dlist: list of degrees for coefficient functions
  #mlist: list of basis matrices
  #lambda: smoothing tunning parameter
  L=length(ylist)
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  dz=ncol(Z)
  D0=NULL
  D1=NULL
  D1_null=NULL
  gamma.old0=NULL
  gamma.old1=NULL
  alpha.old=NULL
  for(l in 1:L){
    K0=Klist[[l]][1]
    d0=dlist[[l]][1]  
    K1=Klist[[l]][2]
    d1=dlist[[l]][2]  
    null_factor_matrix=matrix(c(1,rep(0,K1+d1)),ncol=1)
    temp0=c(rep(0,d0+1),rep(1,K0))
    temp1=c(rep(0,d1+1),rep(1,K1))
    temp1_null=0
    D0=c(D0,temp0)
    D1=c(D1,temp1)
    D1_null=c(D1_null,temp1_null)
    B0=myspline(t,K0,d0)
    G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
    B1G=myspline(t,K1,d1)*G_ext
    B1G_null=B1G%*%null_factor_matrix
    B_temp = cbind(B0, Z, B1G_null)
    gamma.l=ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%ylist[[l]]
    gamma.old0=c(gamma.old0,gamma.l[1:(K0+d0+1)])
    alpha.old=c(alpha.old,gamma.l[K0+d0+1+(1:dz)])
    gamma.old1=c(gamma.old1,gamma.l[(K0+d0+2+dz):length(gamma.l)])
  }
  
  D = diag(c(D0, rep(0,L*dz), D1))  
  D_null=diag(c(D0,rep(0,L*dz), D1_null))
  gamma.old=c(gamma.old0,alpha.old,gamma.old1)##simple LS, initial value 
  ngamma=length(gamma.old)
  
  run = 0  
  while(run <= 50){
    run = run+1
    res=Qstat.null.interaction(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma.old)
    arqif2dev1=res$Qpp/N+ 2*lambda*D_null	
    invarqif2dev1 = ginv(arqif2dev1)
    arqif1dev1=res$Qp/N+2*lambda*D_null%*%gamma.old
    gamma.new = gamma.old -  invarqif2dev1%*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    gamma.old = gamma.new
    #print(gamma.old)
  } #loop for gamma
  res=Qstat.null.interaction(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma.old) 
  QNdotdot = res$Qpp
  QN = res$Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D_null)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  
  list(GCV=GCV, gamma=gamma.old, QN=QN)
  
} 

joint_test_interaction <- function(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist){
  result1 = estimation_full(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,0)
  #gamma1 = result1$gamma
  Q1 = result1$QN#Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma1)
  
  # Klist0=list(Klist[[1]][1],Klist[[2]][1])
  # dlist0=list(dlist[[1]][1],dlist[[2]][1])
  result0 = estimation_null_interaction(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist, 0)
  Q0=result0$QN
  # gamma_temp = result0[[2]]
  # gamma0 = c(gamma_temp, rep(0,2*(K1+d1+1)))
  # Q0 = Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma0)
  df_test =Klist[[1]][2]+Klist[[2]][2]+dlist[[1]][2]+dlist[[2]][2] #2*(K1+d1+1)-2
  
  # gamma0=c(result0$gamma,rep(0,df_test))
  # result00=Qstat(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist,gamma0)
  # #print(result00)
  # Q0=result00$Q
  
  test = Q0 - Q1
  
  
  #cutoff = qchisq(0.95, df_test)
  
  #rej = (test > cutoff)*1
  pvalue=pchisq(test,df_test,lower.tail = FALSE)
  return(list(test=test,pvalue=pvalue,df_test=df_test))
}


########################################################
########## Marginal Testing  ###########################
########################################################
Qstat_marginal_full = function(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1, mlist,lambda){
  dz=ncol(Z)
  D1 = c(rep(0,d0+1), rep(1,K0))
  D2 = c(rep(0,d1+1), rep(1,K1))
  D = diag(c(D1, rep(0,dz), D2))    
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  
  ##step 0: initial values
  
  B0 = myspline(t, K0, d0)
  B1 = myspline(t, K1, d1)               
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B_temp = cbind(B0, Z, B1*G_ext)
  
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value 
  
  run = 0  
  while(run <= 50){
    
    run = run+1
    
    ngamma = (K0+d0+1+K1+d1+1+dz)
    arsumg = matrix(rep(0,J*ngamma),nrow=J*ngamma)
    arsumc = matrix(rep(0,J*ngamma*2*ngamma),nrow=J*ngamma)
    #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumgfirstdev1 = matrix(rep(0,J*ngamma*ngamma),nrow=J*ngamma)
    #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi =Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)
      B1i = myspline(ti, K1, d1) 
      
      G_exti = G_ext[seq,]   		 
      
      Bi = cbind(B0i, Zi, B1i*G_exti)
      
      mui = Bi%*%gamma.old
      mudoti = Bi
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% mlist[[j]][[i]]
        gi=c(gi,wi %*% (yi-mui))
        firstdev=rbind(firstdev,- wi %*% mudoti)#first dev of gij
      }
      
      arsumc = arsumc + gi %*% t(gi) 
      arsumg = arsumg + gi 
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc/N)
    
    Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N) 
    
    arqif1dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N))/N + 2*lambda*D%*%gamma.old
    arqif2dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1)/N)/N + 2*lambda*D	
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  QN = Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  k=length(gamma.old)
  BIC = QN+(J-1)*k*log(N)
  
  list(GCV=GCV, QN=QN,gamma=gamma.old,BIC=BIC)
  
} 

Qstat_marginal_null_interaction = function(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1, mlist,lambda){
  dz=ncol(Z)
  D1 = c(rep(0,d0+1), rep(1,K0))
  D2 = 0
  D = diag(c(D1, rep(0,dz), D2))
  Nobs=length(t)
  N=length(mvec)
  J=length(mlist)
  
  ##step 0: initial values
  null_factor_matrix=matrix(c(1,rep(0,K1+d1)),ncol=1)
  B0 = myspline(t, K0, d0)
  B1 = myspline(t, K1, d1)               
  G_ext = matrix(rep(G,K1+d1+1), nrow=Nobs, ncol=K1+d1+1)
  B_temp = cbind(B0, Z, (B1*G_ext)%*%null_factor_matrix)
  
  
  gamma.old = ginv(t(B_temp)%*%B_temp)%*%t(B_temp)%*%y  ##simple LS, initial value 
  ngamma = length(gamma.old)
  ngamma_ex=(K0+d0+1+K1+d1+1+dz)
  run = 0  
  while(run <= 50){
    run = run+1
    arsumg = matrix(rep(0,J*ngamma_ex),nrow=J*ngamma_ex)
    arsumc = matrix(rep(0,J*ngamma_ex*J*ngamma_ex),nrow=J*ngamma_ex)
    #gi = matrix(rep(0,2*ngamma),nrow=2*ngamma)
    arsumgfirstdev1 = matrix(rep(0,J*ngamma_ex*ngamma),nrow=J*ngamma_ex)
    #firstdev = matrix(rep(0,2*ngamma*ngamma),nrow=2*ngamma)
    
    for(i in 1:N){
      
      seq = c(mvec0[i]:mvec1[i])
      ni = mvec[i]
      yi = y[seq]
      Zi =Z[seq,,drop=FALSE]
      ti=t[seq]
      B0i = myspline(ti, K0, d0)
      B1i = myspline(ti, K1, d1) 
      
      G_exti = G_ext[seq,]   		 
      
      Bi = cbind(B0i, Zi, B1i*G_exti)
      
      #mui = Bi%*%gamma.old
      mudoti = Bi
      mudoti_null = cbind(B0i, Zi, (B1i*G_exti)%*%null_factor_matrix)
      mui = mudoti_null%*%gamma.old
      
      gi=NULL
      firstdev=NULL
      for(j in 1:J){
        wi = t(mudoti) %*% mlist[[j]][[i]]
        gi=c(gi,wi %*% (yi-mui))
        firstdev=rbind(firstdev,- wi %*% mudoti_null)#first dev of gij
      }
      
      arsumc = arsumc + gi %*% t(gi) 
      arsumg = arsumg + gi 
      arsumgfirstdev1 = arsumgfirstdev1 + firstdev
      
      
    }
    
    arcinv1 = ginv(arsumc/N)
    
    Q = N*t(arsumg/N) %*% arcinv1 %*% (arsumg/N) 
    
    arqif1dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumg/N))/N + 2*lambda*D%*%gamma.old
    arqif2dev1 = (2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1)/N)/N + 2*lambda*D	
    
    invarqif2dev1 = ginv(arqif2dev1)
    
    gamma.new = gamma.old - invarqif2dev1 %*% arqif1dev1
    gammadiff = max(abs(gamma.new - gamma.old))
    if(gammadiff<1e-6){break}
    
    gamma.old = gamma.new
    #print(gamma.old)
    
  } #loop for gamma
  
  # based on the estimation, calculate Q, covariance
  QNdotdot = 2*t(arsumgfirstdev1/N) %*% arcinv1 %*% (arsumgfirstdev1/N)
  QN = Q
  df = sum(diag(ginv(QNdotdot + 2*N*lambda*D)%*%QNdotdot))
  GCV = (QN/N)/(1-df/N)^2
  
  k=length(gamma.old)
  BIC = QN+(J-1)*k*log(N)
  
  list(GCV=GCV, QN=QN,gamma=gamma.old,BIC=BIC)
  
} 






marginal_test_interaction <- function(t, Z, y,mvec, mvec0,mvec1, G, K0, d0, K1, d1,mlist){
  result1 = Qstat_marginal_full(t, Z, y, mvec, mvec0,mvec1,G, K0, d0, K1, d1,mlist, 0)
  #gamma1 = result1$gamma
  Q1 = result1$QN#Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma1)
  
  result0 = Qstat_marginal_null_interaction(t, Z, y, mvec, mvec0,mvec1,G,K0, d0,K1,d1,mlist, 0)
  Q0=result0$QN
  # gamma_temp = result0[[2]]
  # gamma0 = c(gamma_temp, rep(0,2*(K1+d1+1)))
  # Q0 = Qstat(t, Z, y1_0, y2_0, G, K0, d0, K1, d1, gamma0)
  cat("1st version:")
  print(Q0)
  
  df_test = (K1+d1+1)-1
  # gamma0=c(result0$gamma,rep(0,df_test))
  # Q0=Qstat.marginal(t, Z, y, mvec, mvec0,mvec1, G, K0, d0, K1, d1, mlist, gamma0)
  # cat("2nd version:")
  # print(Q0)
  test = Q0 - Q1
  
  # cutoff = qchisq(0.95, df_test)
  # 
  # rej = (test > cutoff)*1
  # return(rej)
  pvalue=pchisq(test,df_test,lower.tail = FALSE)
  return(list(test=test,pvalue=pvalue,df_test=df_test))
  
}


obs_sim=function(N,m,pr,rho_vec,tau,sigma_vec,alpha_vec,nul,nul2){
  # N = 200
  # m = 10 
  
  id = rep(seq(1,N,1),rep(m,N))
  mvec = rep(m,N)
  Nobs = sum(mvec)
  mvec0 = cumsum(c(1, mvec[-N])) #start loc for each id
  mvec1 = cumsum(mvec) #end loc for each id
  
  Z = matrix(runif(Nobs,0,1),ncol=1)
  
  G0 = sample(c(0,1,2), size = N, replace = T, prob=c(pr^2, 2*pr*(1-pr), (1-pr)^2))
  G = rep(G0,mvec) 
  
  # t0 = seq(1,m,1)/m
  # t = rep(t0,N)    
  t=c(replicate(N,sort(runif(m))))
  
  rho1=rho_vec[1]
  rho2=rho_vec[2]
  rho12=rho_vec[3]
  sigma1=sigma_vec[1]
  sigma2=sigma_vec[2]
  
  
  R1 = matrix(rho1,m,m)
  R2 = matrix(rho2,m,m) 
  R12 = matrix(rho12,m,m) 
  diag(R1)=1
  diag(R2)=1
  diag(R12)=tau
  
  sigma11 = (sigma1^2)*R1
  sigma22 = (sigma2^2)*R2
  sigma12 = sigma1*sigma2*R12
  
  block1 = cbind(sigma11,sigma12)
  block2 = cbind(sigma12,sigma22)
  VAR = rbind(block1,block2)
  
  
  # covdata=huge.generator(N,d=2*m,graph="cluster")
  # temp=covdata$data
  # cat("covariance of the error:")
  # print(covdata$sigma)
  temp = rmvnorm(N,mean=rep(0,2*m), sigma=VAR)
  eps1 = as.vector(t(temp[,1:m]))
  eps2 = as.vector(t(temp[,(m+1):(2*m)]))
  
  beta1_0 = 0.5*cos(2*pi*t)
  beta1_1 = nul*sin(pi*(t-0.2))+nul2
  
  beta2_0 = sin(pi*t)-0.5
  beta2_1 = nul*cos(pi*t-0.8)+nul2
  
  
  
  alpha1=alpha_vec[1]
  alpha2=alpha_vec[2]
  
  #eps_add=rlnorm(2*N*m, 0,1)-exp(1^2/2)
  
  y1 = beta1_0 + beta1_1*G + alpha1*Z + eps1#+eps_add[1:(N*m)]    # full model
  y2 = beta2_0 + beta2_1*G + alpha2*Z + eps2#+eps_add[-(1:(N*m))]    
  
  
  ylist=list(y1,y2)
  return(list(ylist=ylist,t=t,G=G,Z=Z,mvec=mvec,mvec0=mvec0,mvec1=mvec1))
}
