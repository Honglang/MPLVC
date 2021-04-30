#########################################################################
####################### Data Generation #################################
#########################################################################

#######################################################################
#### For one snp, the longitudinal data if formated into long form  ###
#######################################################################

library(tidyverse)
library(stringr)
library("MASS")
library("mvtnorm")
library("Matrix")
library(huge)


####################################################
### Covariance Structure ###########################
####################################################
rho1=0.5
rho2=-0.5
rho12=0.1
tau=-0.5
sigma1=sqrt(0.1)
sigma2=sqrt(0.1)
m=10

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

eigen(VAR)$values

####################################################
### Conduct testing for interaction#################
####################################################
## With the null that both of the two slopes are constants, then if fail to reject the null, 
## it indicate there is no interaction effect between the genetic and the covariate for the joint response



###############################################
### Data Generation  ##########################
###############################################
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

N=200
m=5
pr=0.3
nul=0.001
nul2=0


data=obs_sim(N,m,pr,rho_vec=c(0.5,0.5,0.2),tau=0.2,sigma_vec=c(sqrt(0.1),sqrt(0.1)),alpha_vec=c(0.2,0.3),nul,nul2)
ylist=data$ylist
t=data$t
G=data$G
Z=data$Z
mvec=data$mvec
mvec0=data$mvec0
mvec1=data$mvec1

source("./mplvc_functions.R")

###################################################
####  Model selection    ##########################
###################################################
cor=1
method=0
df=2

if(cor==1){
  ### corr matrix: AR(1) #######################################
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    m0[[i]] = diag(ni)
    m1[[i]] = matrix(rep(0,ni*ni),ni)
    for (k in 1:(ni)) {
      for (l in 1:(ni)) {
        if (abs(k-l)==1) m1[[i]][k,l] =1
      }
    }
  }
  mlist1=list(m0,m1)
  
  ### corr matrix: ar(1) #####
  m0 = list()
  m1 = list()
  m2 = list()
  m3 = list()
  for (i in 1:N){
    ni = mvec[i]
    I=diag(ni)
    O=matrix(0,ni,ni)
    U2=O
    for (k in 1:(ni)) {
      for (l in 1:(ni)) {
        if (abs(k-l)==1) U2[k,l] =1
      }
    }
    
    ## corr matrix: exchangable  
    m0[[i]] = diag(2*ni)           
    m1[[i]] = rbind(cbind(O,I),cbind(I,O))
    m2[[i]] = rbind(cbind(U2,O),cbind(O,U2))
    m3[[i]] = rbind(cbind(O,U2),cbind(U2,O))
  }
  mlist2=list(m0,m1,m2,m3) 
}

if(cor==2){
  ### corr matrix: exchangable #####
  m0 = list()
  m1 = list()
  for (i in 1:N){
    ni = mvec[i]
    ## corr matrix: exchangable  
    m0[[i]] = diag(ni)           
    m1[[i]] = matrix(rep(1,ni*ni),ni) - m0[[i]]
  }
  mlist1=list(m0,m1) 
  
  ### corr matrix: exchangable #####
  m0 = list()
  m1 = list()
  m2 = list()
  m3 = list()
  for (i in 1:N){
    ni = mvec[i]
    I=diag(ni)
    O=matrix(0,ni,ni)
    U1=matrix(1,ni,ni)-I
    ## corr matrix: exchangable  
    m0[[i]] = diag(2*ni)           
    m1[[i]] = rbind(cbind(O,I),cbind(I,O))
    m2[[i]] = rbind(cbind(U1,O),cbind(O,U1))
    m3[[i]] = rbind(cbind(O,U1),cbind(U1,O))
  }
  mlist2=list(m0,m1,m2,m3) 
  
}



# method 0 is just use the rule of thumb
if(method==0){
  Klist=vector("list",2)
  dlist=vector("list",2)
  Klist[[1]]=rep(floor(0.6*N^(1/5)),2)  
  Klist[[2]]=rep(floor(0.6*N^(1/5)),2)  
  dlist[[1]]=rep(df,2)
  dlist[[2]]=rep(df,2)
}

# method 1 is selecting the tunning parameters for each model separately
# and for the intercept and slope functions, we use different orders and numbers of knots
if(method==1){
  d0_list = df
  d1_list = df
  K1_list = 1:5
  K0_list = 1:5
  all_list = expand.grid(d0_list,d1_list,K0_list)
  Klist=vector("list",2)
  dlist=vector("list",2)
  
  BIC=apply(all_list,1,function(x) BIC_marginal(t, Z, ylist[[1]], G, mvec, mvec0,mvec1,mlist1,x[3], x[1], x[3], x[2]))
  
  min.row = which.min(BIC)
  
  d0 = all_list[min.row,1]
  d1 = all_list[min.row,2]
  K0 = all_list[min.row,3]
  K1 = all_list[min.row,3]
  Klist[[1]]=c(K0,K1)
  dlist[[1]]=c(d0,d1)
  
  BIC=apply(all_list,1,function(x) BIC_marginal(t, Z, ylist[[2]], G, mvec, mvec0,mvec1,mlist1,x[3], x[1], x[3], x[2]))
  
  min.row = which.min(BIC)
  
  d0 = all_list[min.row,1]
  d1 = all_list[min.row,2]
  K0 = all_list[min.row,3]
  K1 = all_list[min.row,3]
  Klist[[2]]=c(K0,K1)
  dlist[[2]]=c(d0,d1)
  
  Klist
  dlist
}

# method 2 selecting the tunning parameters for each model separately
# and for each model, we first select the order and number of knots for the intercept function 
# and then based on this to select the order and number of knots for the slope function.
if(method==2){
  d0_list = df
  d1_list = df
  K1_list = 0:3
  K0_list = 0:3
  all_list_0 = expand.grid(d0_list,K0_list)
  all_list_1 = expand.grid(d1_list,K1_list)
  Klist=vector("list",2)
  dlist=vector("list",2)
  
  BIC_0=apply(all_list_0,1,function(x) BIC_marginal_null(t, Z, ylist[[1]], mvec, mvec0,mvec1,mlist1,x[2], x[1]))
  min.row = which.min(BIC_0)
  
  d0 = all_list_0[min.row,1]
  K0 = all_list_0[min.row,2]
  
  
  BIC_1=apply(all_list_1,1,function(x) BIC_marginal(t, Z, ylist[[1]], G, mvec, mvec0,mvec1,mlist1,K0, d0, x[2], x[1]))
  
  min.row = which.min(BIC_1)
  
  d1 = all_list_1[min.row,1]
  K1 = all_list_1[min.row,2]
  Klist[[1]]=c(K0,K1)
  dlist[[1]]=c(d0,d1)
  
  BIC_0=apply(all_list_0,1,function(x) BIC_marginal_null(t, Z, ylist[[2]], mvec, mvec0,mvec1,mlist1,x[2], x[1]))
  min.row = which.min(BIC_0)
  
  d0 = all_list_0[min.row,1]
  K0 = all_list_0[min.row,2]
  
  
  BIC_1=apply(all_list_1,1,function(x) BIC_marginal(t, Z, ylist[[2]], G, mvec, mvec0,mvec1,mlist1,K0, d0, x[2], x[1]))
  
  min.row = which.min(BIC_1)
  
  d1 = all_list_1[min.row,1]
  K1 = all_list_1[min.row,2]
  
  Klist[[2]]=c(K0,K1)
  dlist[[2]]=c(d0,d1)
  
  Klist
  dlist  
}

K01=Klist[[1]][1];K11=Klist[[1]][2];K02=Klist[[2]][1];K12=Klist[[2]][2]
d01=dlist[[1]][1];d11=dlist[[1]][2];d02=dlist[[2]][1];d12=dlist[[2]][2]



GCV = function(lam){estimation_full(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist2,lam)$GCV}
(lamb_est = golden(GCV, 0.01, 2))
#lamb_est=1
(est_result=estimation_full(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist2,lamb_est))

cov_gamma=estimation_cov(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist2,lamb_est,est_result$gamma)




w = seq(0.01,0.99,0.01)
gamma.old=est_result$gamma

gamma_index_s=cumsum(c(1,K01+d01+1,K02+d02+1,2,K11+d11+1))
gamma_index_e=cumsum(c(K01+d01+1,K02+d02+1,2,K11+d11+1,K12+d12+1))

cov_gamma01 = cov_gamma[gamma_index_s[1]:gamma_index_e[1],gamma_index_s[1]:gamma_index_e[1]]
cov_gamma02 = cov_gamma[gamma_index_s[2]:gamma_index_e[2],gamma_index_s[2]:gamma_index_e[2]]
cov_gamma11 = cov_gamma[gamma_index_s[4]:gamma_index_e[4],gamma_index_s[4]:gamma_index_e[4]]
cov_gamma12 = cov_gamma[gamma_index_s[5]:gamma_index_e[5],gamma_index_s[5]:gamma_index_e[5]]


B01=myspline(w,K01,d01)
B02=myspline(w,K02,d02)
B11=myspline(w,K11,d11)
B12=myspline(w,K12,d12)
beta01.est = B01%*%gamma.old[gamma_index_s[1]:gamma_index_e[1]]
beta02.est = B02%*%gamma.old[gamma_index_s[2]:gamma_index_e[2]]
beta11.est = B11%*%gamma.old[gamma_index_s[4]:gamma_index_e[4]]
beta12.est = B12%*%gamma.old[gamma_index_s[5]:gamma_index_e[5]]


m01SE = sqrt(diag(B01%*%cov_gamma01%*%t(B01)))
m02SE = sqrt(diag(B02%*%cov_gamma02%*%t(B02)))
m11SE = sqrt(diag(B11%*%cov_gamma11%*%t(B11)))
m12SE = sqrt(diag(B12%*%cov_gamma12%*%t(B12)))

CBl01 = beta01.est - 1.96*m01SE
CBu01 = beta01.est + 1.96*m01SE
CBl02 = beta02.est - 1.96*m02SE
CBu02 = beta02.est + 1.96*m02SE
CBl11 = beta11.est - 1.96*m11SE
CBu11 = beta11.est + 1.96*m11SE
CBl12 = beta12.est - 1.96*m12SE
CBu12 = beta12.est + 1.96*m12SE


par(mfrow=c(2,2), mar=c(2,2,0.5,0.5), oma=c(3,3,1,1))
plot(w,beta01.est,col="red",type = "l",ylim=c(min(CBl01),max(CBu01)))
lines(w,CBl01,col="blue",lty=2)
lines(w,CBu01,col="blue",lty=2)
legend("topleft", legend=expression(hat(beta)["01"](t)),cex=0.8,bg='lightblue')


plot(w,beta11.est,col="red",type = "l",ylim=c(min(CBl11),max(CBu11)))
lines(w,CBl11,col="blue",lty=2)
lines(w,CBu11,col="blue",lty=2)
legend("topleft", legend=expression(hat(beta)["11"](t)),cex=0.8,bg='lightblue')

plot(w,beta02.est,col="red",type = "l",ylim=c(min(CBl02),max(CBu02)))
lines(w,CBl02,col="blue",lty=2)
lines(w,CBu02,col="blue",lty=2)
legend("bottomleft", legend=expression(hat(beta)["02"](t)),cex=0.8,bg='lightblue')

plot(w,beta12.est,col="red",type = "l",ylim=c(min(CBl12),max(CBu12)))
lines(w,CBl12,col="blue",lty=2)
lines(w,CBu12,col="blue",lty=2)
legend("bottomleft", legend=expression(hat(beta)["12"](t)),cex=0.8,bg='lightblue')


test_result_joint=joint_test(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist,mlist2)
test_result_marginal1=marginal_test(t, Z, ylist[[1]],mvec, mvec0,mvec1, G, K01, d01, K11, d11,mlist1)
test_result_marginal2=marginal_test(t, Z, ylist[[2]],mvec, mvec0,mvec1, G, K02, d02, K12, d12,mlist1)

test_result_joint2=joint_test2(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist,mlist2)
test_result_marginal12=marginal_test2(t, Z, ylist[[1]],mvec, mvec0,mvec1, G, K01, d01, K11, d11,mlist1)
test_result_marginal22=marginal_test2(t, Z, ylist[[2]],mvec, mvec0,mvec1, G, K02, d02, K12, d12,mlist1)


test_result_joint_interaction=joint_test_interaction(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist,mlist2)
test_result_marginal1_interaction=marginal_test_interaction(t, Z, ylist[[1]],mvec, mvec0,mvec1, G, K01, d01, K11, d11,mlist1)
test_result_marginal2_interaction=marginal_test_interaction(t, Z, ylist[[2]],mvec, mvec0,mvec1, G, K02, d02, K12, d12,mlist1)


test_results=c(test_result_joint$pvalue,test_result_marginal1$pvalue,test_result_marginal2$pvalue,
                        test_result_joint2$pvalue,test_result_marginal12$pvalue,test_result_marginal22$pvalue,
                        test_result_joint_interaction$pvalue,test_result_marginal1_interaction$pvalue,test_result_marginal2_interaction$pvalue)

