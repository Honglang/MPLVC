# MPLVC: Multivariate Partial Linear Varying Coefficients Model for Gene-Environmental Interactions with Multiple Longitudinal Traits

This code implements the genetic association studeis via QIF with multiple longitudinal traits using partial linear VC models.

## Example

* Prepare your data into the format: multivariate traits in a list, variable of the varying coefficients, genetic variant, and other covariates

* Choose the model choice: number of knots as well as degree for B-spline; working correlation for QIF; GCV to select the penalty parameter for overfitting

* Apply our proposed method for modeling estimation and testing


```r
#########################################################################
####################### Generate Demo Data  #############################
#########################################################################
source("./mplvc_functions.R")

N=200
m=5
pr=0.3
nul=0.001
nul2=0

data=obs_sim(N,m,pr,rho_vec=c(0.5,0.5,0.2),tau=0.2,sigma_vec=c(sqrt(0.1),sqrt(0.1)),alpha_vec=c(0.2,0.3),nul,nul2)

ylist=data$ylist #multivariate triat data
t=data$t #varying in t
G=data$G #genetic variant
Z=data$Z #other covariates

mvec=data$mvec # number of repeated measurements for each subject
mvec0=data$mvec0 # starting position for the data block for each subject
mvec1=data$mvec1 # ending position for the data block for each subject



###################################################
####  Model selection    ##########################
###################################################

#choose the working correlation such as AR(1) below
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




# choose the number of knots as well as the order for Bspline
Klist=vector("list",2)
dlist=vector("list",2)
Klist[[1]]=rep(floor(0.6*N^(1/5)),2)  
Klist[[2]]=rep(floor(0.6*N^(1/5)),2)  
dlist[[1]]=rep(df,2)
dlist[[2]]=rep(df,2)
K01=Klist[[1]][1];K11=Klist[[1]][2];K02=Klist[[2]][1];K12=Klist[[2]][2]
d01=dlist[[1]][1];d11=dlist[[1]][2];d02=dlist[[2]][1];d12=dlist[[2]][2]


#Using GCV to choose the penalty parameter for overfitting
GCV = function(lam){estimation_full(t, Z, ylist,mvec, mvec0,mvec1,G, Klist,dlist, mlist2,lam)$GCV}
(lamb_est = golden(GCV, 0.01, 2))


###################################################
#### Estimation and Testing    ####################
###################################################

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




```