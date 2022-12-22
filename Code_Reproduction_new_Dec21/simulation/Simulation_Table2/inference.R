#simulation for Test 2
#--------------Load package-------------
library("SIS")
library("MASS")
#-------------Compute the critical value-------------
quant_chisq<-c()
for (q in 1:10){
  boostrap<-c()
  for(j in 1:500){
    boostrap<-c(boostrap,rchisq(1,3))
  }
  quant_chisq<-c(quant_chisq,quantile(boostrap,0.95))
}
c_alpha<-mean(quant_chisq)

#---------------Main function:selectfeature: use first part of data to conduct sure screening and use the second part of data to construct test statistics---------
#input Y1: response
#X: covariate
#n: number of observations
selectfeature<-function(Y1,X,n){
  ###In the following step, we use PCA of sample covariance matrix to estimate latent factors via the first part of data################ 
  cov<-(X[1:floor(n^(0.8)),]-rowMeans(X[1:floor(n^(0.8)),]))%*%t((X[1:floor(n^(0.8)),]-rowMeans(X[1:floor(n^(0.8)),])))/floor(n^(0.8)) #covariance matrix
  eigvec<-eigen(cov)$vectors
  eigval<-eigen(cov)$values
  K1<-which.min(diff(log(eigval[1:30]))) #select K: number of factors
  hatf<-sqrt(floor(n^(0.8)))*eigvec[,1:K1] #estimated factors
  hatB<-t(1/floor(n^(0.8))*t(hatf)%*%X[1:floor(n^(0.8)),]) #estimated factor loadings
  hatU<-X[1:floor(n^(0.8)),]-hatf%*%t(hatB) #estimated idiosyncratic component
Y_new<-Y1[1:floor(n^(0.8))]-1/floor(n^(0.8))*hatf%*%t(hatf)%*%Y1[1:floor(n^(0.8))] #
sis<-SIS(hatU,Y_new) #use sure screening to select variables
select1<-sis$ix #ISIS
###########In the following step, we use screened variables to construct part of test statistics via the second part of data###########
cov2<-(X[(floor(n^(0.8))+1):n,]-rowMeans(X[(floor(n^(0.8))+1):n,]))%*%t(X[(floor(n^(0.8))+1):n,]-rowMeans(X[(floor(n^(0.8))+1):n,]))/(n-(floor(n^(0.8))))
eigvec2<-eigen(cov2)$vectors
eigval2<-eigen(cov2)$values
K<-which.min(diff(log(eigval2[1:30]))) ##select number of factors
hatf2<-sqrt(n-(floor(n^(0.8))))*eigvec2[,1:K] #estimate factors through principal component
hatB2<-t(1/(n-floor(n^(0.8)))*t(hatf2)%*%X[(floor(n^(0.8))+1):n,]) #estimate factor loading
hatU2<-X[(floor(n^(0.8))+1):n,]-hatf2%*%t(hatB2) #estimate idiosyncratic components
P_f<-hatf2%*%t(hatf2)/(n-(floor(n^(0.8)))) #construct projection matrix P_f
P_x<-X[(floor(n^(0.8))+1):n,select1]%*%solve(t(X[(floor(n^(0.8))+1):n,select1])%*%(X[(floor(n^(0.8))+1):n,select1]))%*%t(X[(floor(n^(0.8))+1):n,select1]) #projection matrix Px
P_u<-hatU2[,select1]%*%solve(t(hatU2[,select1])%*%hatU2[,select1])%*%t(hatU2[,select1])#projection matrix Pu
Q1<-t(Y1[(floor(n^(0.8))+1):n])%*%(P_f+P_u-P_x)%*%(P_f+P_u-P_x)%*%Y1[(floor(n^(0.8))+1):n,] #test statistics
c(Q1) #Output: test statistics
}
#---------------helper function rcv: return the estimated variance of Y given X------------
#Input  Y1: variance
#X: covariate
#n: number of observations
rcv<-function(Y1,X,n){
  Idx1=1:(n/2)
  Idx2=(n/2+1):n
  lis=list(Idx1,Idx2)
  var_ISIS<-c()
  var_SIS<-c()
  for (Idx in lis){
  cov<-2*(X[Idx,]-rowMeans(X[Idx,]))%*%t((X[Idx,]-rowMeans(X[Idx,])))/n #construct factors
  eigvec<-eigen(cov)$vectors
  eigval<-eigen(cov)$values
  K1<-which.min(diff(log(eigval[1:100]))) #select K
  hatf<-sqrt(n/2)*eigvec[,1:K1] #hatF
  hatB<-t(2/n*t(hatf)%*%X[Idx,]) #hat B
  hatU<-X[Idx,]-hatf%*%t(hatB) #hat U
  Y_new<-Y1[Idx]-2/n*hatf%*%t(hatf)%*%Y1[Idx] #
  sis<-SIS(hatU,Y_new)
  select1<-sis$ix #ISIS  
  select2<-sis$ix0 #SIS
  cov2<-(X[-Idx,]-rowMeans(X[-Idx,]))%*%t(X[-Idx,]-rowMeans(X[-Idx,]))/(n/2)
  eigvec2<-eigen(cov2)$vectors
  eigval2<-eigen(cov2)$values
  K<-which.min(diff(log(eigval2[1:100]))) #select K
  hatf2<-sqrt(n/2)*eigvec2[,1:K] #hat F using second half of data
  hatB2<-t(2/n*t(hatf2)%*%X[-Idx,]) #hat B
  hatU2<-X[-Idx,]-hatf2%*%t(hatB2) #hat U  
X0=cbind(hatf2,hatU2)
X0.ISIS = X0[,c(1:K,select1+K)]
X0.SIS = X0[,c(1:K,select2+K)]
lm0.ISIS = lm(Y1[-Idx]~X0.ISIS-1) 
lm0.SIS = lm(Y1[-Idx]~X0.SIS-1)
Q0.ISIS.H1 = sum((resid(lm0.ISIS))^2)
sigma.hat0.ISIS = Q0.ISIS.H1/((n/2) - K - length(select1))
Q0.SIS.H1 = sum((resid(lm0.SIS))^2)
sigma.hat0.SIS = Q0.SIS.H1/((n/2) - K - length(select2))
#X2=cbind(F.hat2,U.hat2)
var_ISIS<-c(var_ISIS,sigma.hat0.ISIS)
var_SIS<-c(var_SIS,sigma.hat0.SIS)
  }
c(mean(var_ISIS),mean(var_SIS)) #return the estimated variance of Y, using (ISIS method, SIS method), respectively.
}


#------------- Main codes -----------------------
choose_p=2 #we can let choose_p be 1 or 2; if it is 1, then $p=250$, otherwise it $p=600$.
choose_noise=2 #we can let choose_noise be 1 or 2; if it is 1, then noise follows Gaussian, otherwise it follows Uniform distribution.
choose_mix=2 #we can let choose_mix be 1 or 2; if it is 1, then covariate are i.i.d, otherwise it follows strong mixing.
n=250 #observations
if (choose_p==1){
  p=250}else{
    p=600
  }
K=3 #number of factors
Phi<-matrix(nrow=K,ncol=K)
for (i in c(1:K)){
  for (j in c(1:K)){
    Phi[i,j]<-0.5^(abs(i-j)+1)
  }
}
F_start=rnorm(K)
F_origin=F_start
for (t in c(1:200)){
  F_start=Phi%*%F_origin+rnorm(K)
  F_origin=F_start
}

top_mat<-matrix(nrow=p,ncol=p)
for (i in c(1:p)){
  for (j in c(1:p)){
    if(i==j){ top_mat[i,j]=1}
    else
    {top_mat[i,j]=0}
  }
}
beta=c(rep(0.9,4),rep(0,p-4)) #true beta
sigma=0.5
B<-matrix(runif(p*K,-1,1),nrow=p) #factor loading

#-----------Record the size or power of test when $m$ takes 6 values------------#
mat1<-matrix(nrow=4,ncol=6)
mat2<-matrix(nrow=4,ncol=6)
mat3<-matrix(nrow=4,ncol=6)
mat4<-matrix(nrow=4,ncol=6)
mat5<-matrix(nrow=4,ncol=6)
mat6<-matrix(nrow=4,ncol=6)
mat_isis<-matrix(nrow=4,ncol=6)
mat_sis<-matrix(nrow=4,ncol=6)
for (m in c(0,0.4,0.8,1.2,1.6,2.0)){
  idx<-c() #record p-values for every replication
  idx2<-c()
  idx3<-c()
  idx4<-c()
  idx5<-c()
  idx6<-c()
  mean_isis<-c()
  mean_sis<-c()
  gamma1=c(rep(m/10,K))
  for(f in c(1:4)){
    isis1<-c()
    isis2<-c()
    isis3<-c()
    isis5<-c()
    isis4<-c()
    isis6<-c()
    sigma_isis<-c()
    sigma_sis<-c()
    ##select features and do inference###
    for (i in 1:1000){ #1000 replications
      if (choose_mix==1){ #Generate factors; 1:
        F<-matrix(nrow=n,ncol=K) #factor
        for (r in c(1:n)){
          F_start=rnorm(K)
          F_origin=F_start
          F[r,]<-F_start
        }
        U<-mvrnorm(n, rep(0,p),top_mat) #noise
      }else{
      F<-matrix(nrow=n,ncol=K) #factor
      for (r in c(1:n)){
        F_start=Phi%*%F_origin+rnorm(K)
        F_origin=F_start
        F[r,]<-F_start
      }
      U<-mvrnorm(n, rep(0,p),top_mat) #noise
      }
      X=F%*%t(B)+U #design matrix
      if (choose_noise==1){ #generate noise distribution 1:Gaussian; 2: Uniform
      Y2=F%*%gamma1+X%*%beta+rnorm(n,mean=0,sd=0.5)}else{
        Y2=F%*%gamma1+X%*%beta+runif(n,-sqrt(3)/2,sqrt(3)/2)
      }
      output1<-try(selectfeature(Y2,X,n),silent=TRUE) #estimate the test statistics
      if ('try-error' %in% class(output1)) {
        output1<-c(0)}
      sig_est<-try(rcv(Y2,X,n),silent=TRUE) #estimate variance through function rcv
      if ('try-error' %in% class(sig_est)) {
        sig_est<-c(0.25,0.25)}
      var_est1<-sig_est[1] #estimated variance via isis.
      sigma_isis<-c(sigma_isis,var_est1) #keep the variance
      isis1<-c(isis1,output1[1]/var_est1)#sigma_isis)
    }
    idx<-c(idx,sum(isis1>c_alpha)/1000) #record p-values
    write.csv(idx,"idx.csv")
  }
  mat1[,(floor(m/0.4)+1)]<-idx #record p_values for 1000 replications
  write.csv(mat1,"mat1.csv")
}
