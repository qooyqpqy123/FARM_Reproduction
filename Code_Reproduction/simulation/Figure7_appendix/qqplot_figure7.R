#############simulation for Test 2 Using No-splitted Data#############
#######################################################################
set.seed(100)
#-------------Load Packages------------
library("SIS")
library("MASS")
#------------Generate Random Numbers from \chisquare distribution with freedom 3 ---------
quant_chisq<-c()
for (q in 1:10){
  boostrap<-c()
  for(j in 1:500){
    boostrap<-c(boostrap,rchisq(1,3))
  }
  quant_chisq<-c(quant_chisq,quantile(boostrap,0.95))
}
c_alpha<-mean(quant_chisq) #compute the critical value

#--------Main function selectfeature: construct test statistics---------------
#Input Y1: response; X: covariate, n: number of observations
selectfeature<-function(Y1,X,n){
  #-------Use covariance matrix to estimate factors -------------------
  cov<-(X[1:n,]-rowMeans(X[1:n,]))%*%t((X[1:n,]-rowMeans(X[1:n,])))/floor(n) #construct covariance matrix
  eigvec<-eigen(cov)$vectors #eigenvectors
  eigval<-eigen(cov)$values #eigenvalues
  K1<-which.min(diff(log(eigval[1:30]))) #select number of factors
  hatf<-sqrt(n)*eigvec[,1:K1] #estimated factors
  hatB<-t(1/n*t(hatf)%*%X[1:n,]) #estimated factor loading
  hatU<-X[1:n,]-hatf%*%t(hatB)  #estimated idiosyncratic component
  Y_new<-Y1[1:n]-1/floor(n)*hatf%*%t(hatf)%*%Y1[1:n] 
  sis<-SIS(hatU,Y_new) #use sure screening to select variable
  #S.hat0.ISIS = SIS0$ix; S.hat0.SIS = SIS0$ix0
  #-----------use covariance matrix to construct test statistics---------
  cov2<-(X[1:n,]-rowMeans(X[1:n,]))%*%t(X[1:n,]-rowMeans(X[1:n,]))/(n) # covariance matrix
  eigvec2<-eigen(cov2)$vectors #eigevectors
  eigval2<-eigen(cov2)$values  #eigenvalues
  K<-which.min(diff(log(eigval2[1:30]))) #select the number of factors
  hatf2<-sqrt(n)*eigvec2[,1:K] #estimated factors
  hatB2<-t(1/n*t(hatf2)%*%X[1:n,]) #factor loading
  hatU2<-X[1:n,]-hatf2%*%t(hatB2) #idiosyncratic component
  P_f<-hatf2%*%t(hatf2)/(n) #construct projection matrix Pf
  select2<-sis$ix0 #SIS selected variable
  P_x<-X[1:n,select2]%*%solve(t(X[1:n,select2])%*%(X[1:n,select2]))%*%t(X[1:n,select2]) #projection matrix Px
  P_u<-hatU2[,select2]%*%solve(t(hatU2[,select2])%*%hatU2[,select2])%*%t(hatU2[,select2]) #projection matrix Pu
  Q2<-t(Y1[1:n])%*%(P_f+P_u-P_x)%*%(P_f+P_u-P_x)%*%Y1[1:n,] #test statistics
  #X2=cbind(F.hat2,U.hat2)
  c(Q2) #output test statistics in section 4.
}

#-------------Helper Function rcv: compute the conditional variance of Y1 given X---------------------
#Input Y1: response; X: covariate; n: number of observations
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
    K1<-which.min(diff(log(eigval[1:100]))) #number of factors
    hatf<-sqrt(n/2)*eigvec[,1:K1] #estimated factors
    hatB<-t(2/n*t(hatf)%*%X[Idx,]) #factor loading
    hatU<-X[Idx,]-hatf%*%t(hatB) #idiosyncratic component
    Y_new<-Y1[Idx]-2/n*hatf%*%t(hatf)%*%Y1[Idx] 
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
    var_ISIS<-c(var_ISIS,sigma.hat0.ISIS)
    var_SIS<-c(var_SIS,sigma.hat0.SIS)
  }
  print(var_ISIS)
  print(var_SIS)
  c(mean(var_ISIS),mean(var_SIS)) #output estimate variance using ISIS method and SIS method, respectively.
}
################The main procedure####################
######################################################
#-----------Data Generation Setting ---------------
n=250 #number of observation
p=250 #dimension
K=3 #number of factors
#---------Basic setting -----------
top_mat<-matrix(nrow=p,ncol=p)
for (i in c(1:p)){
  for (j in c(1:p)){
    if(i==j){ top_mat[i,j]=1}
    else
    {top_mat[i,j]=0}
  }
}
beta=c(rep(0.9,4),rep(0,p-4)) #true beta
sigma=0.5 #noise level
B<-matrix(runif(p*K,-1,1),nrow=p) #factor loading



#-------------inference procedure----------------
#X: design matrix, select1: selected indexes from SIS/ISIS,Y1: responses



mat4<-matrix(nrow=4,ncol=6)
mat_isis<-matrix(nrow=4,ncol=6)
mat_sis<-matrix(nrow=4,ncol=6)
for (m in c(0,0.6,1.2,1.8,2.4,3)){
  idx4<-c()
  mean_isis<-c()
  mean_sis<-c()
  gamma1=c(rep(m/10,K))
  for(f in c(1:4)){
    isis4<-c()
    sigma_isis<-c()
    sigma_sis<-c()
    ##select features and do inference, replicate 1000 times###
    for (i in 1:1000){
      F<-matrix(nrow=n,ncol=K) #generate factor
      for (r in c(1:n)){
        F_start=rnorm(K)
        F_origin=F_start
        F[r,]<-F_start
      }
      U<-mvrnorm(n, rep(0,p),top_mat)  #generate idiosyncratic component
      X=F%*%t(B)+U                          #design matrix
      Y2=F%*%gamma1+X%*%beta+rnorm(n,mean=0,sd=0.5) #generate response
      output1<-try(selectfeature(Y2,X,n),silent=TRUE)
      if ('try-error' %in% class(output1)) { #estimate test statistics
        output1<-c(0)}
      sig_est<-try(rcv(Y2,X,n),silent=TRUE)
      if ('try-error' %in% class(sig_est)) { #estimate variance
        sig_est<-c(0.25,0.25)}
      var_est1<-sig_est[1]
      isis4<-c(isis4,output1/var_est1) #normalized test statistics
    }
    if ((f==1)&(m==0)){
      write.csv(isis4,"sure.csv")  #record the data when null hypothesis holds.
    }
    mean_isis<-c(mean_isis,mean(sigma_isis))  #record variance
    mean_sis<-c(mean_sis,mean(sigma_sis))
    idx4<-c(idx4,sum(isis4>c_alpha)/1000) #record the result of the thest
  }
  mat4[,(floor(m/0.6)+1)]<-idx4
  mat_isis[,(floor(m/0.6)+1)]<-mean_isis
  mat_sis[,(floor(m/0.6)+1)]<-mean_sis
}


############## the QQ-Plot###############
no_split<-read.csv("sure.csv") #read sure.csv
quantile_null<-no_split[,2] 
length(quantile_null)

quantile_true<-rchisq(1000,3) #benchmark

qq<-qqplot(quantile_null, quantile_true,xlab="",ylab="") #plot qq-plot
lines(c(0,2000),c(0,2000),col='red')