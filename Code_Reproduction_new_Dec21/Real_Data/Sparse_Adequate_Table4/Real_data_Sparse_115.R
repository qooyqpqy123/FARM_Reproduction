library(SIS)
library(MASS)
#input Data############################
set.seed(100)
real_data<-read.csv("realdata_115_126.csv")
#########Preprocessing #################
X<-data.frame(real_data)
X=X[,-1]        #Delete the index column
X<-as.matrix(X)
X<-X[,-c(81,83)] #Delete high-correlated columns
n=nrow(X)        #number of rows
p=ncol(X)-1      #number of columns
X = t(t(X)-colMeans(X))  #normalize X
X = t(t(X)/apply(X,2,sd))
#########################################

#---------------------Main Function: selectfeature--------------------------
#In this main function, we first use a small size of data to conduct sure screening
#Then we construct test statistics based on screened features
#------------------------------------------------------------------------------
selectfeature<-function(Y1,X,n){ #input: Response Y1, covariate X, number of samples n.
   ####In the following step, we use PCA of sample covariance matrix to estimate latent factors ################ 
  cov<-(X[1:floor(n^(0.8)),]-rowMeans(X[1:floor(n^(0.8)),]))%*%t((X[1:floor(n^(0.8)),]-rowMeans(X[1:floor(n^(0.8)),])))/floor(n^(0.8)) #construct covariance matrix
  eigvec<-eigen(cov)$vectors                      #compute the eigenvectors of covariance
  eigval<-eigen(cov)$values                      #compute the eigenvalues of covariance
  K1<-max(2,which.min(diff(log(eigval[1:10]))))  #select the number of factors
  hatf<-sqrt(floor(n^(0.8)))*eigvec[,1:K1]       #estimate the factors using eigenvectors
  hatB<-t(1/floor(n^(0.8))*t(hatf)%*%X[1:floor(n^(0.8)),]) #Estimate the Factor Loadings
  hatU<-X[1:floor(n^(0.8)),]-hatf%*%t(hatB)      #Estimate the Idiosyncratic component
  ##########In the following step, we use iterative sure screening to select Features #################
  Y_new<-Y1[1:floor(n^(0.8))]-1/floor(n^(0.8))*hatf%*%t(hatf)%*%Y1[1:floor(n^(0.8))] #We first remove the factor effect from the response
  sis<-SIS(hatU,Y_new) #Use sure screening to screen the variables
  select2<-sis$ix0 #Selected variables
  ###########In the following step, we use screened variables to construct part of test statistics ###########
  cov2<-(X[(floor(n^(0.8))+1):n,]-rowMeans(X[(floor(n^(0.8))+1):n,]))%*%t(X[(floor(n^(0.8))+1):n,]-rowMeans(X[(floor(n^(0.8))+1):n,]))/(n-(floor(n^(0.8))))
  eigvec2<-eigen(cov2)$vectors
  eigval2<-eigen(cov2)$values
  K<-max(2,which.min(diff(log(eigval2[1:10]))))            #select K
  hatf2<-sqrt(n-(floor(n^(0.8))))*eigvec2[,1:K]            #hat F using second half of data
  hatB2<-t(1/(n-floor(n^(0.8)))*t(hatf2)%*%X[(floor(n^(0.8))+1):n,]) #Factor loading
  hatU2<-X[(floor(n^(0.8))+1):n,]-hatf2%*%t(hatB2)          #Idiosyncratic Component
  P_f<-hatf2%*%t(hatf2)/(n-(floor(n^(0.8))))
  P_x<-X[(floor(n^(0.8))+1):n,select2]%*%solve(t(X[(floor(n^(0.8))+1):n,select2])%*%(X[(floor(n^(0.8))+1):n,select2]))%*%t(X[(floor(n^(0.8))+1):n,select2])
  P_u<-hatU2[,select2]%*%solve(t(hatU2[,select2])%*%hatU2[,select2])%*%t(hatU2[,select2])
  Q2<-t(Y1[(floor(n^(0.8))+1):n])%*%(P_f+P_u-P_x)%*%(P_f+P_u-P_x)%*%Y1[(floor(n^(0.8))+1):n] #Test Statistics
  c(Q2,K)                          #Output: Test Statistics, number of selected factors
}
#---------------------------------------------------------------------------------------
###Introduce rcv: Helper function for estimating the variance $\hat\sigma$ using refitted-validation.
#----------------------------------------------------------------------------------------
rcv<-function(Y1,X,n){###Input Response Y, covariate X, sample size $n$. 
  Idx1=1:(n/2)
  Idx2=(n/2+1):n
  lis=list(Idx1,Idx2)
  var_ISIS<-c()
  var_SIS<-c()
  for (Idx in lis){
    cov<-2*(X[Idx,]-rowMeans(X[Idx,]))%*%t((X[Idx,]-rowMeans(X[Idx,])))/n #construct factors
    eigvec<-eigen(cov)$vectors
    eigval<-eigen(cov)$values
    K1<-which.min(diff(log(eigval[1:10])))
    hatf<-sqrt(n/2)*eigvec[,1:K1] #Estimated Factods using first half of data
    hatB<-t(2/n*t(hatf)%*%X[Idx,]) #Factor loading
    hatU<-X[Idx,]-hatf%*%t(hatB) #Idiosyncratic Component
    Y_new<-Y1[Idx]-2/n*hatf%*%t(hatf)%*%Y1[Idx] 
    sis<-SIS(hatU,Y_new)
    select1<-sis$ix 
    select2<-sis$ix0 
    cov2<-(X[-Idx,]-rowMeans(X[-Idx,]))%*%t(X[-Idx,]-rowMeans(X[-Idx,]))/(n/2)
    eigvec2<-eigen(cov2)$vectors
    eigval2<-eigen(cov2)$values
    K<-which.min(diff(log(eigval2[1:10]))) 
    hatf2<-sqrt(n/2)*eigvec2[,1:K]    #Estimated Factods using second half of data
    hatB2<-t(2/n*t(hatf2)%*%X[-Idx,]) #Factor loading
    hatU2<-X[-Idx,]-hatf2%*%t(hatB2)  #Idiosyncratic Component
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
  c(mean(var_ISIS),mean(var_SIS)) ##Output estimated variances, the first entry is using ISIS method and the second entry is using SIS method.
}


p_value<-c() #record p-values
##############Main Codes##################
i=49 #Choose i=49 represent HOUSENE, choose i=81 to be GS5
Y1=X[,i] #response
X1=X[,-i] #covariate
output1<-try(selectfeature(Y1,X1,n),silent=TRUE)
if ('try-error' %in% class(output1)) {
  output1<-c(0,1)} #Use main function to compute test statistics
sig_est<-try(rcv(Y1,X1,n),silent=TRUE) #estimate variance using helper function rcv
var_est1<-sig_est[1] #variance
est<-output1[1]/var_est1 #Values of test statistics
p_value<-c(p_value,1-pchisq(est,output1[2])) #p-values of the test
print(p_value) #output p_value

