#-------Load Package ---------
library(glmnet)
#-------Basic Setting ---------
p=1000 #dimension
K=2 #number of factors
N=7 #points on x-axis 
gamma1=0.5*rep(1,K) #true gamma
beta=c(rep(0.5,3),rep(0,p-3)) #true beta
B<-matrix(runif(p*K,-1,1),nrow=p) #factor loading
ratio = seq(0.4,0.7,length=N)
mat<-matrix(nrow=7,ncol=500) #record 500 times
mat2<-matrix(nrow=7,ncol=500)
i=0
for (r in ratio){
  i=i+1
  dis<-c()
  dis2<-c()
  for (j in 1:500){
#-----------------Generate Data -------------
n=ceiling(log(p)/(r^2/9))
F<-matrix(rnorm(n*K),nrow=n) #factors
U<-matrix(rnorm(n*p),nrow=n) #idiosyncratic component
X=F%*%t(B)+U                  #design matrix
Y=F%*%gamma1+U%*%beta+rt(n, 3) #response, the noise is t_3 distributed
#----------------Use covariance matrix to estimate factors--------
cov<-X%*%t(X)/n   #covariance matrix
eigvec<-eigen(cov)$vectors #eigenvector
eigval<-eigen(cov)$values #eigenvalue
K_est<-which.min(diff(log(eigval[1:10]))) #select number of factors
hatf<-sqrt(n)*eigvec[,1:K_est] #estimated factors
hatB<-t(1/n*t(hatf)%*%X) #factor loading
hatU<-X-hatf%*%t(hatB) #Idiosyncratic component
Y_1<-Y-1/n*hatf%*%t(hatf)%*%Y 
fit1 = glmnet(hatU, Y_1, intercept=FALSE,
              lambda=cv.glmnet(hatU,Y_1,intercept=FALSE)$lambda.1se) #Fit FARM model
lambda_beta1<-fit1$lambda #fitted beta
beta_hat1<-as.vector(fit1$beta)
dis<-c(dis,sum(abs(beta_hat1-beta))) #record the difference between fitted value and true value
write.csv(dis,"dis.csv")
}
mat[i,]<-dis #record 500 times.
write.csv(mat,"mat.csv")
}



