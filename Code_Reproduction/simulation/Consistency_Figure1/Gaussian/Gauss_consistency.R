#-----------Loading Packages ------------#
library(glmnet)
#----------Basic Setting to generate data----------#
p=1000   #dimension
K=2 #true number of factors
gamma1=0.5*rep(1,K) #generate true gamma in FARM
beta=c(rep(0.5,3),rep(0,p-3)) #generate true beta
B<-matrix(runif(p*K,-1,1),nrow=p) #generate factor loading
ratio<-c(0.20,0.25,0.30,0.35,0.40,0.45,0.50) #the x-axis of figure 1.
mat<-matrix(nrow=7,ncol=500) #record results for FARM
mat2<-matrix(nrow=7,ncol=500) #record results for Lasso
i=0
for (r in ratio){
  i=i+1
  dis<-c()
  dis2<-c()
  for (j in 1:500){ #500 replications
n=ceiling(log(p)/(r^2/9))
F<-matrix(rnorm(n*K),nrow=n) #factor
U<-matrix(rnorm(n*p),nrow=n) #idiosyncratic component
X=F%*%t(B)+U                  #design matrix
Y=F%*%gamma1+U%*%beta+rnorm(n,mean=0,sd=0.5) #Data is generated following Gaussian noise 
#-----------------Use eigevector of covariance matrix to estimate latent factors -------------#
cov<-X%*%t(X)/n
eigvec<-eigen(cov)$vectors #eigenvectors
eigval<-eigen(cov)$values  #eigenvalues
K_est<-which.min(diff(log(eigval[1:10]))) #estimate factor numbers using eigenvalue ratios.
hatf<-sqrt(n)*eigvec[,1:K_est] #estimated factors
hatB<-t(1/n*t(hatf)%*%X) #Estimated Factor Loading
hatU<-X-hatf%*%t(hatB)   #Estimated Idiosyncratic component
Y_1<-Y-1/n*hatf%*%t(hatf)%*%Y
fit1 = glmnet(hatU, Y_1, intercept=FALSE,
              lambda=cv.glmnet(hatU,Y_1,intercept=FALSE)$lambda.1se) #Estimate FARM model
fit2=glmnet(X, Y,intercept=FALSE,
            lambda=cv.glmnet(hatU,Y_1,intercept=FALSE)$lambda.1se)  #Estimate Lasso model
lambda_beta1<-fit1$lambda    #fitted beta for FARM
beta_hat1<-as.vector(fit1$beta)

lambda_beta2<-fit2$lambda   #fitted beta for Lasso
beta_hat2<-as.vector(fit2$beta)


dis<-c(dis,sum(abs(beta_hat1-beta))) #record the difference between estimated FARM beta with true one
dis2<-c(dis2,sum(abs(beta_hat2-beta))) ##record the difference between estimated LASSO beta with true one
write.csv(dis,"dis.csv") #record FARM
write.csv(dis2,"dis2.csv") #record FARM
}
mat[i,]<-dis    #Repeat 500 times and record as a matrix.
mat2[i,]<-dis2
write.csv(mat,"mat.csv")
write.csv(mat2,"mat2.csv")
}


