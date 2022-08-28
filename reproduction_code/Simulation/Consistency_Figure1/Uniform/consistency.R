library(glmnet)
#n=200
p=1000
K=2
gamma1=0.5*rep(1,K)
beta=c(rep(0.5,3),rep(0,p-3))
B<-matrix(runif(p*K,-1,1),nrow=p)
#generate data
ratio<-c(0.20,0.25,0.30,0.35,0.40,0.45,0.50)
mat<-matrix(nrow=7,ncol=500)
mat2<-matrix(nrow=7,ncol=500)
i=0
for (r in ratio){
  i=i+1
  dis<-c()
  dis2<-c()
  for (j in 1:500){
n=ceiling(log(p)/(r^2/9))
F<-matrix(rnorm(n*K),nrow=n) #factor
U<-matrix(rnorm(n*p),nrow=n) 
X=F%*%t(B)+U #design matrix
#Y1=X%*%beta+rnorm(n,mean=0,sd=0.5) #Null
Y=F%*%gamma1+U%*%beta+rnorm(n,mean=0,sd=0.5)#runif(n,-sqrt(3)/2,sqrt(3)/2)#rnorm(n,mean=0,sd=0.5) #alternative

cov<-X%*%t(X)/n
eigvec<-eigen(cov)$vectors
eigval<-eigen(cov)$values
K_est<-which.min(diff(log(eigval[1:10]))) #select K
hatf<-sqrt(n)*eigvec[,1:K_est]
hatB<-t(1/n*t(hatf)%*%X)
hatU<-X-hatf%*%t(hatB)
Y_1<-Y-1/n*hatf%*%t(hatf)%*%Y
fit1 = glmnet(hatU, Y_1, intercept=FALSE,
              lambda=cv.glmnet(hatU,Y_1,intercept=FALSE)$lambda.1se)
fit2=glmnet(X, Y,intercept=FALSE,
            lambda=cv.glmnet(hatU,Y_1,intercept=FALSE)$lambda.1se)
lambda_beta1<-fit1$lambda
beta_hat1<-as.vector(fit1$beta)

lambda_beta2<-fit2$lambda
beta_hat2<-as.vector(fit2$beta)


dis<-c(dis,sum(abs(beta_hat1-beta)))
dis2<-c(dis2,sum(abs(beta_hat2-beta)))
write.csv(dis,"dis.csv")
write.csv(dis2,"dis2.csv")
}
mat[i,]<-dis
mat2[i,]<-dis2
write.csv(mat,"mat.csv")
write.csv(mat2,"mat2.csv")
}



