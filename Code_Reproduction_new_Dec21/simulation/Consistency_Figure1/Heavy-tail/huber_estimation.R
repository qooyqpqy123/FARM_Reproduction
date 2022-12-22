#-----------Load Packages ------------#
library(tfHuber)
#Function Ahubert: Generate data and Do model estimation; 
#------------------
#Input: r: ratio, a number between 0-1
#p: dimension
#K number of factors
#theta: true parameter
#s: sparsity level
#----------------
AHubert = function(r,p,K,theta,s){
	gamma0 = 0.5*rep(1,K) #true gamma
	beta0 = c(rep(0.5,s),rep(0,p-s)) #true beta
	B0 = matrix(runif(p*K,-1,1),nrow=p) #factor loading 
##-------------- Generate Data---------------------------
	n = ceiling(log(p)/(r/(K+s))^(1+1/theta))#number of observations
	F = matrix(rnorm(n*K),nrow=n) #generate factors
	U = matrix(rnorm(n*p),nrow=n) #generate idiosyncratic component
	X = F%*%t(B0)+U #generate design matrix
	E = rt(n,df=(2+theta)) #generate noise distribution
	Y = F%*%gamma0+U%*%beta0+E #generate response variable
##--------------- Estimate Factor Model --------------------------
	Sigma = tcrossprod(X)/n #covariance matrix
	Eig = eigen(Sigma,only.values=TRUE) 
	Eigval = Eig$values #eigenvales
	K.est = which.min(diff(log(Eigval[1:20]))) #estimate number of factors
	Svd = svd(X,nu=K.est,nv=0) 
	Eigvec = Svd$u #eigenvectors
	F.hat = sqrt(n)*Eigvec[,1:K.est] #estimate factors
	mol = lm(X~F.hat-1) #fit model
	U.hat = resid(mol) #estimate idiosyncratic component
	X.hat = cbind(F.hat,U.hat) 
##------------Estimate FARM Model------------------------------
	Huber.est = cvHuberLasso(X.hat,Y) #Estimate FARM mode using Huber regression
	Huber.beta = Huber.est$theta 
	Huber.beta = Huber.beta[-(1:(K.est+1))] #estimated beta
	z = sum(abs(Huber.beta-beta0)) #record the difference
	return(z) #return the difference
}

#----------------Main setting---------------------------------
p = 1000; #dimension
K = 2; #number of factors
s = 3; #sparsity level
theta = 1 #parameter
N = 7 #
ratio = seq(0.4,0.7,length=N)
L = matrix(0,N,500)#record the results 500 times
#
for(i in 1:N){
	print(i)
	r = ratio[i] #for every r in ratio, use AHubert 500 times to compute the parameter differences.
	system.time(AHubert(r,p,K,theta,s))
	tmp = replicate(500,AHubert(r,p,K,theta,s))
	tmp = as.numeric(tmp)
	L[i,] = tmp #record
}
write.csv(L,"p1000_K2_s3_t3.csv") #record 500 replication results.
	

	
	
	



















