library(tfHuber)
AHubert = function(r,p,K,theta,s){
	gamma0 = 0.5*rep(1,K)
	beta0 = c(rep(0.5,s),rep(0,p-s))
	B0 = matrix(runif(p*K,-1,1),nrow=p)
##-----------------------------------------
	n = ceiling(log(p)/(r/(K+s))^(1+1/theta))
	F = matrix(rnorm(n*K),nrow=n)
	U = matrix(rnorm(n*p),nrow=n)
	X = F%*%t(B0)+U
	E = rt(n,df=(2+theta))
	Y = F%*%gamma0+U%*%beta0+E
##-----------------------------------------
	Sigma = tcrossprod(X)/n
	Eig = eigen(Sigma,only.values=TRUE)
	Eigval = Eig$values
	K.est = which.min(diff(log(Eigval[1:20])))
	Svd = svd(X,nu=K.est,nv=0)
	Eigvec = Svd$u
	F.hat = sqrt(n)*Eigvec[,1:K.est]
	mol = lm(X~F.hat-1)
	U.hat = resid(mol)
	X.hat = cbind(F.hat,U.hat)
##------------------------------------------
	Huber.est = cvHuberLasso(X.hat,Y)
	Huber.beta = Huber.est$theta
	Huber.beta = Huber.beta[-(1:(K.est+1))]
	z = sum(abs(Huber.beta-beta0))
	return(z)
}

p = 1000; K = 2; s = 3; theta = 1
N = 7
ratio = seq(0.4,0.7,length=N)
L = matrix(0,N,500)
for(i in 1:N){
	print(i)
	r = ratio[i]
	system.time(AHubert(r,p,K,theta,s))
	tmp = replicate(500,AHubert(r,p,K,theta,s))
	tmp = as.numeric(tmp)
	L[i,] = tmp
}
meanL = rowMeans(L)
plot(ratio,meanL,type='b')
abline(0,1)
write.csv(L,"p1000_K2_s3_t3.csv")
	

	
	
	



















