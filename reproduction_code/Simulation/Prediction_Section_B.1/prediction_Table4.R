library(tfHuber)
##zzzzz = rnorm(77777)
Huber.predict = function(n,p,K,w,u,v){
	F = matrix(rnorm(n*K),n,K)
	U = matrix(rnorm(n*p),n,p)
	B = matrix(runif(p*K,-1,1),p,K)*v
	X = F%*%t(B) + U
	beta.0 = c(rep(w,20),rep(0,p-20))
	varphi.0 = rep(1,K)*u
	E = rt(n,df=3)
	Y = F%*%varphi.0 + X%*%beta.0 + E
	SigmaX = tcrossprod(X)/n
	eigenX = eigen(SigmaX)
	eigvec = eigenX$vectors
	eigvalue = eigenX$values
	K.hat = which.min(diff(log(eigvalue[1:100])))
	F.hat = eigvec[,1:K.hat]*sqrt(n)
	B.hat = t(t(F.hat)%*%X)/n
	U.hat = X - F.hat%*%t(B.hat)
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	fit.lasso.X = cvHuberLasso(X,Y)
	beta.hat.X = as.vector(fit.lasso.X$theta)                                ## Huber LASSO
	beta.hat.X = beta.hat.X[-1]
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	lmY.F = huberReg(F.hat,Y)
	gamma.hat = as.vector(lmY.F$theta)                                       ## Huber PCR
	gamma.hat = gamma.hat[-1]
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	X.hat = cbind(F.hat,U.hat)
	fit.lasso.U = cvHuberLasso(X.hat,Y)
	beta.hat.U = as.vector(fit.lasso.U$theta)                                ## Huber F-LASSO
	beta.hat.U = beta.hat.U[-1]
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	F = rnorm(K)
	U = rnorm(p)
	X = B%*%F + U                                                            
	E = rt(1,df=3)
	Y = sum(F*varphi.0) + sum(X*beta.0) + E
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	lmX.B = lm(X~B.hat-1)
	F.new = coef(lmX.B)
	U.new = resid(lmX.B)
	X.new = c(F.new,U.new)
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	X.predict = Y - sum(X*beta.hat.X)
	U.predict = Y - sum(X.new*beta.hat.U)
	F.predict=  Y - sum(F.new*gamma.hat)
	tmp.X = X.predict^2
	tmp.U = U.predict^2
	tmp.F= F.predict^2
	return(c(tmp.F,tmp.X,tmp.U))
}
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
n_seq = round(seq(200,400,length=8))
w = 0.8; u = 0.8; vv = 2; K = 5; p = 1000
system.time(Huber.predict(400,p,K,w,u,vv))
MEAN1 = matrix(0,length(n_seq),3)
SD1 = matrix(0,length(n_seq),3)
for(i in c(1:length(n_seq))){
	print(i)  
	n = n_seq[i]  
	A.tmp = replicate(500,Huber.predict(n,p,K,w,u,vv))
	MEAN1[i,] = rowMeans(A.tmp)
	SD1[i,]= apply(A.tmp,1,sd)
	write.csv(MEAN1,'MEAN1.csv')
	write.csv(SD1,'SD1.csv')
}














