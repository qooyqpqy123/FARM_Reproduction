library(tfHuber)
#Main Function: Huber.predict
##Input n: number of observations
#p dimension
#K number of factors
#w: signal strength for non-zero beta
#u: signal strength for \phi
#v: signal strength for factor loading
Huber.predict = function(n,p,K,w,u,v){
  #------------Generate Data ---------------
  F = matrix(rnorm(n*K),n,K) #factor
  U = matrix(rnorm(n*p),n,p) #idiosyncratic component
  B = matrix(runif(p*K,-1,1),p,K)*v #factor loading
  X = F%*%t(B) + U #covariate
  beta.0 = c(rep(w,20),rep(0,p-20)) #true signal beta
	varphi.0 = rep(1,K)*u #true signal \phi
	E = rt(n,df=3)
	Y = F%*%varphi.0 + X%*%beta.0 + E #response variable
	#-------------Use covariance matrix to estimate factors -----------
	SigmaX = tcrossprod(X)/n #covariance matrix
	eigenX = eigen(SigmaX)
	eigvec = eigenX$vectors
	eigvalue = eigenX$values
	K.hat = which.min(diff(log(eigvalue[1:100]))) #estimate number of factors
	F.hat = eigvec[,1:K.hat]*sqrt(n) #estimated factors
	B.hat = t(t(F.hat)%*%X)/n #factor loading
	U.hat = X - F.hat%*%t(B.hat) #idiosyncratic component
##----------Model Estimation-------------------------------------------------------------------
##-----------------------------------------------------------------------------
	fit.lasso.X = cvHuberLasso(X,Y)
	beta.hat.X = as.vector(fit.lasso.X$theta)                                ## Fitted Huber LASSO
	beta.hat.X = beta.hat.X[-1]
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	lmY.F = huberReg(F.hat,Y)
	gamma.hat = as.vector(lmY.F$theta)                                       ## Fitted Huber PCR
	gamma.hat = gamma.hat[-1]
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
	X.hat = cbind(F.hat,U.hat)
	fit.lasso.U = cvHuberLasso(X.hat,Y)
	beta.hat.U = as.vector(fit.lasso.U$theta)                                ## Fitted Huber FARM
	beta.hat.U = beta.hat.U[-1]
##-----------------------------------------------------------------------------
##---------------------Prediction: new sample--------------------------------------------------------
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
##----------------------Model prediction results--------------------------------------------------------
	X.predict = Y - sum(X*beta.hat.X) #prediction of Lasso
	U.predict = Y - sum(X.new*beta.hat.U) #prediction of FARM
	F.predict=  Y - sum(F.new*gamma.hat) #prediction of PCR
	tmp.X = X.predict^2
	tmp.U = U.predict^2
	tmp.F= F.predict^2
	return(c(tmp.F,tmp.X,tmp.U)) #record the prediction error (PCR, LASSO, FARM)
}
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
##-----------------------------------------------------------------------------
#-----------------Main codes -----------------
n_seq = round(seq(200,400,length=8))
w = 0.8; u = 0.8; vv = 2; K = 5; p = 1000
system.time(Huber.predict(400,p,K,w,u,vv))
MEAN1 = matrix(0,length(n_seq),3)
for(i in c(1:length(n_seq))){ #for different n, compute the prediction results of FARM and LASSO
	print(i)  
	n = n_seq[i]  
	A.tmp = replicate(500,Huber.predict(n,p,K,w,u,vv))
	MEAN1[i,] = rowMeans(A.tmp)
	write.csv(MEAN1,'MEAN1.csv') # Table 4 in appendix.
}














