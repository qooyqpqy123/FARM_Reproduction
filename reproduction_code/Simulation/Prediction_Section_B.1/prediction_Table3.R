library(glmnet)
PCR.predict = function(n,p,K,w,u,v){
  F = matrix(rnorm(n*K),n,K)
  U = matrix(rnorm(n*p),n,p)
  B = matrix(runif(p*K,-1,1),p,K)*v
  X = F%*%t(B) + U
  beta.0 = c(rep(w,20),rep(0,p-20))
  varphi.0 = rep(1,K)*u
  sigma.0 = 0.5
  E = rnorm(n)*sigma.0
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
  cv.fit.X = cv.glmnet(X,Y,intercept=FALSE)	
  lambda.fit.X = cv.fit.X$lambda.min
  fit.X = glmnet(X,Y,intercept=FALSE,lambda=lambda.fit.X)
  beta.hat.X = as.vector(fit.X$beta)
  ##----------------------------------------------------------------------------
  lmY.F = lm(Y~F.hat-1)
  gamma.hat = coef(lmY.F)
  Y.tilde = resid(lmY.F)
  cv.fit.U = cv.glmnet(U.hat,Y.tilde,intercept=FALSE)	
  lambda.fit.U = cv.fit.U$lambda.min
  fit.U = glmnet(U.hat,Y.tilde,intercept=FALSE,lambda=lambda.fit.U)
  beta.hat.U = as.vector(fit.U$beta)
  ##----------------------------------------------------------------------------
  F = rnorm(K)
  U = rnorm(p)
  X = B%*%F + U
  E = rnorm(1)*sigma.0
  Y = sum(F*varphi.0) + sum(X*beta.0) + E
  ##----------------------------------------------------------------------------
  lmX.B = lm(X~B.hat-1)
  F.new = coef(lmX.B)
  U.new = resid(lmX.B)
  ##----------------------------------------------------------------------------
  X.predict = Y - sum(X*beta.hat.X)
  U.predict = Y - sum(F.new*gamma.hat) - sum(U.new*beta.hat.U)
  F.predict=  Y - sum(F.new*gamma.hat)
  tmp.X = X.predict^2
  tmp.U = U.predict^2
  tmp.F= F.predict^2
  return(c(tmp.F,tmp.X,tmp.U))
}
#n = 100; p = 200; K = 3; w = 3; u = 0
#system.time(PCR.predict(n,p,K,w,u,1))
#vv = seq(1,2,0.2)
#m = length(vv)
#VV = matrix(0,m,2)
#for(i in 1:m){
#	v.tmp = vv[i]
##	A.tmp = replicate(2000,PCR.predict(n,p,K,w,u,v.tmp))
#	VV[i,] = rowMeans(A.tmp)
#}
#vv1 = VV
#vv1
n_seq=round(seq(200,400,length=8))
K = 5
w = 0.8
vv = 2
p=1000
u=0.8
#vv_seq=c(0.2,0.4,0.6,0.8,1.0,1.2,1.4)
#seq_p=c(50,100,200,400,800,1600,3200)
MEAN1=matrix(0,length(n_seq),3)
SD1=matrix(0,length(n_seq),3)
for (i in c(1:length(n_seq))){
  print(i)
  #n = 100; p = 1000; K = 3; w = 3; u = 2
  #system.time(PCR.predict(n,p,K,w,u,1))
  #u=vv_seq[i]
  n=n_seq[i]
  #VV = matrix(0,1,2)
  A.tmp = replicate(500,PCR.predict(n,p,K,w,u,vv))
  MEAN1[i,] = rowMeans(A.tmp)
  SD1[i,]=apply(A.tmp,1,sd)
  #vv2 = VV
  write.csv(MEAN1,'~/Desktop/MEAN.csv')
  write.csv(SD1,'~/Desktop/SD.csv')
}

#n = 100; p = 200; K = 3; w = 3; u = 4
#system.time(PCR.predict(n,p,K,w,u,1))
#vv = seq(1,2,0.2)
#m = length(vv)
#VV = matrix(0,m,2)
#for(i in 1:m){
#	v.tmp = vv[i]
#	A.tmp = replicate(2000,PCR.predict(n,p,K,w,u,v.tmp))
#	VV[i,] = rowMeans(A.tmp)
#}
#vv3 = VV
#vv3












