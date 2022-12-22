library(glmnet)
#-----------Main Function: PCR.predict: Compute the prediction performance of factor regression and latent factor regression respectively ----
#Input n: number of observations
#p dimension
#K number of factors
#w: signal strength for non-zero beta
#u: signal strength for \phi
#v: signal strength for factor loading
PCR.predict = function(n,p,K,w,u,v){
  #------------Generate Data ---------------
  F = matrix(rnorm(n*K),n,K) #factors
  U = matrix(rnorm(n*p),n,p) #idiosyncratic component
  B = matrix(runif(p*K,-1,1),p,K)*v #factor loading matrix
  X = F%*%t(B) + U #covariate
  beta.0 = c(rep(w,20),rep(0,p-20)) #true signal beta
  varphi.0 = runif(K,-1,1)*u #true signal \phi
  sigma.0 = 0.5 #noise level
  E = rnorm(n)*sigma.0 
  Y = F%*%varphi.0 + U%*%beta.0 + E #response variable 
  #-------------Use covariance matrix to estimate factors -----------
  SigmaX = tcrossprod(X)/n #covariance matrix
  eigenX = eigen(SigmaX) 
  eigvec = eigenX$vectors #eigenvectors
  eigvalue = eigenX$values #eigenvalues
  K.hat = which.min(diff(log(eigvalue[1:100]))) #select the number of factors
  F.hat = eigvec[,1:K.hat]*sqrt(n) #estimate factos
  B.hat = t(t(F.hat)%*%X)/n #estimate factor loadings
  U.hat = X - F.hat%*%t(B.hat) #idiosyncratic component
  lmY.F = lm(Y~F.hat-1) 
  gamma.hat = coef(lmY.F) #estimate gamma for PCR
  Y.tilde = resid(lmY.F) 
  cv.fit.U = cv.glmnet(U.hat,Y.tilde,intercept=FALSE)	
  lambda.fit.U = cv.fit.U$lambda.min
  fit.U = glmnet(U.hat,Y.tilde,intercept=FALSE,lambda=lambda.fit.U)
  beta.hat.U = as.vector(fit.U$beta) #fit beta for FARM
  ##-------------New Data for prediction---------------------------------------------------------------
  F = rnorm(K)
  U = rnorm(p)
  X = B%*%F + U
  E = rnorm(1)*sigma.0
  Y = sum(F*varphi.0) + sum(U*beta.0) + E
  ##----------------------------------------------------------------------------
  lmX.B = lm(X~B.hat-1)
  F.new = coef(lmX.B)
  U.new = resid(lmX.B)
  ##--------------Prediction results--------------------------------------------------------------
  F.predict = Y - sum(F.new*gamma.hat) #prediction of factor regression
  U.predict = Y - sum(F.new*gamma.hat) - sum(U.new*beta.hat.U)  #prediction for FARM
  tmp.F = F.predict^2  
  tmp.U = U.predict^2
  return(c(tmp.F,tmp.U)) #return prediction errors for (PCR, FARM)
}
#-----------------Main codes -----------------
n = 200; K = 5; u=0.8;p=1000; vv=2 #setting on dimention, signal strengths
w_seq=c(0.2,0.4,0.6,0.8,1.0,1.2,1.4)


MEAN=matrix(0,length(w_seq),2) 
for (i in c(1:length(w_seq))){ #for different signal strength, compute the prediction results of FARM and factor regression
  w=w_seq[i]
  A.tmp = replicate(500,PCR.predict(n,p,K,w,u,vv))
  MEAN[i,] = rowMeans(A.tmp)
  write.csv(MEAN,'MEAN_2.csv') # Table 1 in appendix.
}

