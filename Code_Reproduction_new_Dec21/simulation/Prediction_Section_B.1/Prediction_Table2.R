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
  SigmaX = tcrossprod(X)/n#covariance matrix
  eigenX = eigen(SigmaX)
  eigvec = eigenX$vectors
  eigvalue = eigenX$values
  K.hat = which.min(diff(log(eigvalue[1:100])))  #select the number of factors
  F.hat = eigvec[,1:K.hat]*sqrt(n) #estimated factors
  B.hat = t(t(F.hat)%*%X)/n #factor loading
  U.hat = X - F.hat%*%t(B.hat) #idiosyncratic component
  ##------------------Model Estimation for Lasso -----------------------------------------------------------
  cv.fit.X = cv.glmnet(X,Y,intercept=FALSE)	 #Fit lass
  lambda.fit.X = cv.fit.X$lambda.min
  fit.X = glmnet(X,Y,intercept=FALSE,lambda=lambda.fit.X)
  beta.hat.X = as.vector(fit.X$beta) #fitted beta for Lasso
  ##----------------- Model Estimation for FARM -----------------------------------------------------
  lmY.F = lm(Y~F.hat-1) #Fit factor loading
  gamma.hat = coef(lmY.F) #fitted gamma
  Y.tilde = resid(lmY.F)
  cv.fit.U = cv.glmnet(U.hat,Y.tilde,intercept=FALSE)	
  lambda.fit.U = cv.fit.U$lambda.min
  fit.U = glmnet(U.hat,Y.tilde,intercept=FALSE,lambda=lambda.fit.U) #fit FARM
  beta.hat.U = as.vector(fit.U$beta) #fitted value
  ##-------------------Model prediction: New data---------------------------------------------------------
  F = rnorm(K) 
  U = rnorm(p)
  X = B%*%F + U
  E = rnorm(1)*sigma.0
  Y = sum(F*varphi.0) + sum(X*beta.0) + E
  ##----------------------------------------------------------------------------
  lmX.B = lm(X~B.hat-1)
  F.new = coef(lmX.B)
  U.new = resid(lmX.B)
  ##------------------ Model prediction results---------------------------------------------------------
  X.predict = Y - sum(X*beta.hat.X) #prediction of Lasso
  U.predict = Y - sum(F.new*gamma.hat) - sum(U.new*beta.hat.U) #prediction for FARM
  tmp.X = X.predict^2
  tmp.U = U.predict^2
  return(c(tmp.X,tmp.U))  #return prediction errors for (LASSO, FARM)
}

#-----------------Main codes -----------------
n = 200; K = 5; w = 0.8; vv = 2;p=1000 #setting on dimention, signal strengths
vv_seq=c(0.2,0.4,0.6,0.8,1.0,1.2,1.4)
MEAN1=matrix(0,length(vv_seq),2)
for (i in c(1:length(vv_seq))){ #for different signal strength, compute the prediction results of FARM and LASSO
  # print(i)
  u=vv_seq[i]
  u=0
  A.tmp = replicate(500,PCR.predict(n,p,K,w,u,vv))
  MEAN1[i,] = rowMeans(A.tmp)
  write.csv(MEAN1,'MEAN_p.csv') # Table 2 in appendix.
}













