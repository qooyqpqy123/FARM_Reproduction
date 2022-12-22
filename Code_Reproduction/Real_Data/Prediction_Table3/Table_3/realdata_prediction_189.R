###Read in Data####
D = read.csv('realdata_189_126.csv')
D = D[,-1]
D<-data.frame(D)
D=D[,-c(81,83)]
nnn = colnames(D)
##-----------------------------------------------------------------
#Loading required packages
##------------------------------------------------------------------
library(glmnet)
library(elasticnet)
library(randomForest)
set.seed(100)
T = 90  ## moving window approach, window size                         
M = nrow(D)-T  ## predict sample size 
N = ncol(D) #number of columns
#---------------------------------------------------------------------------------------
#########################################################################################################
#### The following values are initialized for storing prediction difference, R^2, predicted values, obtained from different methods
## like FARM, FLasso, Lasso, Sample Mean, PCR, Ridge Elastic Net, Random Forest, respectively. Details are in Section 5.4 of the paper.##########
#########################################################################################################
Pred.FARM = matrix(0,M,N)                       ## FARM for estimation and prediction
Pred.FARMLASSO = matrix(0,M,N)                  ## Flasso for prediction
Pred.LASSO = matrix(0,M,N)                      ## Lasso for prediction
Pred.MEAN = matrix(0,M,N)                       ## MEAN for prediction
Pred.PCR = matrix(0,M,N)				                ## PCR for prediction
Pred.RIDGE=matrix(0,M,N)                        ##Ridge Regression for prediction
Pred.ElA=matrix(0,M,N)                          ##Elastic Net for prediction
Pred.RF=matrix(0,M,N)                           ##Random Forest for prediction
#-----------------------------------
#Store out-of-sample R^2 for several methods mentioned above
#------------------------------------
R2.FARM = numeric(N)           #FARM  
R2.LASSO = numeric(N)          #LASSO
R2.PCR = numeric(N)           #PCR
R2.FARMLASSO = numeric(N)     #FLasso
R2.RIDGE = numeric(N)          #Ridge
R2.ELA=numeric(N)              #Elastic Net
R2.RF=numeric(N)               # Random Forest
#-----------------------------------
#Store prediction values for several methods mentioned above
#------------------------------------
Prvalue_FARM=matrix(0,M,N)     #FARM
Prvalue_Lasso=matrix(0,M,N)    #LASSO
Prvalue_FarmLasso=matrix(0,M,N) #FLasso
Prvalue_PCR=matrix(0,M,N) #PCR
Prvalue_MEAN=matrix(0,M,N) #Sample MEAN
Prvalue_RIDGE=matrix(0,M,N)#RIdge Regression
Prvalue_ELA=matrix(0,M,N) #Elastic Net
Prvalue_RF=matrix(0,M,N) #Random Forest
#----------------------------------------------
#Predictions is conducted Next
#-----------------------------------------------
j=49 #set j=81 is the prediction result for GS5
X = D[,-j] #Use rest of the data to be covariate
Y = D[,j] #Response variable
for(i in 1:M){
  idx = i:(i+T-1) #Moving window prediction, every window has length T
  x = X[idx,] #Training Data Covariate
  y = Y[idx]  #Training Data Response
  x.new = X[i+T,] #Prediction Covariate
  y.new = Y[i+T] #Prediction Response
  ##--------------------------------------------------------------------
  #In the following we first do data normalization##
  ##-------------------------------------------------------------------
  x.mu = colMeans(x)
  x.sd = as.numeric(apply(x,2,sd))
  y.mu = mean(y)
  Prvalue_MEAN[i,j]=	y.mu 
  y.sd = sd(y)
  x = t((t(x)-x.mu)/x.sd)                               ## Data normalization, we standardize every column of X to be mean zero and sd 1
  y = (y-y.mu)/y.sd
  x.new = (x.new-x.mu)/x.sd
  x.new<-unlist(x.new)
  y.new = (y.new-y.mu)/y.sd
  ##-------------------------------------------------------------------
  ###In the following, we conduct several prediction comparison using sample Mean, Lasso, Ridge, Elastic Net, FARM, PCR, FLasso
  ##--------------------------------------------------------------------
  Pred.MEAN[i,j] = (y.new*y.sd)^2   #Sample Mean Prediction
  ##--------------------------------------------------------------------
  ##Lasso
  ##--------------------------------------------------------------------
  cv.fit.x = cv.glmnet(x,y,intercept=FALSE)	
  lambda.fit.x = cv.fit.x$lambda.min
  fit.x = glmnet(x,y,intercept=FALSE,lambda=lambda.fit.x)  ## Lasso Estimation
  beta.hat.x = as.vector(fit.x$beta)
  Pred.LASSO[i,j] = (y.new-sum(x.new*beta.hat.x))^2*y.sd^2
  Prvalue_Lasso[i,j]=sum(x.new*beta.hat.x)*y.sd+y.mu ## Lasso Prediction
  ##--------------------------------------------------------------------
  ## Ridge estimation
  ##--------------------------------------------------------------------
  cv.fit.x = cv.glmnet(x,y,intercept=FALSE,alpha=0)	
  lambda.fit.x = cv.fit.x$lambda.1se
  fit.x = glmnet(x,y,intercept=FALSE,lambda=lambda.fit.x,alpha=0)  ## Ridge estimation
  beta.hat.xr = as.vector(fit.x$beta)
  Pred.RIDGE[i,j] = (y.new-sum(x.new*beta.hat.xr))^2*y.sd^2
  Prvalue_RIDGE[i,j]=sum(x.new*beta.hat.xr)*y.sd+y.mu #Ridge Prediction
  ##-------------------------------------------------------------------
  ## Elastic Net
  ##--------------------------------------------------------------------
  alpha_list=seq(0,1,0.1)
  min_list=c()
  lambda_list=c()
  for (k in c(1:length(alpha_list))){
    cv.fit.x1 = cv.glmnet(x,y,intercept=FALSE,alpha=alpha_list[k])	
    min_list=c(min_list,min(cv.fit.x1$cvm))                    
    lambda_list=c(lambda_list,cv.fit.x1$lambda.min)
  }
  fit.x1 = glmnet(x,y,intercept=FALSE,lambda=lambda_list[which.min(min_list)],alpha=alpha_list[which.min(min_list)])  ## Elastic Net estimation
  beta.hat.xe = as.vector(fit.x1$beta)
  Pred.ElA[i,j] = (y.new-sum(x.new*beta.hat.xe))^2*y.sd^2
  Prvalue_ELA[i,j]=sum(x.new*beta.hat.xe)*y.sd+y.mu #Elastic Net Prediction
  ## --------------------------------------------------------------------
  ##Factor Estimation
  ## --------------------------------------------------------------------
  Sigma.x = tcrossprod(x)/T      #covariance matrix
  eigenx = eigen(Sigma.x)        
  eigvec = eigenx$vectors        #Eigen vector
  eigvalue = eigenx$values       #Eigen values
  K.hat = max(2,which.min(diff(log(eigvalue[1:10]))))      ## Factor estimation
  F.hat = eigvec[,1:K.hat]*sqrt(T)        #Estimated Factor
  B.hat = t(t(F.hat)%*%x)/T          #Estimated Factor Loading
  U.hat = x-F.hat%*%t(B.hat)        #Estimated idiosyncratic component
  ##--------------------------------------------------------------------
  #FARM
  ##--------------------------------------------------------------------
  lmY.F = lm(y~F.hat-1)           
  gamma.hat = coef(lmY.F)    #Estimate \gamma vector in the FARM model and PCR
  Y.tilde = resid(lmY.F)                             
  cv.fit.U = cv.glmnet(U.hat,Y.tilde,intercept=FALSE)	 
  lambda.fit.U = cv.fit.U$lambda.min #tuning parameter selection
  fit.U = glmnet(U.hat,Y.tilde,intercept=FALSE,lambda=lambda.fit.U) ##FARM Estimation
  beta.hat.U = as.vector(fit.U$beta) #Obtain the fitted beta for FARM
  ##------------------------------------------------------------------------------
  #FLasso and PCR
  ##------------------------------------------------------------------------------
  lmx.B = lm(x.new~B.hat-1)
  f.new = coef(lmx.B)
  u.new = resid(lmx.B)
  Pred.FARM[i,j] = (y.new-sum(f.new*gamma.hat)-sum(u.new*beta.hat.U))^2*y.sd^2 #FARM prediction Difference with the true value
  Prvalue_FARM[i,j]=(sum(f.new*gamma.hat)+sum(u.new*beta.hat.U))*y.sd+y.mu  #FARM predicted value
  Pred.FARMLASSO[i,j] = (y.new-sum(x.new*beta.hat.U))^2*y.sd^2 #FLasso Prediction Difference with the true value
  Prvalue_FarmLasso[i,j]=(sum(x.new*beta.hat.U))*y.sd+y.mu  #FLasso predicted value
  Pred.PCR[i,j] = (y.new-sum(f.new*gamma.hat))^2*y.sd^2     #PCR prediction Difference with the true value
  Prvalue_PCR[i,j]=sum(f.new*gamma.hat)*y.sd+y.mu           #PCR predicted value
}
#————————————————————————————————————————————————————————————
##Output R^2 for Table 3
#------------------------------------------------------------
R2.FARM[j] = 1-sum(Pred.FARM[,j])/sum(Pred.MEAN[,j])	
print(R2.FARM[j])                    #Output R^2 value for FARM
R2.FARMLASSO[j] = 1-sum(Pred.FARMLASSO[,j])/sum(Pred.MEAN[,j])	
print(R2.FARMLASSO[j])              #Output R^2 value for Lasso
R2.PCR[j] = 1-sum(Pred.PCR[,j])/sum(Pred.MEAN[,j])
print(R2.PCR[j])                   #Output R^2 value for PCR
R2.LASSO[j] = 1-sum(Pred.LASSO[,j])/sum(Pred.MEAN[,j])
print(R2.LASSO[j])                 #Output R^2 value for LASSO
R2.RIDGE[j] = 1-sum(Pred.RIDGE[,j])/sum(Pred.MEAN[,j])
print(R2.RIDGE[j])                 #Output R^2 value for RIDGE
R2.ELA[j] = 1-sum(Pred.ElA[,j])/sum(Pred.MEAN[,j])
print(R2.ELA[j])                 #Output R^2 value for Elastic Net

###------------------------------------------
##Finally, we compare it with random forest
###----------------------------------------

for(i in 1:M){
  idx = i:(i+T-1)          #Moving Window, window size
  x = X[idx,]              #Trainind Data, Covariate
  y = Y[idx]              #Train Data, Response
  x.new = X[i+T,]         #Test Data
  y.new = Y[i+T]          #Test Data
  ##--------------------------------------------------------------------
  x.mu = colMeans(x)
  x.sd = as.numeric(apply(x,2,sd))
  y.mu = mean(y)
  Prvalue_MEAN[i,j]=	y.mu 
  y.sd = sd(y)
  x = t((t(x)-x.mu)/x.sd)                               ## Do the same Data normalization as above
  y = (y-y.mu)/y.sd
  x.new = (x.new-x.mu)/x.sd
  x.new<-unlist(x.new)
  y.new = (y.new-y.mu)/y.sd
  ##-----------------------------------------------------------------------
  rf<- randomForest(x,y,ntree=100)                          #Train random Forest
  Pred.RF[i,j]=(y.new-predict(rf,newdata=x.new))^2*y.sd^2    ##Random Forest Predicted difference with true value
  Prvalue_RF[i,j]=predict(rf,newdata=x.new)*y.sd+y.mu }      #Random Forest predicted values
R2.RF[j]=1-sum(Pred.RF[,j])/sum(Pred.MEAN[,j])            #Random Forest out-of-sample R^2
print(R2.RF[j])                                           #print the R^2 of Random Forest



