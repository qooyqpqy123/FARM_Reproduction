set.seed(400)
library(glmnet)
library(SIS)
library(MASS)
#library("clime")
#n=300 p=300

n=200
p=200
K=2

Phi<-matrix(nrow=K,ncol=K)
for (i in c(1:K)){
  for (j in c(1:K)){
    Phi[i,j]<-0.4^(abs(i-j)+1)
  }
}



top_mat<-matrix(nrow=p,ncol=p)
for (i in c(1:p)){
  for (j in c(1:p)){
    top_mat[i,j]=0.5^abs(i-j)
  }
}
beta_0=c(rep(0,3),rep(0,p-3)) #true beta
beta_1=c(rep(0.05,3),rep(0,p-3))
beta_2=c(rep(0.1,3),rep(0,p-3))
beta_3=c(rep(0.15,3),rep(0,p-3))
beta_4=c(rep(0.2,3),rep(0,p-3))
gamma1=0.5*rep(1,K)#alternative
sigma=0.5
B<-matrix(runif(p*K,-1,1),nrow=p) #factor loading

ratio<-c()
idx1<-c()
idx2<-c()
idx3<-c()
idx4<-c()
idx5<-c()
idx12<-c()
idx22<-c()
idx32<-c()
idx42<-c()
idx52<-c()
idx13<-c()
idx23<-c()
idx33<-c()
idx43<-c()
idx53<-c()
var_isis<-c()
var_sis<-c()
sum_idx1<-c()
sum_idx12<-c()
sum_idx13<-c()
  for (i in c(1:1000)){
  #norm(gamma1-t(B)%*%beta_0)
  #F<-matrix(rnorm(n*K),nrow=n) #factor
    F_start=rnorm(K)
    F_origin=F_start
    for (t in c(1:200)){
      #F_origin=F_start
      F_start=Phi%*%F_origin+rnorm(K)
      F_origin=F_start
    }
  F<-matrix(nrow=n,ncol=K) #factor
  for (r in c(1:n)){
      F_start=Phi%*%F_origin+rnorm(K)
      F_origin=F_start
      F[r,]<-F_start
    }
  U<-mvrnorm(n, rep(0,p),top_mat) #noise
  
  X=F%*%t(B)+U #design matrix
  X<-X-rowMeans(X)
  
  cov<-X%*%t(X)/n
  eigvec<-eigen(cov)$vectors
  eigval<-eigen(cov)$values
  K_est<-which.min(diff(log(eigval[1:10]))) #select K
  hatf<-sqrt(n)*eigvec[,1:K_est]
  hatB<-t(1/n*t(hatf)%*%X)
  hatU<-X-hatf%*%t(hatB)
  
  C<-as.matrix(diag(rep(1,p)))
  T<-c()
  for( j in 1:ncol(hatU)){
    fit_u = glmnet(hatU[,-j], hatU[,j], intercept=FALSE,
                   lambda=cv.glmnet(hatU[,-j], hatU[,j],intercept=FALSE)$lambda.1se)
    beta<-as.vector(fit_u$beta)
    C[j,-j]<--beta
    T<-c(T,1/n*sum((hatU[,j]-hatU[,-j]%*%beta)^2)+fit_u$lambda/2*sum(abs(beta)))
    if (j%%10==0){
      print(j)
    }
  }
  T1<-diag(1/T)
  Theta<-T1%*%C
  
  quant<-c() #threshold
  for (m in 1:5){
    boostrap<-c()
    for(j in 1:200){
      boostrap<-c(boostrap,max(abs(Theta%*%t(hatU)%*%rnorm(n,mean=0,sd=1))))
    }
    quant<-c(quant,quantile(boostrap,0.95))
  }
  c_alpha<-1/sqrt(n)*mean(quant)#0.95 threshold
  
  
  Idx = 1:(n/2)
  X0 = X[Idx,]; X2 = X[-Idx,]
  Sigma0 = tcrossprod(X0)*2/n; Sigma2 = tcrossprod(X2)*2/n
  eigen0 = eigen(Sigma0); eigen2 = eigen(Sigma2)
  eigvec0 = eigen0$vectors; eigvalue0 = eigen0$values
  eigvec2 = eigen2$vectors; eigvalue2 = eigen2$values
  K0 = which.min(diff(log(eigvalue0[1:100])))
  K2 = which.min(diff(log(eigvalue2[1:100])))
  F.hat0 = eigvec0[,1:K0]*sqrt(n/2); F.hat2 = eigvec2[,1:K2]*sqrt(n/2)
  B.hat0.T = t(F.hat0)%*%X0*2/n; B.hat2.T = t(F.hat2)%*%X2*2/n
  U.hat0 = X0 - F.hat0%*%B.hat0.T; U.hat2 = X2 - F.hat2%*%B.hat2.T
  
  rcv<-function(Y){
    Y0 = Y[Idx]; Y2 = Y[-Idx]
    
    tmp0 = tcrossprod(F.hat0); tmp2 = tcrossprod(F.hat2)
    Y0.new = Y0 - tmp0%*%Y0*2/n
    Y2.new = Y2 - tmp2%*%Y2*2/n
    SIS0 = SIS(U.hat0,Y0.new)
    SIS2 = SIS(U.hat2,Y2.new)
    S.hat0.ISIS = SIS0$ix; S.hat0.SIS = SIS0$ix0
    S.hat2.ISIS = SIS2$ix; S.hat2.SIS = SIS2$ix0 
    X0=cbind(F.hat0,U.hat0)
    X2=cbind(F.hat2,U.hat2)
    X0.ISIS = X0[,c(1:K0,S.hat2.ISIS+K0)]
    #print(c(1:K2,S.hat2.ISIS+K2))
    X0.SIS = X0[,c(1:K0,K0+S.hat2.SIS)]
    X2.ISIS = X2[,c(1:K2,K2+S.hat0.ISIS)]; X2.SIS = X2[,c(1:K2,K2+S.hat0.SIS)]
    lm0.ISIS = lm(Y0~X0.ISIS-1) 
    lm0.SIS = lm(Y0~X0.SIS-1)
    lm2.ISIS = lm(Y2~X2.ISIS-1) 
    lm2.SIS = lm(Y2~X2.SIS-1)
    Q0.ISIS.H1 = sum((resid(lm0.ISIS))^2)
    sigma.hat0.ISIS = Q0.ISIS.H1/(n/2 - K0 - length(S.hat2.ISIS))
    Q0.SIS.H1 = sum((resid(lm0.SIS))^2)
    sigma.hat0.SIS = Q0.SIS.H1/(n/2 - K0 - length(S.hat2.SIS))
    Q2.ISIS.H1 = sum((resid(lm2.ISIS))^2)
    sigma.hat2.ISIS = Q2.ISIS.H1/(n/2 - K2 - length(S.hat0.ISIS))
    Q2.SIS.H1 = sum((resid(lm2.SIS))^2)
    sigma.hat2.SIS = Q2.SIS.H1/(n/2 - K2 - length(S.hat0.SIS))
    sigma.hat.ISIS = mean(c(sigma.hat0.ISIS, sigma.hat2.ISIS))
    sigma.hat.SIS = mean(c(sigma.hat0.SIS, sigma.hat2.SIS))
  c(sigma.hat.ISIS,sigma.hat.SIS)}
  
    Y=F%*%gamma1+runif(n,-sqrt(3)/2,sqrt(3)/2)#runif(n,-sqrt(3)/2,sqrt(3)/2)#rnorm(n,0,sigma)#runif(n,-sqrt(3)/2,sqrt(3)/2)#rnorm(n,0,sigma)#runif(n,-sqrt(3)/2,sqrt(3)/2)#rnorm(n,0,sigma)#runif(n,-sqrt(3)/2,sqrt(3)/2) #null
    Y1=Y+U%*%beta_0 #alternative
    Y2=Y+U%*%beta_1
    Y3=Y+U%*%beta_2
    Y4=Y+U%*%beta_3
    Y5=Y+U%*%beta_4
    sd_alter1<-rcv(Y1)
    var_isis<-c(var_isis,sd_alter1[1])
    write.csv(var_isis,"isis.csv")
    var_sis<-c(var_sis,sd_alter1[2])
    write.csv(var_sis,"sis.csv")
    sd_alter2<-rcv(Y2)
    sd_alter3<-rcv(Y3)
    sd_alter4<-rcv(Y4)
    sd_alter5<-rcv(Y5)
    Y_new1<-Y1-1/n*hatf%*%t(hatf)%*%Y1
    Y_new2<-Y2-1/n*hatf%*%t(hatf)%*%Y2
    Y_new3<-Y3-1/n*hatf%*%t(hatf)%*%Y3
    Y_new4<-Y4-1/n*hatf%*%t(hatf)%*%Y4
    Y_new5<-Y5-1/n*hatf%*%t(hatf)%*%Y5
    fit1 = glmnet(hatU, Y_new1, intercept=FALSE,
                  lambda=cv.glmnet(hatU,Y_new1,intercept=FALSE)$lambda.1se)
    fit2 = glmnet(hatU, Y_new2, intercept=FALSE,
                  lambda=cv.glmnet(hatU,Y_new2,intercept=FALSE)$lambda.1se)
    fit3 = glmnet(hatU, Y_new3, intercept=FALSE,
                  lambda=cv.glmnet(hatU,Y_new3,intercept=FALSE)$lambda.1se)
    fit4 = glmnet(hatU, Y_new4, intercept=FALSE,
                  lambda=cv.glmnet(hatU,Y_new4,intercept=FALSE)$lambda.1se)
    fit5 = glmnet(hatU, Y_new5, intercept=FALSE,
                  lambda=cv.glmnet(hatU,Y_new5,intercept=FALSE)$lambda.1se)
    lambda_beta1<-fit1$lambda
    lambda_beta2<-fit2$lambda
    lambda_beta3<-fit3$lambda
    lambda_beta4<-fit4$lambda
    lambda_beta5<-fit5$lambda
    beta_hat1<-as.vector(fit1$beta) #fitted value
    beta_hat2<-as.vector(fit2$beta)
    beta_hat3<-as.vector(fit3$beta)
    beta_hat4<-as.vector(fit4$beta)
    beta_hat5<-as.vector(fit5$beta)
    d_beta1<-beta_hat1+1/n*Theta%*%t(hatU)%*%(Y1-hatU%*%beta_hat1) #debiased beta under the alternative
    d_beta2<-beta_hat2+1/n*Theta%*%t(hatU)%*%(Y2-hatU%*%beta_hat2)
    d_beta3<-beta_hat3+1/n*Theta%*%t(hatU)%*%(Y3-hatU%*%beta_hat3)
    d_beta4<-beta_hat4+1/n*Theta%*%t(hatU)%*%(Y4-hatU%*%beta_hat4)
    d_beta5<-beta_hat5+1/n*Theta%*%t(hatU)%*%(Y5-hatU%*%beta_hat5)
    idx1<-c(idx1,sum(sqrt(n)*max(abs(d_beta1))/sqrt(sd_alter1[1])>c_alpha))
    sum_idx1<-c(sum_idx1,sum(idx1)/i)
    write.csv(sum_idx1,"sum_idx1.csv")
    idx2<-c(idx2,sum(sqrt(n)*max(abs(d_beta2))/sqrt(sd_alter2[1])>c_alpha))
    idx3<-c(idx3,sum(sqrt(n)*max(abs(d_beta3))/sqrt(sd_alter3[1])>c_alpha))
    idx4<-c(idx4,sum(sqrt(n)*max(abs(d_beta4))/sqrt(sd_alter4[1])>c_alpha))
    idx5<-c(idx5,sum(sqrt(n)*max(abs(d_beta5))/sqrt(sd_alter5[1])>c_alpha))
    write.csv(idx1,"idx1.csv")
    write.csv(idx2,"idx2.csv")
    write.csv(idx3,"idx3.csv")
    write.csv(idx4,"idx4.csv")
    write.csv(idx5,"idx5.csv")
    idx12<-c(idx12,sum(sqrt(n)*max(abs(d_beta1))/0.5>c_alpha))
    sum_idx12<-c(sum_idx12,sum(idx12)/i)
    write.csv(sum_idx12,"sum_idx12.csv")
    idx22<-c(idx22,sum(sqrt(n)*max(abs(d_beta2))/0.5>c_alpha))
    idx32<-c(idx32,sum(sqrt(n)*max(abs(d_beta3))/0.5>c_alpha))
    idx42<-c(idx42,sum(sqrt(n)*max(abs(d_beta4))/0.5>c_alpha))
    idx52<-c(idx52,sum(sqrt(n)*max(abs(d_beta5))/0.5>c_alpha))
    write.csv(idx12,"idx12.csv")
    write.csv(idx22,"idx22.csv")
    write.csv(idx32,"idx32.csv")
    write.csv(idx42,"idx42.csv")
    write.csv(idx52,"idx52.csv")
    idx13<-c(idx13,sum(sqrt(n)*max(abs(d_beta1))/sqrt(sd_alter1[2])>c_alpha))
    sum_idx13<-c(sum_idx13,sum(idx13)/i)
    write.csv(sum_idx13,"sum_idx13.csv")
    idx23<-c(idx23,sum(sqrt(n)*max(abs(d_beta2))/sqrt(sd_alter1[2])>c_alpha))
    idx33<-c(idx33,sum(sqrt(n)*max(abs(d_beta3))/sqrt(sd_alter1[2])>c_alpha))
    idx43<-c(idx43,sum(sqrt(n)*max(abs(d_beta4))/sqrt(sd_alter1[2])>c_alpha))
    idx53<-c(idx53,sum(sqrt(n)*max(abs(d_beta5))/sqrt(sd_alter1[2])>c_alpha))
    write.csv(idx13,"idx13.csv")
    write.csv(idx23,"idx23.csv")
    write.csv(idx33,"idx33.csv")
    write.csv(idx43,"idx43.csv")
    write.csv(idx53,"idx53.csv")
  }

