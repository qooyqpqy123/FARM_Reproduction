#simulation for Test 2
set.seed(100)
#generate data
#set.seed(200)
library("SIS")
quant_chisq<-c()
for (q in 1:10){
  boostrap<-c()
  for(j in 1:500){
    boostrap<-c(boostrap,rchisq(1,3))
  }
  quant_chisq<-c(quant_chisq,quantile(boostrap,0.95))
}
c_alpha<-mean(quant_chisq)

#Y1 responses, method:1==ISIS;2==SIS
selectfeature<-function(Y1,X){
  cov<-2*(X[1:(n/2),]-rowMeans(X[1:(n/2),]))%*%t((X[1:(n/2),]-rowMeans(X[1:(n/2),])))/n #construct factors
  eigvec<-eigen(cov)$vectors
  eigval<-eigen(cov)$values
  K1<-which.min(diff(log(eigval[1:100]))) #select K
  hatf<-sqrt(n/2)*eigvec[,1:K1] #hatF
  hatB<-t(2/n*t(hatf)%*%X[1:(n/2),]) #hat B
  hatU<-X[1:(n/2),]-hatf%*%t(hatB) #hat U
Y_new<-Y1[1:(n/2)]-2/n*hatf%*%t(hatf)%*%Y1[1:(n/2)] #
sis<-SIS(hatU,Y_new)
select1<-sis$ix #ISIS
#S.hat0.ISIS = SIS0$ix; S.hat0.SIS = SIS0$ix0
cov2<-2*(X[((n/2)+1):n,]-rowMeans(X[((n/2)+1):n,]))%*%t(X[((n/2)+1):n,]-rowMeans(X[((n/2)+1):n,]))/n
eigvec2<-eigen(cov2)$vectors
eigval2<-eigen(cov2)$values
K<-which.min(diff(log(eigval2[1:100]))) #select K
hatf2<-sqrt(n/2)*eigvec2[,1:K] #hat F using second half of data
hatB2<-t(2/n*t(hatf2)%*%X[((n/2)+1):n,]) #hat B
hatU2<-X[((n/2)+1):n,]-hatf2%*%t(hatB2) #hat U
P_f<-2*hatf2%*%t(hatf2)/n
P_x<-X[(((n/2)+1):n),select1]%*%solve(t(X[(((n/2)+1):n),select1])%*%(X[((n/2)+1):n,select1]))%*%t(X[(((n/2)+1):n),select1])
P_u<-hatU2[,select1]%*%solve(t(hatU2[,select1])%*%hatU2[,select1])%*%t(hatU2[,select1])
Q1<-t(Y1[((n/2)+1):n])%*%(P_f+P_u-P_x)%*%(P_f+P_u-P_x)%*%Y1[(n/2+1):n,]
select2<-sis$ix0 #SIS
P_x<-X[(((n/2)+1):n),select2]%*%solve(t(X[(((n/2)+1):n),select2])%*%(X[((n/2)+1):n,select2]))%*%t(X[(((n/2)+1):n),select2])
P_u<-hatU2[,select2]%*%solve(t(hatU2[,select2])%*%hatU2[,select2])%*%t(hatU2[,select2])
Q2<-t(Y1[((n/2)+1):n])%*%(P_f+P_u-P_x)%*%(P_f+P_u-P_x)%*%Y1[(n/2+1):n,]
#X2=cbind(F.hat2,U.hat2)
c(Q1,Q2)
}
rcv<-function(Y1,X){
  Idx1=1:(n/2)
  Idx2=(n/2+1):n
  lis=list(Idx1,Idx2)
  var_ISIS<-c()
  var_SIS<-c()
  for (Idx in lis){
  cov<-2*(X[Idx,]-rowMeans(X[Idx,]))%*%t((X[Idx,]-rowMeans(X[Idx,])))/n #construct factors
  eigvec<-eigen(cov)$vectors
  eigval<-eigen(cov)$values
  K1<-which.min(diff(log(eigval[1:100]))) #select K
  hatf<-sqrt(n/2)*eigvec[,1:K1] #hatF
  hatB<-t(2/n*t(hatf)%*%X[Idx,]) #hat B
  hatU<-X[Idx,]-hatf%*%t(hatB) #hat U
  Y_new<-Y1[Idx]-2/n*hatf%*%t(hatf)%*%Y1[Idx] #
  sis<-SIS(hatU,Y_new)
  select1<-sis$ix #ISIS  
  select2<-sis$ix0 #SIS
  
  cov2<-(X[-Idx,]-rowMeans(X[-Idx,]))%*%t(X[-Idx,]-rowMeans(X[-Idx,]))/(n/2)
  eigvec2<-eigen(cov2)$vectors
  eigval2<-eigen(cov2)$values
  K<-which.min(diff(log(eigval2[1:100]))) #select K
  hatf2<-sqrt(n/2)*eigvec2[,1:K] #hat F using second half of data
  hatB2<-t(2/n*t(hatf2)%*%X[-Idx,]) #hat B
  hatU2<-X[-Idx,]-hatf2%*%t(hatB2) #hat U  
  
  
X0=cbind(hatf2,hatU2)
X0.ISIS = X0[,c(1:K,select1+K)]
X0.SIS = X0[,c(1:K,select2+K)]
lm0.ISIS = lm(Y1[-Idx]~X0.ISIS-1) 
lm0.SIS = lm(Y1[-Idx]~X0.SIS-1)
Q0.ISIS.H1 = sum((resid(lm0.ISIS))^2)
sigma.hat0.ISIS = Q0.ISIS.H1/((n/2) - K - length(select1))
Q0.SIS.H1 = sum((resid(lm0.SIS))^2)
sigma.hat0.SIS = Q0.SIS.H1/((n/2) - K - length(select2))
#X2=cbind(F.hat2,U.hat2)
var_ISIS<-c(var_ISIS,sigma.hat0.ISIS)
var_SIS<-c(var_SIS,sigma.hat0.SIS)
  }
  print(var_ISIS)
  print(var_SIS)
c(mean(var_ISIS),mean(var_SIS))
}

n=200
p=500
K=3
beta=c(rep(0.7,4),rep(0,p-4)) #true beta
#gamma1=c(rep(m/10,K))#runif(K,-m/10,m/10)#alternative
sigma=0.5
B<-matrix(runif(p*K,-1,1),nrow=p) #factor loading



#inference using second half of data
#X: design matrix, select1: selected indexes from SIS/ISIS,Y1: responses


mat1<-matrix(nrow=4,ncol=6)
mat2<-matrix(nrow=4,ncol=6)
mat3<-matrix(nrow=4,ncol=6)
mat4<-matrix(nrow=4,ncol=6)
mat5<-matrix(nrow=4,ncol=6)
mat6<-matrix(nrow=4,ncol=6)
mat_isis<-matrix(nrow=4,ncol=6)
mat_sis<-matrix(nrow=4,ncol=6)
for (m in c(0,0.6,1.2,1.8,2.4,3)){
  idx<-c()
  idx2<-c()
  idx3<-c()
  idx4<-c()
  idx5<-c()
  idx6<-c()
  mean_isis<-c()
  mean_sis<-c()
  gamma1=c(rep(m/10,K))
  for(f in c(1:4)){
    isis1<-c()
    isis2<-c()
    isis3<-c()
    isis5<-c()
    isis4<-c()
    isis6<-c()
    sigma_isis<-c()
    sigma_sis<-c()
    ##select features and do inference###
    for (i in 1:1000){
      
      F<-matrix(rnorm(n*K),nrow=n) #factor
      U<-matrix(rnorm(n*p),nrow=n) #noise
      
      X=F%*%t(B)+U #design matrix
      #Y1=X%*%beta+rnorm(n,mean=0,sd=0.5) #Null
      Y2=F%*%gamma1+X%*%beta+runif(n,-sqrt(3)/2,sqrt(3)/2)#rnorm(n,mean=0,sd=0.5)#runif(n,-sqrt(3)/2,sqrt(3)/2)#rnorm(n,mean=0,sd=0.5) #alternative
      output1<-try(selectfeature(Y2,X),silent=TRUE)
      if ('try-error' %in% class(output1)) {
        output1<-c(0,0)}
      sig_est<-try(rcv(Y2,X),silent=TRUE)
      if ('try-error' %in% class(sig_est)) {
        sig_est<-c(0.25,0.25)}
      var_est1<-sig_est[1]#(output1[3]+output2[3])/2
      sigma_isis<-c(sigma_isis,var_est1)
      write.csv(sigma_isis,"sigma_isis.csv")
      var_est2<-sig_est[2]
      sigma_sis<-c(sigma_sis,var_est2)
      write.csv(sigma_sis,"sigma_sis.csv")
      isis1<-c(isis1,output1[1]/var_est1)#sigma_isis)
      isis2<-c(isis2,output1[1]/0.25)#sigma_sis) #ISIS_SIS
      isis3<-c(isis3,output1[1]/var_est2)#sigma_sis)
      isis4<-c(isis4,output1[2]/var_est1)
      isis5<-c(isis5,output1[2]/0.25)
      isis6<-c(isis6,output1[2]/var_est2)
      #isis42<-c(isis42,output1[2]/0.25)#sigma_isis)#SIS_ISIS
    }
    mean_isis<-c(mean_isis,mean(sigma_isis))
    write.csv(mean_isis,"mean_isis.csv")
    mean_sis<-c(mean_sis,mean(sigma_sis))
    write.csv(mean_sis,"mean_sis.csv")
    idx<-c(idx,sum(isis1>c_alpha)/1000)
    write.csv(idx,"idx.csv")
    idx2<-c(idx2,sum(isis2>c_alpha)/1000)
    write.csv(idx2,"idx2.csv")
    idx3<-c(idx3,sum(isis3>c_alpha)/1000)
    write.csv(idx3,"idx3.csv")
    idx4<-c(idx4,sum(isis4>c_alpha)/1000)
    write.csv(idx4,"idx4.csv")
    idx5<-c(idx5,sum(isis5>c_alpha)/1000)
    write.csv(idx5,"idx5.csv")
    idx6<-c(idx6,sum(isis6>c_alpha)/1000)
    write.csv(idx6,"idx6.csv")
    
  }
  mat1[,(floor(m/0.6)+1)]<-idx
  mat2[,(floor(m/0.6)+1)]<-idx2
  mat3[,(floor(m/0.6)+1)]<-idx3
  mat4[,(floor(m/0.6)+1)]<-idx4
  mat5[,(floor(m/0.6)+1)]<-idx5
  mat6[,(floor(m/0.6)+1)]<-idx6
  mat_isis[,(floor(m/0.6)+1)]<-mean_isis
  mat_sis[,(floor(m/0.6)+1)]<-mean_sis
  write.csv(mat1,"mat1.csv")
  write.csv(mat2,"mat2.csv")
  write.csv(mat3,"mat3.csv")
  write.csv(mat4,"mat4.csv")
  write.csv(mat5,"mat5.csv")
  write.csv(mat6,"mat6.csv")
}
