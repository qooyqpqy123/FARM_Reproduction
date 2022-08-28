X0 = read.csv('current.csv')
nnn = colnames(X0)
nnn = nnn[-1]
tcode = as.numeric(X0[1,-1]) ## length 127
X0 = X0[399:735,-1] ## dim 337*127
X0 = X0[221:337,]  ## dim 117*127
##X0 = X0[1:150,]
##-----------------------------------------------------------------------------
D = matrix(0,117-2,127) ## data transformation
for(j in 1:127){
	tmp.code = tcode[j]
	if(tmp.code==1){
		D[,j] = X0[-c(1,2),j]
	}
	if(tmp.code==2){
		tmp.v = X0[-1,j]-head(X0[,j],-1)
		D[,j] = tmp.v[-1]
	}
	if(tmp.code==3){
		tmp.v = X0[-1,j]-head(X0[,j],-1)
		tmp.v2 = tmp.v[-1]-head(tmp.v,-1)
		D[,j] = tmp.v2
	}
	if(tmp.code==4){
		D[,j] = log(X0[-c(1,2),j])
	}
	if(tmp.code==5){
		tmp.v = log(X0[-1,j])-log(head(X0[,j],-1))
		D[,j] = tmp.v[-1]
	}
	if(tmp.code==6){
		tmp.v = log(X0[-1,j])-log(head(X0[,j],-1))
		tmp.v2 = tmp.v[-1]-head(tmp.v,-1)
		D[,j] = tmp.v2
	}
	if(tmp.code==7){
		tmp.v = X0[-1,j]/head(X0[,j],-1)
		D[,j] = tmp.v[-1]-head(tmp.v,-1)
	}
}
D = D[,-69] ## dim 115*126
nnn = nnn[-69]