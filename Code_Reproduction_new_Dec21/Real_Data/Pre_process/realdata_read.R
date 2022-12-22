Z0 = read.csv('current.csv')
nnn = colnames(Z0)
nnn = nnn[-1]
choose_time=1 #set choose_time=1: we choose dataset in August 2010 to February 2020; set choose_time=2: we choose time February 1992 to October 2007
tcode = as.numeric(Z0[1,-1]) ## length 127
Y0 = Z0[399:735,-1] ## dim 337*127 #Choose the time table that works best for us.
if (choose_time==1){
X0 = Y0[221:337,]
D = matrix(0,117-2,127)}else
{X0=Y0[1:191,]
D = matrix(0,191-2,127)}
#If we choose X0=Y0[1:189,], instead and conduct the remaining steps, we obtain the dataset realdata_189_126.csv
##-------------------Do Data Transformation----------------------------------------------------------
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
D = D[,-69] ## output the data, we obtain realdata_115_126.csv, realdata_189_126.csv can be achieved similarly.
nnn = nnn[-69]#name of every line
