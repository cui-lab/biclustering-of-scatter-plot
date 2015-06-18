
############################################################################################
MhD=function(w,x) {
  mu=apply(x,1,mean)
  sigma=cov(t(x))
  D=1/(1+colSums(solve(sigma)%*%(w-mu)*(w-mu)))
  return(D)
}

QualityIndex=function(F,G,dep) {
  Dfx=dep(F,F)
  Dfy=dep(G,F)
  Q=mean(rep(1,length(Dfy))%*%t(Dfx)-Dfy<=0)
  return(Q)
}

QTest=function(F,G,dep,nperm) {
  Qfg=QualityIndex(F,G,dep)
  Qgf=QualityIndex(G,F,dep)
  if (abs(Qfg-1/2)>abs(Qgf-1/2)) {
    TSO=Qfg
  } else {
    TSO=Qgf
  }
  Pool=cbind(F,G)
  TS=c()
  for (b in 1:nperm) {
    resample=sample(seq(1,(ncol(F)+ncol(G))),ncol(F),replace=FALSE)
    FP=Pool[,resample]
    GP=Pool[,-resample]
    QfgP=QualityIndex(FP,GP,dep)
    QgfP=QualityIndex(GP,FP,dep)
    if (abs(QfgP-1/2)>abs(QgfP-1/2)) {
      TS[b]=QfgP
    } else {
      TS[b]=QgfP
    }
  }
  plower=mean(TS<TSO)
  pupper=mean(TS>TSO)
  pvalue=2*min(plower,pupper)  
  return(pvalue)
}
############################################################################################
 
set.seed(0)
DEP=MhD
nPerm=500
k=1

############################################################################################

Protein=read.table("Data3/Protein.txt",header=TRUE,na.strings='.')
Microbe=read.table("Data3/Microbe.txt",header=TRUE)
Sample=read.table("Data3/Sample.txt",header=TRUE)
Pro=Protein[,-1]
Mic=Microbe[,-1]
nPro=dim(Pro)[2]
nMic=dim(Mic)[2]
nSam=dim(Pro)[1]

Dist=matrix(0,nrow=nPro,ncol=nPro)
for (i in 1:(nPro-1)) {
  for (j in (i+1):nPro) {
    FSam=rbind(Mic[,k],Pro[,i])
    GSam=rbind(Mic[,k],Pro[,j])
    Dist[j,i]=1-QTest(FSam,GSam,DEP,nPerm)
  }
}
 
write(t(Dist),"Dist1.txt",ncolumns=nPro)
