library(multtest)

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
############################################################################

DEP=MhD
deltaSeed=0.2
deltaAdd=0.2
minPro=20
minMic=5
linkageMethod="complete"

############################################################################

Protein=read.table("Data3/Protein.txt",header=TRUE,na.strings='.')
Microbe=read.table("Data3/Microbe.txt",header=TRUE)
Sample=read.table("Data3/Sample.txt",header=TRUE)
Pro=Protein[,-1]
Mic=Microbe[,-1]
nPro=dim(Pro)[2]
nMic=dim(Mic)[2]
nSam=dim(Pro)[1]

cocluster=list()
for (k in 1:nMic) {
  distChar=as.character(k)
  distName=paste(c("Dist",distChar,".txt"),collapse="")
  distSeed=read.table(file=paste("matrixMhD",distName,sep="/"))
  branch=cutree(hclust(as.dist(distSeed),method=linkageMethod),h=1-deltaSeed)
  branchTab=table(branch)
  branchCount=rbind(as.numeric(names(branchTab)),as.vector(branchTab))
  branchCount=matrix(branchCount[,branchCount[2,]>=minPro],nrow=2)
  nSeed=dim(branchCount)[2]
  if (nSeed>0) {
    for (indSeed in 1:nSeed) {
      clusterPro=as.vector(which(branch==branchCount[1,indSeed]))
      sampleRef=rbind(rep(Mic[,k],branchCount[2,indSeed]),c(as.matrix(Pro[,clusterPro])))
      pvalue=c()
      for (pj in c(1:nMic)[-k]) {
        for (pi in clusterPro) {
          sampleSca=rbind(Mic[,pj],Pro[,pi])
          quality=QualityIndex(sampleRef,sampleSca,DEP)
          pvalue=c(pvalue,2*pnorm(-abs(quality-1/2),0,sqrt((1+1/branchCount[2,indSeed])/(12*nSam))))
        }
      }
      pvalueAdj=mt.rawp2adjp(pvalue,"Holm")
      checkMat=cbind(matrix(as.numeric(pvalueAdj$adjp[order(pvalueAdj$index),2]>=deltaAdd),ncol=(nMic-1)),rep(1,length(clusterPro)))
      indexMat=c(c(1:nMic)[-k],k)
      checkMat=checkMat[,order(indexMat)]
      clusterMic=which(colSums(checkMat)==length(clusterPro))
      if (length(clusterMic)>=minMic) {
        cocluster[[length(cocluster)+1]]=list(clusterPro,clusterMic)
      }
    }
  }
}

cocluster

for (i in 1:length(cocluster)) {
  write.table(names(Pro)[cocluster[[i]][[1]]],"cocluster.xls",append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(t(names(Mic)[sort(cocluster[[i]][[2]])]),"cocluster.xls",append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE)
}
