library(truncreg)
library(gbm)
library(neuralnet)
library(rpart)
generateData<-function(n,p,alpha){
  X = matrix(runif(n * p), n, p)
  nu<-rnorm(n)
  t<-runif(n,min=rep(0,n),max = (alpha*(X[,1])))
  for(i in 1:p){
    for(j in 1:p){
      X<-cbind(X,X[,i]*X[,j])  
    }
  }
  X<-cbind(X,t)
  X<-cbind(X,t^2)
  X<-cbind(X,t^3)
  return(X)
}

generateCoefficients<-function(treatStrength, sparsity,X){
  coef<-dim(X)[2]
  beta<-rnorm(coef)
  if(sparsity>0){
    indexzero<-sample(1:(coef-3))[1:floor(sparsity*(coef-3))]
    beta[indexzero]=0
  }
  beta[coef-2]<-treatStrength
  beta[coef-1]<-(1/8)*treatStrength
  beta[coef]<-(1/32)*treatStrength
  return(beta)
}

generateOutcome<-function(X,beta){
  n<-dim(X)[1]
  eps<-rnorm(n)
  Y= X%*%beta + eps
  return(Y)
}
n=5000
alpha=5
treatStrength=2
p=20
trtstr=treatStrength
sparsity=0.8

X<-generateData(n,p,alpha = alpha)
coef<-dim(X)[2]
beta<-generateCoefficients(treatStrength,sparsity,X)
trtstr<-treatStrength
Y<-generateOutcome(X,beta)
t<-X[,(p+p^2+1)]
plot(X[,1],t)
plot(t,Y)
maxT<-max(t)
NoTreatOutcome<-sum(beta[1:(p+p^2)]*(1/2))
tseq<-seq(0,maxT,0.01)
actual_response<-tseq*beta[coef-2]+beta[coef-1]*tseq^2+beta[coef]*tseq^3+NoTreatOutcome
plot(tseq,actual_response)


#AVERAGE TREE PREDICTIONS:
avgLossWgtd<-rep(0,1000)
avgLossNonWgtd<-rep(0,1000)

for(i in 1:1000){
  index<-sample(n,replace=TRUE)
  Xsamp<-X[index,]
  Ysamp<-Y[index]
  tsamp<-Xsamp[,(p+p^2+1)]
  propensities<-(1/(alpha*Xsamp[,1]))
  uncondense<- (-log(tsamp/alpha))/alpha
  Y.star<-Ysamp*uncondense/propensities
  treemodel<-rpart(Y.star~tsamp,cp=0.002,method='anova')
  preds.tree<-(predict(treemodel,newdata = data.frame(tsamp=tseq)))
  cutoffs<-as.vector(sort(treemodel$splits[,4]))
  predsAvg<-rep(0,length(tseq))
  for(k in 1:length(cutoffs)){
    if(k==1){
      predsAvg[tseq<cutoffs[k]]<- mean(Ysamp[tseq<cutoffs[k]])
    }else if((k>1)&(k<length(cutoffs))){
      predsAvg[tseq>cutoffs[k-1]&tseq<cutoffs[k]]<- mean(Ysamp[tseq>cutoffs[k-1]&tseq<cutoffs[k]])
    }else{
      predsAvg[tseq>cutoffs[k]]<-mean(Ysamp[tseq>cutoffs[k]])
    }
  }
  Y.star<-Ysamp/propensities
  treemodel<-rpart(Y.star~tsamp,cp=0.002,method='anova')
  preds.tree.unweighted<-(predict(treemodel,newdata = data.frame(tsamp=tseq)))
  avgLossNonWgtd[i]<-mean((actual_response-preds.tree.unweighted)^2)
  avgLossWgtd[i]<-mean((actual_response-preds.tree)^2)
  if(i==1){
    avgWeightedTree<-data.frame(preds=preds.tree)
    avgNaiveTree<-data.frame(preds=preds.tree.unweighted)
    avgOutcomeAvgs<-data.frame(preds=predsAvg)
  }else{
    avgWeightedTree$preds<- avgWeightedTree$preds + (preds.tree-avgWeightedTree$preds)/i
    avgNaiveTree$preds<- avgNaiveTree$preds + (preds.tree-avgNaiveTree$preds)/i
    avgOutcomeAvgs$preds<-avgOutcomeAvgs$preds + (predsAvg-avgOutcomeAvgs$preds)/i
  }
  
  
}

ggplot()+geom_point(aes(x=tseq,y=avgWeightedTree$preds))+geom_point(aes(x=tseq,y=actual_response),col="red")+geom_point(aes(x=tseq,y=avgNaiveTree$preds),col="blue")








