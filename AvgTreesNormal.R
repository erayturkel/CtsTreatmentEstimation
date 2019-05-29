library(np)
library(np)
n=50000
p=20
X = matrix(rnorm(n * p), n, p)
eps<-rnorm(n)
nu<-rnorm(n)
alpha<-2
t<-rnorm(n,mean=X[,1]*alpha)
for(i in 1:p){
  for(j in 1:p){
    X<-cbind(X,X[,i]*X[,j])  
  }
}

for(i in 1:p){
  X<-cbind(X,X[,i]*t)  
}
X<-cbind(X,t)
coef<-dim(X)[2]
indexzero<-sample(1:p^2)[1:floor(0.5*p^2)]
beta<-rnorm(coef)
beta[coef]<-10
beta[indexzero]=0
Y= X%*%beta + eps
plot(t,Y)
hist(t)
tvals<-seq(-5,5,0.01)
tseq<-seq(-5,5,0.01)
actual_response<-beta[coef]*tvals

ggplot()+geom_point(aes(x=tvals$tsamp,y=actual_response),col="red")


avgLossWgtdNorm<-rep(0,1000)
avgLossNonWgtdNorm<-rep(0,1000)

for(i in 1:1000){
  index<-sample(n,replace=TRUE)
  Xsamp<-X[index,]
  Ysamp<-Y[index]
  tsamp<-Xsamp[,(2*p+p^2+1)]
  meanVector<-Xsamp[,1]*alpha
  propensities<-rep(0,length(tsamp))
  for(j in 1:length(tsamp)){
    mu<-meanVector[j]
    propensities[j]<- dnorm(tsamp[j]-mu)
  }
  densityofT<- exp(-tsamp^2/(2*alpha^2+2))/(sqrt(2*pi)*sqrt(alpha^2+1))
  Y.star<-Ysamp*densityofT/propensities
  treemodel<-rpart(Y.star~tsamp,cp=0.005,method='anova')
  tvals<-seq(-5,5,0.01)
  tvals<-data.frame(tsamp=tvals)
  preds.tree<-(predict(treemodel,newdata = tvals))
  predsAvg<-rep(0,nrow(tvals))
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
  treemodel<-rpart(Y.star~tsamp,cp=0.005,method='anova')
  preds.tree.unweighted<-(predict(treemodel,newdata = tvals))
  avgLossNonWgtdNorm[i]<-mean((actual_response-preds.tree.unweighted)^2)
  avgLossWgtdNorm[i]<-mean((actual_response-preds.tree)^2)
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






ggplot()+geom_point(aes(x=tvals$tsamp,y=avgWeightedTree$preds))+geom_point(aes(x=tvals$tsamp,y=actual_response),col="red")+geom_point(aes(x=tvals$tsamp,y=avgNaiveTree$preds),col="blue")+coord_cartesian(ylim=c(-150,150))



#ggplot()+geom_point(aes(x=tvals$tsamp,y=avgWeightedTree$preds))+geom_point(aes(x=tvals$tsamp,y=actual_response),col="red")+geom_point(aes(x=tvals$tsamp,y=avgNaiveTree$preds),col="blue")+geom_point(aes(x=tvals$tsamp,y=avgOutcomeAvgs$preds),col="yellow")+coord_cartesian(ylim=c(-150,150))

