library(nnet)
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
  beta[coef-1]<-1/2*treatStrength
  beta[coef]<-1/8*treatStrength
  return(beta)
}

generateOutcome<-function(X,beta){
  n<-dim(X)[1]
  eps<-rnorm(n)
  Y= X%*%beta + eps
  return(Y)
}


n=10000
alpha=5
treatStrength=3
p=20
trtstr=3
sparsity=0.2
X<-generateData(n,p,alpha = alpha)
beta<-generateCoefficients(treatStrength,sparsity,X)
trtstr<-treatStrength
Y<-generateOutcome(X,beta)
t<-X[,(p+p^2+1)]
plot(X[,1],t)
plot(t,Y)
maxT<-max(t)
coef=length(beta)
NoTreatOutcome<-sum(beta[1:p]*(1/2)) + sum(beta[p:p^2]*0.25)
actual_response<-t*beta[coef-2]+beta[coef-1]*t^2+beta[coef]*t^3+NoTreatOutcome

propensities<-(1/(alpha*X[,1]))
uncondense<- (-log(t/alpha))/alpha



Y.star<-Y*uncondense/propensities

treemodel<-rpart(Y.star~t,cp=0.005)

preds.tree<-(predict(treemodel))


ggplot()+geom_point(aes(x=t,y=preds.tree))+geom_point(aes(x=t,y=actual_response),col="red")

classes<-treemodel$where
classpred<-multinom(classes~X[,(1:p)])
propscores<-as.vector(apply(fitted(classpred), 1, max))
predsdiscrete=rep(0,length(t))
for(c in classes){
  predsdiscrete[classes==c]=sum(Y[classes==c]/propscores[classes==c])/sum((1/propscores)[classes==c])
}


ggplot()+geom_point(aes(x=t,y=preds.tree))+geom_point(aes(x=t,y=actual_response),col="red")+geom_point(aes(x=t,y=predsdiscrete),col="blue")



lasso.model<-cv.glmnet(X,Y)

condmeanLassoY<-rep(0,length(t))

for(i in 1:length(t)){
  treatmentvec<-rep(t[i],length(t))
  Xpred<-cbind(X[,(1:(p+p^2))],treatmentvec)
  Xpred<-cbind(Xpred,treatmentvec^2)
  Xpred<-cbind(Xpred,treatmentvec^3)
  LassoYpred<-predict.cv.glmnet(lasso.model,newx=Xpred)
  condmeanLassoY[i]<-mean(LassoYpred)
}

ggplot()+geom_point(aes(x=predsoutcome$t,y=predsoutcome$pred))+geom_point(aes(x=t,y=actual_response),col="red")+geom_point(aes(x=t,y=condmeanLassoY),col="blue")