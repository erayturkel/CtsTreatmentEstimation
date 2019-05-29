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


n=20000
alpha=5
treatStrength=3
p=20
trtstr=5
sparsity=0.2
X<-generateData(n,p,alpha = alpha)
beta<-generateCoefficients(treatStrength,sparsity,X)
trtstr<-treatStrength
Y<-generateOutcome(X,beta)
t<-X[,(p+p^2+1)]
plot(X[,1],t)
plot(t,Y)
maxT<-max(t)
NoTreatOutcome<-sum(beta[1:(p+p^2)]*(1/2))
actual_response<-t*treatStrength+0.5*treatStrength*t^2+(1/8)*treatStrength*t^3+NoTreatOutcome
plot(t,actual_response)






propensities<-(1/(alpha*X[,1]))
uncondense<- (-log(t/alpha))/alpha



Y.star<-Y*uncondense/propensities

treemodel<-rpart(Y.star~t,cp=0.001)

preds.tree<-(predict(treemodel))


ggplot()+geom_point(aes(x=t,y=preds.tree))+geom_point(aes(x=t,y=actual_response),col="red")
