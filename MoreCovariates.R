library(np)
n=10000
p=20
X = matrix(rnorm(n * p), n, p)
eps<-rnorm(n)
nu<-rnorm(n)
alpha<-1
t<-exp(X[,1]*alpha+nu)

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

#Setting equal to zero:
indexzero<-sample(1:p^2)[1:floor(0.5*p^2)]

beta<-rnorm(coef)
beta[coef]<-10

beta[indexzero]=0

Y= X%*%beta + eps






t2<-round(t,digits = 2)

tvals<-sort(unique(t2))
condmean<-rep(0,length(tvals))

for(i in 1:length(tvals)){
  condmean[i]<-mean(Y[t2==tvals[i]])
}

actual_response<-beta[coef]*tvals


library(ggplot2)

ggplot()+geom_point(aes(x=tvals,y=actual_response),col="red")+geom_point(aes(x=tvals,y=condmean),col="blue")



lasso.model<-cv.glmnet(X,Y)

condmeanLassoY<-rep(0,length(tvals))

for(i in 1:length(tvals)){
  treatmentvec<-rep(tvals[i],length(t))
  Xpred<-X[,1:(p^2+p)]
  for(j in 1:p){
    Xpred<-cbind(Xpred,Xpred[,j]*treatmentvec)  
  }
  
  Xpred<-cbind(Xpred,t)
  LassoYpred<-predict.cv.glmnet(lasso.model,newx=Xpred)
  condmeanLassoY[i]<-mean(LassoYpred)
}


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="yellow")+geom_point(aes(x=tvals,y=condmeanLassoY),col="red")+geom_point(aes(x=tvals,y=condmean),col="blue")





ridge.model<-cv.glmnet(X,Y,alpha=0)

condmeanRidge<-rep(0,length(tvals))

for(i in 1:length(tvals)){
  treatmentvec<-rep(tvals[i],length(t))
  Xpred<-X[,1:p^2]
  for(j in 1:p){
    Xpred<-cbind(Xpred,Xpred[,j]*treatmentvec)  
  }

  Xpred<-cbind(Xpred,t)
  RidgeYpred<-predict.cv.glmnet(ridge.model,newx=Xpred)
  condmeanRidge[i]<-mean(RidgeYpred)
}


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="yellow")+geom_point(aes(x=tvals,y=condmeanRidge),col="red")+geom_point(aes(x=tvals,y=condmean),col="blue")

