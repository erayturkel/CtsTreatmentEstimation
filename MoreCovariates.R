n=10000
p=20
X = matrix(rnorm(n * p), n, p)
eps<-rnorm(n)
nu<-rnorm(n)
alpha<-0.8
t<-pmax(X1*alpha+eps,0)

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

beta<-rnorm(coef)
beta[coef]<-floor(max(beta))

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
  Xpred<-X[1:p^2]
  for(j in 1:p){
    Xpred<-cbind(Xpred,Xpred[,j]*treatmentvec)  
  }
  Xpred<-cbind(Xpred,t)
  LassoYpred<-predict.cv.glmnet(lasso.model,newx=Xpred)
  condmeanLassoY[i]<-mean(LassoYpred)
}


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="yellow")+geom_point(aes(x=tvals,y=condmeanLassoY),col="red")+geom_point(aes(x=tvals,y=condmean),col="blue")

