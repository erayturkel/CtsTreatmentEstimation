set.seed(2)
n=2000
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
beta[coef]<-5

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


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="red")+geom_point(aes(x=tvals,y=condmeanLassoY),col="yellow")



ridge.model<-cv.glmnet(X,Y,alpha=0)

condmeanRidge<-rep(0,length(tvals))

for(i in 1:length(tvals)){
  treatmentvec<-rep(tvals[i],length(t))
  Xpred<-X[,1:(p+p^2)]
  for(j in 1:p){
    Xpred<-cbind(Xpred,Xpred[,j]*treatmentvec)  
  }

  Xpred<-cbind(Xpred,t)
  RidgeYpred<-predict.cv.glmnet(ridge.model,newx=Xpred)
  condmeanRidge[i]<-mean(RidgeYpred)
}


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="red")+geom_point(aes(x=tvals,y=condmeanRidge),col="green")+geom_point(aes(x=tvals,y=condmeanLassoY),col="yellow")+coord_cartesian(xlim = c(0,30))


meanVector<-X[,1]*alpha
propensities<-rep(0,length(t))
for(j in 1:length(t)){
  mu<-meanVector[j]
  propensities[j]<- dnorm(log(t[j])-mu)/t[j]
}
densityofT<- exp(-t^2/(2*alpha^2+2))/(sqrt(2*pi)*sqrt(alpha^2+1))
Y.star<-Y*densityofT/propensities


treemodel<-rpart(Y.star[propensities>0.05]~t[propensities>0.05],cp=0.005)

preds.tree.2<-(predict(treemodel))
num_leaves <- length(unique(preds.tree.2))
testframe=data.frame(Y=Y.star[propensities>0.05],X=X[propensities>0.05,],t=t[propensities>0.05])
testframe$leaf <- factor(preds.tree.2, labels = seq(num_leaves))



classes<-treemodel$where
classtest<-testframe$leaf
dataclasses=data.frame(c=classes,X[propensities>0.05,(1:(coef-1))])
classpred<-multinom(as.factor(c)~.,dataclasses,MaxNWts =10000000)
if(num_leaves>2){
  propscores<-as.vector(apply(predict(classpred,newdata=data.frame(X[propensities>0.05,(1:(coef-1))]),type="probs"), 1, max))  
  predsdiscrete=rep(0,length(t[propensities>0.05]))
  for(c in unique(classtest)){
    predsdiscrete[classtest==c]=sum(testframe$Y[classtest==c]/propscores[classtest==c])/sum((1/propscores)[classtest==c])
  }
}else{
  propscores<-as.vector(predict(classpred,newdata=data.frame(X[,(1:(coef-1))]),type="probs"))
  base<-unique(testframe$leaf)[1]
  treats<-testframe$leaf
  propscores[treats==base]<-(1-propscores[treats==base])
  predsdiscrete=rep(0,length(t[propensities>0.05]))
  for(c in unique(classtest)){
    predsdiscrete[classtest==c]=sum(testframe$Y[classtest==c]/propscores[classtest==c])/sum((1/propscores)[classtest==c])
  }
}


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="red")+geom_point(aes(x=tvals,y=condmeanRidge),col="green")+geom_point(aes(x=tvals,y=condmeanLassoY),col="yellow")+geom_point(aes(x=t[propensities>0.05],y=predsdiscrete),col="orange")+coord_cartesian(xlim = c(0,30))




ggplot()+geom_point(aes(x=t[propensities>0.05],y=preds.tree.2),col="red")+geom_point(aes(x=tvals,y=condmeanLassoY),col="blue")+geom_point(aes(x=tvals,y=actual_response),col="yellow")+coord_cartesian(ylim = c(0,100),xlim = c(0,10))



Y.star<-Y

datatrain<-data.frame(Y.star=Y.star,t=t,p=propensities,tsq=t^2,psq=propensities^2,int=t*p)
linmodel<-lm(Y.star~.,data=datatrain)

predsoutcome<-data.frame(t=t,pred=0)
for( treat in t){
  datapredict<-data.frame(Y.star=Y.star,t=treat,p=propensities,tsq=treat^2,psq=propensities^2,int=treat*propensities)
  preds.lm<-(predict(linmodel,newdata=datapredict))
  predsoutcome[predsoutcome$t==treat,]$pred<-mean(preds.lm)
  
}



ggplot()+geom_point(aes(x=predsoutcome$t,y=predsoutcome$pred),col="green")+geom_point(aes(x=t[propensities>0.05],y=preds.tree.2),col="red")+geom_point(aes(x=tvals,y=condmeanLassoY),col="blue")+geom_point(aes(x=tvals,y=actual_response),col="yellow")+coord_cartesian(ylim = c(0,100),xlim = c(0,10))





