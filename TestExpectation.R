set.seed(2)
alpha_0<-0
alpha_1<-0.5
alpha_2<-0.2
alpha_3<-0.7
alpha_4<-5
alpha_5<-2
alpha_6<-3

alpha<-2

X1<-rnorm(mean=0,1000)
X2<-rnorm(mean=0,1000)
X3<-rnorm(mean=0,1000)
eps<-rnorm(mean=0,1000)
nu<-rnorm(mean=0,sd=1,1000)
t<-pmax(X2*alpha+eps,0)

Y=alpha_0+alpha_1*X1+alpha_2*t+alpha_3*X2*t+alpha_4*X3+alpha_5*X1*X3+alpha_6*X1^2*X3+nu


t2<-round(t,digits = 2)

tvals<-sort(unique(t2))
condmean<-rep(0,length(tvals))

for(i in 1:length(tvals)){
  condmean[i]<-mean(Y[t2==tvals[i]])
}


plot(tvals,condmean)

real_t<-alpha_0+alpha_2*tvals+alpha_3*(alpha/(alpha^2+1))*tvals^2
actual_response<-alpha_2*tvals

plot(X2,t)
plot(tvals,real_t)

library(ggplot2)

ggplot()+geom_point(aes(x=tvals,y=real_t),col="red")+geom_point(aes(x=tvals,y=condmean),col="blue")+geom_point(aes(x=tvals,y=actual_response),col="yellow")

X.mod<-matrix(c(X1,X2,X3,t,X1*X2,X2*X3,X1*X2*X3,X1*X3,X1*t,X2*t,X3*t,t^2,X1*t^2,X2*t^2,X3*t^2),ncol=15)


lasso.model<-cv.glmnet(X.mod,Y)

condmeanLassoY<-rep(0,length(tvals))

for(i in 1:length(tvals)){
  treatmentvec<-rep(tvals[i],length(t))
  Xpred<-matrix(c(X1,X2,X3,treatmentvec,X1*X2,X2*X3,X1*X2*X3,X1*X3,X1*treatmentvec,X2*treatmentvec,X3*treatmentvec,treatmentvec^2,X1*treatmentvec^2,X2*treatmentvec^2,X3*treatmentvec^2),ncol=15)
  LassoYpred<-predict.cv.glmnet(lasso.model,newx=Xpred)
  condmeanLassoY[i]<-mean(LassoYpred)
}


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="yellow")+geom_point(aes(x=tvals,y=condmeanLassoY),col="red")+geom_point(aes(x=tvals,y=condmean),col="blue")


linear.model<-lm(Y~X.mod)
linmodY<-predict(linear.model)

condmeanLinY<-rep(0,length(tvals))

for(i in 1:length(tvals)){
  condmeanLinY[i]<-mean(linmodY[t2==tvals[i]])
}


ggplot()+geom_point(aes(x=tvals,y=actual_response),col="yellow")+geom_point(aes(x=tvals,y=condmeanLinY),col="red")+geom_point(aes(x=tvals,y=condmeanLassoY),col="blue")

summary(linear.model)
