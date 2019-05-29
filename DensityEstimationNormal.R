library(np)
linmodel<-truncreg(t~ -1+X[,1:p])

linmodel<-lm(log(t)~ -1+ X[,1:p])
# 
# coefsFromModel<-linmodel$coefficients
# sigma<-coefsFromModel[21]
# coefsFromModel<-coefsFromModel[1:20]
# coefsFromModel<-as.vector(coefsFromModel)
# sigma<-as.numeric(sigma)
# meanVector<-X[,1:p]%*%coefsFromModel
# 
meanVector<-X[,1]*alpha

propensities<-rep(0,length(t))
for(i in 1:length(t)){
  if(t[i]!=0){
    mu<-meanVector[i]
    # nom<-dnorm((t[i]-mu)/sigma)
    # denom<- sigma
    propensities[i]<- (1/(t[i]*sqrt(2*pi)))*exp(-(log(t[i])-mu)^2/2)
  }else{
    mu<-meanVector[i]  
    propensities[i]<-pnorm(-mu/sigma)
  }
}
hist(propensities)
hist(t)
plot(t,propensities)

plot(t[t<10],propensities[t<10])
hist(propensities[t==0])

uncondense<-npuniden.boundary(X=t,Y=t,a=0,b=Inf)

#densityofT<-uncondense$f


densityofT<- exp((-1/4)*log(t)^2)/(t*2*sqrt(pi))

Y.star<-Y[propensities>0.025&propensities<0.975]*densityofT[propensities>0.025&propensities<0.975]/propensities[propensities>0.025&propensities<0.975]


tees<-t[propensities>0.025&propensities<0.975]

treemodel<-rpart(Y.star~tees,cp=0.005)

plot(tees,Y.star)

tvals<-seq(0,5,0.01)
tvals<-data.frame(tees=tvals)


preds.tree<-(predict(treemodel,newdata = tvals))

actual_response<-beta[coef]*tvals$tees

predsActual<-rep(0,nrow(tvals))
for(i in 1:nrow(tvals)){
  localpred=preds.tree[i]
  predsActual[i]<-mean(Y[preds.tree==localpred])
}


ggplot()+geom_point(aes(x=tvals$tees,y=preds.tree))+geom_point(aes(x=tvals$tees,y=actual_response),col="red")


Ydata<-Y[propensities>0.05&propensities<0.95]
predDataFrame<-data.frame(Y=Ydata,t=t[propensities>0.05&propensities<0.95],props=propensities[propensities>0.05&propensities<0.95], preds=predict(treemodel))
predDataFrame$means<-NA
for(i in unique(predDataFrame$preds)){
  predDataFrame[predDataFrame$preds==i,]$means<-((predDataFrame[predDataFrame$preds==i,]$Y)*(predDataFrame[predDataFrame$preds==i,]$props)/sum(predDataFrame[predDataFrame$preds==i,]$props))
}


ggplot()+geom_point(aes(x=predDataFrame$t,y=predDataFrame$means))+geom_point(aes(x=tvals$tees,y=actual_response),col="red")
