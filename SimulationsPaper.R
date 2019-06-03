library(nnet)
library(gbm)
library(neuralnet)
library(rpart)
set.seed(1)
generateData<-function(n,p,alpha){
  X = matrix(runif(n * p), n, p)
  nu<-rnorm(n)
  t<-runif(n,min=rep(0,n),max=X[,1]*alpha)
  for(i in 1:p){
    for(j in 1:p){
      X<-cbind(X,X[,i]*X[,j])  
    }
  }
  Xnorm = matrix(rnorm(n * 10), n, 10)
  for(i in 1:10){
    X<-cbind(X,Xnorm[,i]*t)  
  }
  X<-cbind(X,t)
  return(X)
}

generateCoefficients<-function(sparsity,X){
  coef<-dim(X)[2]-1
  beta<-rnorm(coef)
  if(sparsity>0){
    indexzero<-sample(1:(coef))[1:floor(sparsity*(coef))]
    beta[indexzero]=0
  }
  return(beta)
}

generateOutcome<-function(X,beta,t){
  n<-dim(X)[1]
  p<-dim(X)[2]
  eps<-rnorm(n)
  Y= X[,-p]%*%beta + t^2 -0.1*t^3+eps
  return(Y)
}


nsim=1000

SimResultsFrame<-data.frame(lasso=rep(0,nsim),ridge=rep(0,nsim), guido=rep(0,nsim), newmethod=rep(0,nsim))
num=2000 #Number of data poins in a simulated sample
vars=30 #Number of variables 
a=3 #Confounding term
s=0.05 #Sparseness parameter


for(sim in 1:nsim){
  Xtrain<-generateData(num,vars,a)
  ttrain<-Xtrain[,dim(Xtrain)[2]]
  beta<-generateCoefficients(s,Xtrain)
  Ytrain<-generateOutcome(Xtrain,beta,ttrain)
  Xtest<-generateData(num,vars,a)
  ttest<-Xtest[,dim(Xtest)[2]]
  Ytest<-generateOutcome(Xtest,beta,ttest)
  testprop<-(1/(a*Xtest[,1]))
  uncondense.test<- (-log(ttest/a))/a
  Ystar.test<-Ytest*uncondense.test/testprop
  NoTreatOutcome<-sum(beta[1:vars]*(1/2)) + sum(beta[vars:(vars+vars^2)]*0.25)
  actualResponse<-ttest^2-0.1*ttest^3+NoTreatOutcome
  propensities<-(1/(a*Xtrain[,1]))
  uncondense<- (-log(ttrain/a))/a
  Y.star<-Ytrain*uncondense/propensities
  ##Validation set for choosing Complexity parameter
  index<-sample(num,replace=FALSE)
  trainindex<-index[1:(floor(0.7*num))]
  valindex<-(-trainindex)
  Xval<-Xtrain[valindex]
  Yval<-Y.star[valindex]
  tval<-ttrain[valindex]
  optval<-1000000
  optreg<-0
  for(reg in seq(0.0001,0.01,0.0005)){
    cvframe=data.frame(Y=Y.star[trainindex],t=ttrain[trainindex])
    treemodel<-rpart(Y~t,data=cvframe,cp=reg)
    preds.tree<-(predict(treemodel,newdata=data.frame(Y=Yval,t=tval)))
    val<-mean((Y.star[valindex]-preds.tree)^2)
    if(val<optval){
      optval=val
      optreg=reg
    }
  }
  treemodel<-rpart(Y~t,data=cvframe,cp=optreg)
  testframe=data.frame(Y=Ystar.test,t=ttest)
  preds.tree<-(predict(treemodel,newdata=testframe))
  num_leaves <- length(unique(preds.tree))
  testframe$leaf <- factor(preds.tree, labels = seq(num_leaves))
  classes<-treemodel$where
  classtest<-testframe$leaf
  dataclasses=data.frame(c=classes,Xtrain[trainindex,(1:vars)])
  if(num_leaves>2){
    classpred<-multinom(as.factor(c)~.,dataclasses,MaxNWts =10000000)
    propscores<-as.vector(apply(predict(classpred,newdata=data.frame(Xtest[,(1:vars)]),type="probs"), 1, max))  
    predsdiscrete=rep(0,length(ttest))
    for(c in unique(classtest)){
      predsdiscrete[classtest==c]=sum(Ytest[classtest==c]/propscores[classtest==c])/sum((1/propscores)[classtest==c])
    }
     }else if(num_leaves==1){
       predsdiscrete=rep(mean(Ytest),length(ttest))
    }else{
    classpred<-multinom(as.factor(c)~.,dataclasses,MaxNWts =10000000)
    propscores<-as.vector(predict(classpred,newdata=data.frame(Xtest[,(1:vars)]),type="probs"))
    base<-unique(testframe$leaf)[1]
    treats<-testframe$leaf
    propscores[treats==base]<-(1-propscores[treats==base])
    predsdiscrete=rep(0,length(ttest))
    for(c in unique(classtest)){
      predsdiscrete[classtest==c]=sum(Ytest[classtest==c]/propscores[classtest==c])/sum((1/propscores)[classtest==c])
    }
  }
  Y.star<-Ytrain
  datatrain<-data.frame(Y.star=Y.star,t=ttrain,p=propensities,tsq=ttrain^2,psq=propensities^2,int=ttrain*propensities)
  linmodel<-lm(Y.star~.,data=datatrain)
  predsoutcome<-data.frame(t=ttest,pred=0)
  for( treat in ttest){
    datapredict<-data.frame(Y.star=Ytest,t=treat,p=testprop,tsq=treat^2,psq=testprop^2,int=treat*testprop)
    preds.lm<-(predict(linmodel,newdata=datapredict))
    predsoutcome[predsoutcome$t==treat,]$pred<-mean(preds.lm)
    
  }
  # lasso.model<-cv.glmnet(Xtrain,Ytrain,nfolds = 4)
  # condmeanLassoY<-rep(0,length(ttest))
  # for(i in 1:length(ttest)){
  #   treatmentvec<-rep(ttest[i],length(ttest))
  #   Xpred<-cbind(Xtest[,(1:(dim(Xtest)[2]-1))],treatmentvec)
  #   LassoYpred<-predict.cv.glmnet(lasso.model,newx=Xpred)
  #   condmeanLassoY[i]<-mean(LassoYpred)
  # }
  # ridge.model<-cv.glmnet(Xtrain,Ytrain,alpha=0,nfolds = 4)
  # condmeanRidgeY<-rep(0,length(ttest))
  # for(i in 1:length(ttest)){
  #   treatmentvec<-rep(ttest[i],length(ttest))
  #   Xpred<-cbind(Xtest[,(1:(dim(Xtest)[2]-1))],treatmentvec)
  #   RidgeYpred<-predict.cv.glmnet(ridge.model,newx=Xpred)
  #   condmeanRidgeY[i]<-mean(RidgeYpred)
  # }
  # 
  
  # ggplot()+geom_point(aes(x=predsoutcome$t,y=predsoutcome$pred),col="green")+
  #   geom_point(aes(x=ttest,y=actualResponse),col="red")+
  #   geom_point(aes(x=ttest,y=predsdiscrete),col="blue")
    # geom_point(aes(x=ttest,y=condmeanLassoY),col="yellow")+
    # geom_point(aes(x=ttest,y=condmeanRidgeY),col="orange")
  
  # SimResultsFrame[sim,]$lasso<-mean((actualResponse-condmeanLassoY)^2)
  # SimResultsFrame[sim,]$ridge<-mean((actualResponse-condmeanRidgeY)^2)
  SimResultsFrame[sim,]$guido<-mean((actualResponse-predsoutcome$pred)^2)
  SimResultsFrame[sim,]$newmethod<-mean((actualResponse-predsdiscrete)^2)
  
}




results005sparse<-SimResultsFrame