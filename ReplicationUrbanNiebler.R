library(np)
library(rpart)
library(ggplot2)
library(nnet)
cov<-data.frame(dens=x$Pop,HHinc=x$PerCapitaHHInc,college=x$per_collegegrads/100,over65=x$PercentOver65,black=x$PercentBlack,hisp=x$PercentHispanic,white=x$PercentWhite,voted=x$PercentVoted,rural=x$Rural,RepShare=x$RepubShare)
y<-x$Cont
t<-x$TotAds
plot(t,y)
t<-as.vector(t)


uncondenseT<-npuniden.boundary(X=t,Y=t,a=0,b=Inf)
set.seed(1)
ind<-sample(nrow(cov),2000,replace=TRUE)
covtrain<-cov[ind,]
ttrain<-t[ind]
bwold<-npcdensbw(xdat=covtrain,ydat=ttrain,nmulti = 50,bwmethod = "cv.ml",cxkertype = "epanechnikov", cykertype = "epanechnikov",ftol = 0.01, tol=0.01,remin = TRUE)

newbws<-npcdensbw(xdat=covtrain,ydat=ttrain,nmulti = 1,bws = bwold , bwmethod = "cv.ml",cxkertype = "epanechnikov", cykertype = "epanechnikov")
  
  
bws<-newbws

condens<-npcdens(bws = bws,txdat = cov,tydat = t)

conditionalDenseT<-condens$condens


Y.star<- y*(1000*uncondenseT$f)/(1000*conditionalDenseT)
treemodel<-rpart(Y.star~t,cp=0.0005)
preds.tree<-predict(treemodel,newdata = data.frame(t=seq(0,20000,1)))
ggplot()+geom_point(aes(x=seq(0,20000,1),y=preds.tree))
classes<-treemodel$where
dataclasses=data.frame(c=classes,cov)
classpred<-multinom(as.factor(c)~.,dataclasses,MaxNWts =10000000,maxit=100)
propscores<-as.vector(apply(predict(classpred,type="probs"), 1, max))  
predsdiscrete=rep(0,length(seq(0,max(t),1)))
for(c in unique(classes)){
  mint<-min(t[classes==c])
  maxt<-max(t[classes==c])
  predsdiscrete[mint:maxt]<-sum(y[classes==c]/propscores[classes==c])/sum((1/propscores)[classes==c])
}
predsdiscrete[predsdiscrete==0]<-NA
ggplot()+geom_point(aes(x=seq(0,max(t),1),y=predsdiscrete))+coord_cartesian(ylim=c(0,20))





datatrain<-data.frame(Y=y,t=t,p=conditionalDenseT,tsq=t^2,psq=conditionalDenseT^2,int=t*conditionalDenseT)
linmodel<-lm(Y~.,data=datatrain)
predsoutcome<-data.frame(t=seq(0,20000,100),pred=0)
for( treat in seq(0,20000,100)){
  testprop<-predict(condens,exdat=cov,eydat=rep(treat,nrow(cov)),txdat=cov,tydat=t)
  datapredict<-data.frame(Y=y,t=treat,p=testprop,tsq=treat^2,psq=testprop^2,int=treat*testprop)
  preds.lm<-(predict(linmodel,newdata=datapredict))
  predsoutcome[predsoutcome$t==treat,]$pred<-mean(preds.lm)
}



PaperResults<-rep(0,length(seq(0,max(t),1)))
PaperResults[(1000:length(PaperResults))]<-6.3


ggplot()+
  geom_point(aes(x=seq(0,max(t),1),y=PaperResults),col="red")+
  geom_point(aes(x=seq(0,20000,100),y=predsoutcome$pred),col="blue")+
  xlab("t")+ylab("Contributions (thousands)")+
  coord_cartesian(ylim=c(0,30))



propensities<-conditionalDenseT
index<-sample(length(y))
trainindex<-index[1:8000]
valindex<-index[8001:11000]
testindex<-index[11001:16265]


EstT<-t[trainindex]
EstY<-y[trainindex]
EstP<-propensities[trainindex]
EstFrame<-data.frame(EstY=EstY,EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)

ValT<-t[valindex]
ValY<-y[valindex]
ValP<-propensities[valindex]
ValFrame<-data.frame(EstY=ValY,EstT=ValT,EstP=ValP,EstTsq=ValT^2,EstPsq=ValP^2,Intr=ValT*ValP)
MSEmin=10000000
optlayer=0
optlayer2=0
for(i in 5:20){
  for(k in 1:i){
    modelNNET<-neuralnet(EstY~EstT+EstP+EstTsq+EstPsq+Intr, data= EstFrame, hidden = c(i,k), threshold = 0.2,act.fct='tanh',learningrate=0.0001, algorithm="backprop",stepmax = 1e+05,rep=1,lifesign = 'full',act.fct = 'logistic')
    MSE<-mean((ValY-compute(modelNNET,ValFrame[,-1])$net.result)^2)
    if(MSE<MSEmin){
      MSEmin=MSE
      optlayer=i
      optlayer2=j
    }
  }
}



modelNNET<-neuralnet(EstY~EstT+EstP+EstTsq+EstPsq+Intr, data= EstFrame, hidden = c(optlayer,optlayer2), threshold = 5,learningrate=0.00001, algorithm="backprop",stepmax = 2e+05,rep=1,lifesign = 'full',act.fct = 'logistic',linear.output=TRUE)

TestT<-t[testindex]
TestY<-y[testindex]
TestP<-propensities[testindex]
Testcov<-cov[testindex,]
TestFrame<-data.frame(EstY=TestY,EstT=TestT,EstP=TestP,EstTsq=TestT^2,EstPsq=TestP^2,Intr=TestT*TestP)
predsoutcomeNNET<-data.frame(t=seq(0,20000,100),pred=0)
for( treat in seq(0,20000,100)){
  TestP<-fitted(npcdens(exdat=Testcov, eydat=rep(treat,nrow(Testcov)), bws=bws,txdat=cov,tydat=t))
  TestT<-rep(treat,nrow(Testcov))
  TestFrame<-data.frame(EstY=TestY,EstT=TestT,EstP=TestP,EstTsq=TestT^2,EstPsq=TestP^2,Intr=TestT*TestP)
  preds.nnet<-mean(compute(modelNNET,TestFrame[,-1])$net.result)
  predsoutcomeNNET[predsoutcomeNNET$t==treat,]$pred<-mean(preds.nnet)
}



library(gbm)

index<-sample(length(y))
trainindex<-index[1:10000]
testindex<-index[10001:16265]

EstT<-t[trainindex]
EstY<-y[trainindex]
EstP<-propensities[trainindex]
EstFrame<-data.frame(EstY=EstY,EstT=EstT,EstP=EstP,EstTsq=EstT^2,EstPsq=EstP^2,Intr=EstT*EstP)

TestT<-t[testindex]
TestY<-y[testindex]
TestP<-propensities[testindex]
Testcov<-cov[testindex,]
TestFrame<-data.frame(EstY=TestY,EstT=TestT,EstP=TestP,EstTsq=TestT^2,EstPsq=TestP^2,Intr=TestT*TestP)
predsoutcomeGBM<-data.frame(t=seq(0,20000,100),pred=0)


modelGBM<-gbm(EstY~.,data=EstFrame,distribution = "gaussian",n.trees = 15000,interaction.depth = 2,cv.folds = 10)
besttrees<-gbm.perf(modelGBM,method = "cv")
predsoutcomeGBM<-data.frame(t=seq(0,20000,100),pred=0)
for(treat in seq(0,20000,100)){
  TestP<-fitted(npcdens(exdat=Testcov, eydat=rep(treat,nrow(Testcov)), bws=bws,txdat=cov,tydat=t))
  TestT<-rep(treat,nrow(Testcov))
  EstFrame<-data.frame(EstT=TestT,EstP=TestP,EstTsq=TestT^2,EstPsq=TestP^2,Intr=TestT*TestP)
  estResp<-mean(predict(modelGBM,newdata=EstFrame,num.trees=besttrees))
  predsoutcomeGBM[predsoutcomeGBM$t==treat,]$pred<-estResp
}

ggplot()+
  geom_point(aes(x=seq(0,20000,100),y=predsoutcome$pred,colour="Flexible linear regression"))+
  geom_point(aes(x=seq(0,20000,100),y=predsoutcomeGBM$pred,colour="GBM"))+
  xlab("t")+ylab("Contributions (thousands)")

