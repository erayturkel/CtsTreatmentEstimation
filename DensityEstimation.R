
dfkern<-data.frame(X[,1:p], t=t)
n <- names(dfkern)
f <- as.formula(paste("t ~", paste(n[!n %in% "t"], collapse = " + ")))
bandwidths<-npcdensbw(f, data=dfkern,bwmethod= 'normal-reference',memfac=200)
DensEstimate<-npcdens(bws=bandwidths,formula=f,data=dfkern, exdat= X[,(1:p)] , eydat=t)
plot(DensEstimate)
