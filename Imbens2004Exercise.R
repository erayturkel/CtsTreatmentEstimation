library(plotly)
library(hdrcde)
set.seed(1)
X1<-rnorm(mean=1,1000)
X2<-rnorm(mean=1,1000)


alpha=0.7
T.vals<-pmax(0,rnorm(1000)+alpha*X1)


X<-matrix(c(X1,X2),ncol=2)

Y=X1+X2+2*T.vals^2+rnorm(1000)

basedf<-data.frame(Y,X1,X2,T.vals)
linmodel<-lm(Y~X1+X2+T.vals*X2*X1,data=basedf)

sort(unique(T.vals))
y.naive.pred<-data.frame(T.vals=sort(unique(T.vals)),Y=0)
for(i in unique(T.vals)){
  preddf<-data.frame(X1,X2,T.vals=i)
  y.naive.pred[y.naive.pred$T.vals==i,]$Y<-mean(predict(linmodel,preddf))
}
plot(y.naive.pred$T.vals,y.naive.pred$Y)

t=seq(0,15,0.1)
Ey.givent<- 2+2*t*(1+ (alpha^2)*(t-alpha)/(1+alpha^2))

plot(T.vals,Y)
plot(X1,T.vals)
plot(X2,T.vals)


DensityEst<-cde(X,T.vals,x.margin=X,y.margin = T.vals)
plot(DensityEst)
DensityEst$x

plot_ly(x=X1, y=X2, z = T.vals, type = "contour") 
