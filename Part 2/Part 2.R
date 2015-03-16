install.packages("UsingR")
library(UsingR)
?father.son
mean(father.son$sheight)
sub71 = subset(father.son, 70.5<=father.son$fheight&father.son$fheight<71.5)
mean(sub71$sheight)
sub71
father.son$fheight
X= matrix(1:1000,100,10)
X[25,3]
x = 1:10
x2 =  2*x
x3 =  3*x
x4 =  4*x
x5 =  5*x
a= cbind(x,x2,x3,x4,x5)
sum(a[7,])
matrix(1:60,20,3,byrow=TRUE)
x=11:20;rbind(x,2*x,3*x)
x=1:40;
X=matrix(3*x,20,2)
X = matrix(c(3,2,1,5,4,2,-1,0,-5,2,5,0,1,-1,-5,1),4,4)
Y = matrix(c(10,5,7,4),4,1)
solve(X)%*%Y
a <- matrix(1:12, nrow=4)
b <- matrix(1:15, nrow=3)
a%*%b
sum(a[3,]*b[,2])



g = 9.8
n = 25
tt = seq(0,3.4,len=n)
f = 56.67+0*tt-0.5*g*tt^2
y = f+rnorm(n,sd=1)
plot(tt,y,xlab="Time in secs",ylab="Distance in meters")
lines(tt,f, col=2)
rss= function(Beta0,Beta1, Beta2){
  r = y-(Beta0+Beta1*tt+Beta2*tt^2)
  sum(r^2)
}
Beta2s = seq(-10,0,len=100)
RSS = sapply(Beta2s, rss,Beta0=45,Beta1=0)
plot(Beta2s, RSS,type = "l")
tt2 =tt^2
fit = lm(y~tt+tt2)
X =cbind(rep(1,length(tt)),tt,tt^2)
Beta=matrix(c(55,0,5),3,1)
r = y-X %*% Beta
RSS = t(r)%*%r
RSS = crossprod(r)
betahat = solve(t(X)%*%X)%*%t(X)%*%y
QR = qr(X)
Q =qr.Q(QR)
R=qr.R(QR)
backsolve(R, crossprod(Q,y))


g = 9.8 ## meters per second
h0 = 56.67
v0 = 0
n = 25
tt = seq(0,3.4,len=n) ##time in secs, t is a base function
#y = h0 + v0 *tt - 0.5* g*tt^2 + rnorm(n,sd=1)
X = cbind(1,tt,tt^2)
A = solve(crossprod(X))%*%t(X)
Ns = 100000
gs = sapply(1:Ns,function(j){y = h0 + v0 *tt - 0.5* g*tt^2 + rnorm(n,sd=1)})
sd(-2*(A%*%gs)[3,])



library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)
N = 50
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat = lm(y~x)$coef
slopes = sapply(1:10000, function(j){
  N = 50
  index = sample(n,N)
  sampledat = father.son[index,]
  x = sampledat$fheight
  y = sampledat$sheight
  betahat = lm(y~x)$coef
})
sd(slopes[2,])
cov(x,y)



library(UsingR)
x = father.son$fheight
y = father.son$sheight
n = length(y)
N = 50
set.seed(1)
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat = lm(y~x)$coef
fit = lm(y ~ x)
fit$fitted.values
SSR =sum(fit$residuals^2)
sigma2 = SSR/48
X = cbind(rep(1,N), x)
sqrt(diag(solve(t(X)%*%X))*sigma2)
