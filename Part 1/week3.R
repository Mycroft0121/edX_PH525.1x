1-pnorm(2)+pnorm(-2)
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- "mice_pheno.csv"
if (!file.exists(filename)) download(url, destfile=filename)
dat <-read.csv(filename)
controlPopulation <- dat[dat$Sex=="F" & dat$Diet=="chow",3]
length(controlPopulation)
hfPopulation <- dat[dat$Sex=="F"&dat$Diet=="hf",3]
length(hfPopulation)
library(rafalib)
mypar2(1,2)
hist(hfPopulation)
hist(controlPopulation)
mypar2(1,2)
qqnorm(hfPopulation);qqline(hfPopulation)
qqnorm(controlPopulation);qqline(controlPopulation)
pops <- read.csv("mice_pheno.csv")
head(pops)
hf <- pops[pops$Diet=="hf"&pops$Sex=="F",3]
hist(hf)

chow <- pops[pops$Diet=="chow"&pops$Sex=="F",3]
hist(chow)
mean(hf)-mean(chow)

x <-sample(hf,12)
y <-sample(chow,12)
mean(x)-mean(y)
Ns <- c(3,5,10,25)
B <- 10000
res <- sapply(Ns, function(n){
  sapply(1:B, function(j){
    mean(sample(hf,n))-mean(sample(chow,n))
  })
})

library(rafalib)
mypar2(2,2)
for(i in seq(along=Ns)){
  title <- paste("Avg=", signif(mean(res[,i]),3))
  title <- paste(title,"SD=",signif(sd(res[,i]),3))
  qqnorm(res[,i], main=title)
  qqline(res[,i])
}

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
if (!file.exists(filename)) download(url,filename)
dat <- read.csv(filename)
control <- dat[1:12,2]
treatment <- dat[13:24,2]
diff <- mean(treatment)-mean(control)
sd(control)
sd(control)/sqrt(length(control))
se <- sqrt(var(treatment)/length(treatment)+var(control)/length(control))
tstat <-diff/se
1-pnorm(tstat)+pnorm(-tstat)
qqnorm(treatment)
qqline(treatment)
qqnorm(control)
qqline(control)
t.test(treatment, control)

babies = read.table("babies.txt", header=TRUE)
bwt.nonsmoke = babies$bwt[babies$smoke==0]
bwt.smoke = babies$bwt[babies$smoke==1]
mean(bwt.nonsmoke)-mean(bwt.smoke)
sd(bwt.nonsmoke)
sd(bwt.smoke)
dat.ns =  bwt.nonsmoke[1:30]
dat.s =  bwt.smoke[1:30]
X.ns = mean(dat.ns)
sd.ns=sd(dat.ns)
X.s = mean(dat.s)
sd.s = sd(dat.s)
sd.diff = sqrt(sd.ns^2/30+sd.s^2/30)
tval = (X.ns - X.s)/sd.diff
t.test(dat.ns,dat.s)$statistic
pval = 1 - pnorm(abs(tval)) + pnorm(-abs(tval))
pval


library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- "mice_pheno.csv"
if (!file.exists(filename)) download(url,destfile=filename)
dat <- read.csv(filename)
pop <- dat[dat$Sex=="F"& dat$Diet=="chow",3]
mu <- mean(pop)
N<-30
y<-sample(pop, N)
mean(y)
se <- sd(y)/sqrt(N)
se
Q <-qnorm(1-0.05/2)
interval <- c(mean(y)-Q*se, mean(y)+Q*se)
plot(mu+c(-7,7),c(1,1), type="n", xlab="weights", 
     ylab="intervals", ylim=c(1,100))
abline(v=mean(pop))
lines(interval, c(1,1))
for(i in 2:100){
  y <- sample(pop, N)
  se <- sd(y)/sqrt(N)
  interval <- c(mean(y)-Q*se, mean(y)+Q*se)
  color <- ifelse(interval[1]<=mean(pop)&
                    interval[2]>=mean(pop),1,2)
  lines(interval, c(i,i), col=color)
}


N<-5
y<-sample(pop, N)
mean(y)
se <- sd(y)/sqrt(N)
se
Q <-qt(1-0.05/2,4)
interval <- c(mean(y)-Q*se, mean(y)+Q*se)
plot(mu+c(-7,7),c(1,1), type="n", xlab="weights", 
     ylab="intervals", ylim=c(1,100))
abline(v=mean(pop))
for(i in 1:100){
  y <- sample(pop, N)
  se <- sd(y)/sqrt(N)
  interval <- c(mean(y)-Q*se, mean(y)+Q*se)
  color <- ifelse(interval[1]<=mean(pop)&
                    interval[2]>=mean(pop),1,2)
  lines(interval, c(i,i), col=color)
}
babies = read.table("babies.txt", header=TRUE)
bwt.nonsmoke = babies$bwt[babies$smoke==0]
bwt.smoke = babies$bwt[babies$smoke==1]


N<-30
nonsmoke= sample(bwt.nonsmoke,N)
smoke= sample(bwt.smoke,N)
CIs = replicate(1000, t.test(sample(bwt.nonsmoke, 30), sample(bwt.smoke, 30))$conf.int)
mean(CIs[2,] - CIs[1,])
popdiff = mean(bwt.nonsmoke) - mean(bwt.smoke)
mean(CIs[1,] < popdiff & CIs[2,] > popdiff)
dat.ns = sample(bwt.nonsmoke, 30)
dat.s = sample(bwt.smoke, 30)
X.ns = mean(dat.ns)
sd.ns = sd(dat.ns)
X.s = mean(dat.s)
sd.s = sd(dat.s)
sd.diff = sqrt(sd.ns^2/30 + sd.s^2/30)
tval = (X.ns - X.s)/sd.diff
qnorm(1-0.05/2)
ci.upper = (X.ns-X.s) + sd.diff*1.96
ci.lower = (X.ns-X.s) - sd.diff*1.96

babies = read.table("babies.txt", header=TRUE)
bwt.nonsmoke = babies$bwt[babies$smoke==0]
bwt.smoke = babies$bwt[babies$smoke==1]
N<-15
alpha<-0.01
B<-1000
rejections <- sapply(1:B,function(i){
  nonsmoke <-sample(bwt.nonsmoke,N)
  smoke <-sample(bwt.smoke,N)
  t.test(nonsmoke,smoke)$p.value<alpha
})
mean(rejections)
d = read.csv("assoctest.csv")
chisq.test(table(d))
fisher.test(table(d))


#inference III
library(downloader)
url<-"https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- tempfile()
download(url,destfile="babies.txt")
dat <- read.table("babies.txt",header=TRUE)
smokers <- sample(dat$bwt[dat$smoke==1],10)
nonsmokers <- sample(dat$bwt[dat$smoke==0],10)
mean(smokers)-mean(nonsmokers)
babies = read.table("babies.txt", header=TRUE)
bwt.nonsmoke = babies$bwt[babies$smoke==0]
pop.var = var(bwt.nonsmoke)
vars = replicate(1000, var(sample(dat$bwt[dat$smoke==0],50)))
hist(vars)
mean(vars>(1.5*pop.var))
sample.size = 2:400
var.estimate = sapply(sample.size, function(n) var(sample(bwt.nonsmoke, n)))
plot(sample.size, var.estimate)
abline(h=pop.var, col="blue")
set.seed(0)
N <- 50
smokers <- sample(babies$bwt[babies$smoke==1], N)
nonsmokers <- sample(babies$bwt[babies$smoke==0], N)
obs <- median(smokers) - median(nonsmokers)
mediandiff <- replicate(1000, {
  all <- sample(c(smokers,nonsmokers))
  smokersstar <- all[1:N]
  nonsmokersstar <- all[(N+1):(2*N)]
  return(median(smokersstar) - median(nonsmokersstar))
})
median(abs(avgdiff) > abs(obs))
