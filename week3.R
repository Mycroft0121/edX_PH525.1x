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
