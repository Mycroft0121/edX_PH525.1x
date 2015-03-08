dat <-read.csv("Project 1/femaleMiceWeights.csv")
dat[1:12,2]
mean(dat[13:24,2])-mean(dat[1:12,2])
s = split(dat[,2],dat[,1])
stripchart(s,vertical=TRUE,col=1:2)
abline(h=sapply(s,mean),col=1:2)
highfat = s[["hf"]]
sample(highfat,6)
sum(highfat>30)/12
population <-read.csv("Project 1/femaleControlsPopulation.csv")
control <- sample(population[,1],12)
mean(control)
n <- 10000
null <- vector("numeric",n)
for(i in 1:n){
  control <- sample(population[,1],12)
  treatment <- sample(population[,1],12)
  null[i] <- mean(treatment)-mean(control)
}
diff <- mean(dat[13:24,2])-mean(dat[1:12,2])
mean(null>diff)
sampleMean = replicate(10000,mean(sample(population[,1],12)))
head(sampleMean)
plot(sampleMean)
null = replicate(10000,mean(sample(population[,1],12)))-replicate(10000,mean(sample(population[,1],12)))
head(null)                                               
plot(null)
hist(null)
abline(v=diff,col="red")
abline(v=-diff,col="red")
mean((null>diff)|(-diff>null))
library(gapminder)
data(gapminder)
head(gapminder)
x = data.frame(gapminder$country,gapminder$lifeExp )
plot(x)
year_vector =  gapminder$year ==1952
newdf=x[year_vector,]
mean(newdf$gapminder.lifeExp<=40)
mean(newdf$gapminder.lifeExp<=60)-mean(newdf$gapminder.lifeExp<=40)
y = gapminder$pop 
y_new = log10(y[year_vector])
ave = mean(y_new)
stdv = sd(y_new)
qqnorm(y_new)
standarize = function(a,b,c){(a-b)/c}
z = standarize(y_new, ave, stdv)
qqnorm(z)
tail(sort(z),1)
F = function(q) pnorm(q, mean=mean(y_new), sd=sd(y_new))
(F(7)-F(6))*length(y_new)
n = length(y_new)
qs = (1:n-0.5)/n
qnorm(qs[1])
