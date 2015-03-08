load("skew.RData")
par(mfrow=c(1,1))
for (i in 1:9) {
  qqnorm(dat[,i])
}

head(InsectSprays)
boxplot(split(InsectSprays$count, InsectSprays$spray))
install.packages("UsingR")
library(UsingR)
data(father.son)
plot(father.son$fheight, father.son$sheight)
cor(father.son$fheight, father.son$sheight)
identify(father.son$fheight, father.son$sheight)
x = father.son$fheight
y = father.son$sheight
n = nrow(father.son)
plot(scale(x), scale(y))
abline(h=0, v=0)
mean(scale(x)*scale(y))


install.packages("dplyr")
library(dplyr)
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- "msleep_ggplot2.csv"
if (!file.exists(filename)) download(url,filename)
msleep <- read.csv("msleep_ggplot2.csv")
head(msleep)
msleep$rem_proportion =  msleep$sleep_rem/msleep$sleep_total
head(msleep)
msleep %>% 
  mutate(rem_proportion = sleep_rem / sleep_total) %>% 
  # XX Group the animals by their taxonomic order
  group_by(order) %>% 
  # XX Summarise by the median REM proportion
  # Note: for right answer, can't use the option of na.rm=TRUE for the median
  summarise(median.rem = median(rem_proportion), total = n()) %>%
  # head
  # XX Arrange by the median REM proportion
  arrange(median.rem) %>%
  # Take the head() of this to see just the orders 
  #       with smallest median REM proportion
  head


data(ChickWeight)
plot(ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)
head(ChickWeight)
chick = reshape(ChickWeight,idvar=c("Chick","Diet"),timevar="Time",direction="wide")
head(chick)
chick = na.omit(chick)
mean(chick$weight.4)
length(chick$weight.4)
sd(c(chick$weight.4, 3000))/sd(chick$weight.4)
mad(c(chick$weight.4, 3000))/mad(chick$weight.4)
plot(chick$weight.4, chick$weight.21)
cor(c(chick$weight.4, 3000), c(chick$weight.21, 3000))/cor(chick$weight.4, chick$weight.21)
stripchart(chick$weight.4 ~ chick$Diet, method="jitter", vertical=TRUE)
x <- chick[chick$Diet ==1, "weight.4"]
y <- chick[chick$Diet ==4, "weight.4"]
t.test(c(x, 200), y)$p.value
par(mfrow=c(1,3))

boxplot(x,y)

boxplot(x,y+10)

boxplot(x,y+100)
t.test(x,y + 10)$statistic - t.test(x,y + 100)$statistic
wilcox.test(c(1,2,3),c(4,5,6))

wilcox.test(c(1,2,3),c(400,500,600))