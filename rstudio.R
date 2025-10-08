6+9
5+22
6 + 9
x- 1
x-1
x <- 15
x- 1
x + 4 <- 15
x = 5
5*x-> x
x
x <- c(1,-1, 3.5, 2)
x
x + 2
x + 2
x^2
sum((x- mean(x))^2)
1:10
-3:4
9:5
seq(from=2,to=6,by=0.4)
seq(from=-1,to=1,length=6)
rep(5,3)
rep(2:5,each=3)
rep(-1:3,length.out=10)
2^(0:10)
1:3+rep(seq(from=0,by=10,to=30),each=3)
1:10*c(-1,1)
1:7*1:2
1:7*1:2
x <- c(5,9,2,14,-4)
x[3]
x[c(2,3,5)]
x[1:3]
x[3:length(x)]
x > 4
x[x > 4]
x[-1]
x[-c(1,4)]
x <= 2
x == 2
x != 2
(x > 0) & (x < 10)
(x == 5) | (x > 10)
!(x > 5)
x <- c("Hello", "how do you do", "lovely to meet you", 42)
x
x[2:3]
x[-4]
c(x[1:2], "goodbye")
matrix(1:12, nrow=3, ncol=4)
matrix(1:12, nrow=3)
matrix(1:3, nrow=3, ncol=4)
matrix(1:12,nrow=3,byrow=TRUE)
diag(3)
diag(1:3)
1:5%o%1:5
outer(1:3,1:4,"+")
A<-matrix(c(1:8,10),3,3)
x<-c(1,2,3)
A%*%x
A*x
t(A)
det(A)
diag(A)
function (a, b, ...)
  solve(A)
diag(A)
solve(A)
A[2,1]
A[2,2:ncol(A)]
A[,1:2]
A[c(),1:2]
A[2,2:ncol(A),drop=FALSE]
cbind(A,t(A))
rbind(A,1,0)
x <- list(1:3, TRUE, "Hello", list(1:2, 5))
x[[3]]
x[c(1,3)]
x <- list(y=1:3, TRUE, z="Hello")
x
x$y
x[[1]]
names(x)
library(MASS)
head(hills)
class(hills)
is(hills, "data.frame")
hills[3,]
hills[hills$dist >= 12,]
hills$time
View(hills)
hills[1,]
hills[3]
hills %*% c(1,2,4)
mean(hills)
plot(hills$climb[hills$dist < 10], hills$time[hills$dist < 10])
with(hills, plot(climb[dist < 10], time[dist < 10]))
books<-data.frame(author=c("Ripley", "Cox", "Snijders", "Cox"),
                  year=c(1980, 1979, 1999, 2006),
                  publisher=c("Wiley", "Chapman", "Sage", "CUP"))
books
books <- data.frame(author=c("Ripley", "Cox", "Snijders", "Cox"),
                    year=c(1980, 1979, 1999, 2006),
                    publisher=c("Wiley", "Chapman", "Sage", "CUP"))
books
set.seed(1442)
height = round(rnorm(100, mean=rep(c(170,160),each=50), sd=10))
sex = rep(c("M", "F"), each=50)
head(sex)
Sex = as.factor(sex)
head(Sex)
plot(Sex, height)
as.integer(Sex)
attributes(Sex)
as.character(Sex)
Race = factor(birthwt$race)
names(hills)
row.names(hills)
attributes(hills)
attributes(hills) <- c(attributes(hills), list(type="races"))
attributes(hills)
dat <- read.table("smoking.dat", header=True)
head(dat)
dat <- read.table("smoking.dat", header=TRUE)
setdiff
function (x, y)
{
}
x <- as.vector(x)
{
  x <- as.vector(x)
  y <- as.vector(y)
  unique(if (length(x) || length(y))
    x[match(x, y, 0L) == 0L]
    else x)
}
unique(if (length(x) || length(y))
{
  x <- as.vector(x)
  y <- as.vector(y)
  unique(if (length(x) || length(y))
    x[match(x, y, 0L) == 0L]
    else x)
}
args(setdiff)
a<-c(1,4,5,7)
b<-c(1,2,5,9)
setdiff(a,b)
setdiff(a,b)
setdiff(b,a)
setdiff(y=b,a)
x<-rnorm(10)
y<-x+rnorm(10)
lm(y~x)
args(lm)
lm(form = y ~ x)
lm(formula = y ~ x)
square = function(x) {
  x^2
}
square(4)
mean2 <- function(x) {
  n <- length(x)
  sum(x)/n
}
mean2(1:10)
x <- mean2(1:10)
x
factorial2 = function(n) {
  out = 1
  for (i in 1:n) {
    +
      out = out*i
  }
  out
}
factorial2function(n){
  for (i in 1:n) {
    out = out*i
  }
  factorial2 = function(n) {
    out = 1
    for (i in 1:n) {
      out = out*i
    }
    out
  }
  factorial2(10)
  for (sillyname in 1:4) print(sillyname)
  sillyname
  n = 0
  1:n
  seq(n)
  seq_len(n)
  for (i in seq_len(n)) print(i)
  for (i in seq(n)) print(i)
  abs2 = function(x) {
    if (x < 0) out =-x
    else out = x
    out
  }
  abs2(-4)
  abs2(c(1,-3))
  ifelse(TRUE, 94, "hello")
  ifelse(FALSE, 94, "hello")
  isPrime = function(n) {
    i = 2
    if (n < 2) return(FALSE)
    while (i < sqrt(n)) {
      if (n %% i == 0) return(FALSE)
      i = i+1
    }
    return(TRUE)
  }
  isPrime(10)
  isPrime(10)
  isPrime(37)
  system.time(for (i in 1:1e6) i^2)
  system.time(seq_len(1e6)^2)
  out[i] = sum(A[i,]*b)
  mult2 = function(A, b) {
    n1 = nrow(A)
    n2 = ncol(A)
    out = numeric(n1)
    for (i in 1:n1) {
      out[i] = sum(A[i,]*b)
    }
    out
 }
  A = matrix(rnorm(1e6), 1e3, 1e3)
  b = rnorm(1e3)
  system.time(mult(A,b))
  system.time(mult2(A,b))
  system.time(A %*% b)
  fib = function(n) {
    if (n < 2) return(1)
    else return(fib(n-1) + fib(n-2))
  }
  fib(10)
  x <- 3
  f=function(y) {
    x <- 5
    x+y
  }
  f(4)
  x
  x<-3
  g=function(y) {
    x+y
  }
  g(4)
  conv=1.609
  with(hills, mean(dist/conv))
  hist(nlschools$lang, breaks=25, col=2)
  x = cumsum(rnorm(250))
  plot(x, type="l", col=3)
  y = x + rnorm(300)
  x = rnorm(300)
  y = x + rnorm(300)
  plot(x,y, pch=20, col=4, cex=0.5)
  abline(a=0, b=1, lty=4, lwd=1.5)
  abline(lm(y~x),col=2)
  legend(x=-4,y=4,legend=c("y=x","lineofbestfit"),
         lty=c(4,1),lwd=c(1.5,1),col=1:2)
  x ~ a + b*c
  data(genotype)
  head(genotype)
  boxplot(Wt ~ Litter, data=genotype)
  library(lattice)
  head(crabs)
  histogram(~ height | voice.part, data=singer)
  library(MASS)
  densityplot(galaxies)
  plot(sin,-2*pi, 2*pi)
  func = function(x, n=10)
    for (i in 2:10) {
      x <- seq(from=-0.5, to=1.2, length=1000)
      y <- sapply(x, func, n=i)
      points(x, y, type="l")
    }
  func = function(x,y) (x^3-x)*sin(x+y)
  xs = ys = seq(-5, 5, length.out=100)
  out = outer(xs, ys, func)
  wireframe(out)
  x + 2
  x[-1]
  x[-c(1,4)]
  x <- c(5,9,2,14,-4)
  x[3]
  x[c(2,3,5)]
  x[1:3]
  x[3:length(x)]
  x > 4
  x[x > 4]
  x[-1]
  x[-c(1,4)]
  hist(nlschools$lang, breaks=25, col=2)
  hist(nlschools$lang, breaks=25, col=2, xlab="Score", main="Language test scores of Dutch 8th grade pupils")
  x = cumsum(rnorm(250))
  plot(x, type="l", col=3)
  x = rnorm(300)
  y = x + rnorm(300)
  plot(x,y, pch=20, col=4, cex=0.5)
  abline(a=0, b=1, lty=4, lwd=1.5)
  legend(x=-4, y=4, legend=c("y=x","line of best fit"), lty=c(4,1), lwd=c(1.5,1), col=1:2)
  x ~ a + b*c
  data(genotype)
  library(lattice)
  head(crabs)
  form <- log(FL) ~ log(RW) | sp*sex
  form
  xyplot(form, data=crabs)
  histogram(~ height | voice.part, data=singer)
  library(MASS)
  densityplot(galaxies)
  plot(sin, -2*pi, 2*pi)
  idx = 1:n
  plot(function(x) log(1+x), -0.5, 1.2)
  for (i in 2:10) {
    x <- seq(from=-0.5, to=1.2, length=1000)
    y <- sapply(x, func, n=i)
    points(x, y, type="l")
  }
  for (i in 2:10) {
    x <- seq(from=-0.5, to=1.2, length=1000)
    y <- sapply(x, func, n=i)
    points(x, y, type="l")
  }
  for (i in 2:10) {
    x <- seq(from=-0.5, to=1.2, length=1000)
    y <- sapply(x, func, n=i)
    points(x, y, type="l")
  }
  func = function(x,y) (x^3-x)*sin(x+y)
  xs = ys = seq(-5, 5, length.out=100)
  out = outer(xs, ys, func)
  wireframe(out)
  set.seed(1328) # to get the same values as me
  x = rnorm(100) # generate data> plot.new()
  plot.new()
  plot.window(xlim=c(-3,3), ylim=c(-0.1,0.5))
  axis(side=1, pos=-0.1)
  hist(x, breaks=15, add=TRUE, freq=FALSE, col=2)
  plot(dnorm, -3, 3, add=TRUE)
  points(x, rep(-0.05,100), pch="|")
  title(main="Normal random variables")
  pdf("plotfile.pdf")
  plot(hills)
  dev.off()
  A = cbind(1:10, (1:10)^2, (1:10)^3)
  apply(A, 2, sum)
  library(MASS)
  apply(hills, 2, mean)
  apply(hills, 2, sd)
  x = matrix(rnorm(4e6), 2000, 2000)
  system.time(apply(x,1,sum))
  system.time(rowSums(x))
  mu = c(-2,-1,0,1,2)
  out = lapply(mu, function(x) rnorm(100, mean=x))
  lapply(out, mean)
  sapply(out, mean)
  out = replicate(20, rnorm(100), simplify=FALSE)
  library(MASS)
  head(genotype)
  with(genotype, tapply(Wt, Mother, mean))
  with(genotype, tapply(Wt, Mother, summary))
  tapply(genotype$Wt, genotype[,1:2], mean)
  mapply(seq, from=c(1,4,-3), to=c(2,9,0), by=0.5)
  hist(nlschools$lang, breaks=25, col=2)
  hist(nlschools$lang, breaks=25, col=2, xlab="Score", main="Language test scores of Dutch 8th grade pupils")
  x = cumsum(rnorm(250))
  plot(x, type="l", col=3)
  data(genotype)
  library(lattice)
  head(crabs)
  form <- log(FL) ~ log(RW) | sp*sex
  form
  func = function(x,y) (x^3-x)*sin(x+y)
  xs = ys = seq(-5, 5, length.out=100)
  out = outer(xs, ys, func)
  wireframe(out)
  set.seed(1328) # to get the same values as 
  x = rnorm(100) # generate data
  plot.new()
  plot.window(xlim=c(-3,3), ylim=c(-0.1,0.5))
  axis(side=1, pos=-0.1)
  hist(x, breaks=15, add=TRUE, freq=FALSE, col=2)
  plot(dnorm, -3, 3, add=TRUE)
  points(x, rep(-0.05,100), pch="|")
  title(main="Normal random variables")
  pdf("plotfile.pdf")
  plot(hills)
  dev.off()
  out = replicate(20, rnorm(100), simplify=FALSE)
  library(MASS)
  head(genotype)
  with(genotype, tapply(Wt, Mother, mean))
  with(genotype, tapply(Wt, Mother, summary))
  tapply(genotype$Wt, genotype[,1:2], mean)
  mapply(seq, from=c(1,4,-3), to=c(2,9,0), by=0.5)
  occupationalStatus
  arr = array(1:18, dim=c(2,3,3))
  arr
  dim(arr)
  arr[1,2,3]
  arr[,2,]
  house=dget("housing.dat")
  margin.table(house, 1:2)
  library(MASS)
  head(cabbages)
  table(cabbages$Date)
  head(Nile)
  bins=c(0,seq(from=700,to=1300,by=100),Inf)
  bins
  disNile=cut(Nile,bins)
  head(disNile)
  table(disNile)
  cat("Hello")
  myFunc = function(x) {
    cat("Squaring the number\n")
    x^2
  }
  myFunc(3)
  cat("She said \"Hello\" to the man.\n")
  cat("This is a backslash: \\ - how does it make you feel?\n")
  cat(57, "clouds", "\n")
  cat(57, "clouds", "\n", sep="")
  withBox("Hello")
  paste("Hello", "there")
  paste("Plan", LETTERS[1:5])
  paste("x", 1:10, sep="")
  paste(LETTERS[1:10])
  paste(LETTERS[1:10], collapse=" ")
  f(5)
  listfunc(c(1,4,2))
  nchar("How long is this string?")
  strsplit("separate words are fun", split=" ")
  substr("I don't want all of this string", 14, 24)
  x <- rnorm(100)
  y <- x + rnorm(100)
  mod1 <- lm(y ~ x)
  summary(x)
  summary(mod1)
  methods(summary)
  class(mod1)
  summary.lm(mod1)
  attributes(hills)
  attributes(hills,"class")
  attr(hills,"row.names")[4]="Big Hill"
  rows.names(hills)
  attr(hills, "quality") = "Great!"
  attributes(hills)[4]
  x <- list(temp=19.5, wind=12, wind.dir="SSW", rain=20, summ="showers")
  class(x) <- "weather"
  attr(x, "class")
  print.weather = function(object) {
    cat("Weather report\n")
    cat("Temperature ", object$temp, "C\n", sep="")
    cat("Wind ", object$wind, " kts, from ", object$wind.dir, "\n", sep="")
    cat("Rain ", object$rain, " mm: ", object$summ, "\n", sep="")
    object
  }
  "a"+1
  set.seed(241)
  rndSteps <- rnorm(N)
  rndWalk <- cumsum(rndSteps)
  row_sums <- function(x, cols) {
    apply(x[,cols], 1, sum)
  }
  x <- matrix(1:9, 3, 3)
  row_sums(x, 2:3)
  row_sums(x,2)
  if(any(x > 3) && y != 2) {
    print("Hello")
  }
  any(x>3)
  y!=2
  y
  traceback()
  options(error=recover)
  0.3-0.1-0.2
  x=0.3-0.1-0.2
  x==0
  all.equal(x,0)
  if (isTRUE(all.equal(x, 0))) print("We've got nothing!")
  2^-1074
  2^-1075
  (2^-1074)/1.5
  2^2000
  Inf/2
  (-1)*Inf
  0*Inf
  Inf - Inf
  exp-(Inf)
  i1 = 0L
  i2 = 0
  (i1 == i2)
  identical(i1,i2)
  i3=1e40+1
  i3-1e40
  
  