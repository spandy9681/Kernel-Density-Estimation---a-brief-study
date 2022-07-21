#Section 1
#Histogram
#Data Driven
set.seed(020524)
x = c(rnorm(100,mean = 1,sd = 0.5),rnorm(100,mean = 5,sd = 1))
bins = 4*(1:9)

par(mfrow = c(3,3))
sq = seq(0,9,length.out = 1000)
M = max(dnorm(sq,mean = 1,sd = 0.5)+dnorm(sq,mean = 5,sd = 1))/2
for(i in bins)
{
  hist(x,breaks = i,probability = TRUE,main = bquote("No of bins = " ~ .(i)),col = "cyan",ylim = c(0,M+0.1))
  lines(sq,(dnorm(sq,mean = 1,sd = 0.5)+dnorm(sq,mean = 5,sd = 1))/2,col = "darkblue")
}

# Using cross validation to find optimal bandwidth 
par(mfrow = c(1,1))
dat = (x-min(x))/(max(x)-min(x))
par(mfrow = c(1,1))
hist(dat,probability = TRUE,breaks = 40,col = "cyan")
nbins = 1:36

Histo <- function(data,x,x0,h)
{
  n = length(data)
  m = floor((x-x0)/h)
  count = sum((data >= x0+m*h & data < x0+(m+1)*h))
  return(count/(n*h))
}

fnx <- function(x,h)
{
  return(Histo(data = dat,x = x,x0 = mean(dat),h = h))
}

M0h = NULL
n = length(dat)
for(i in 1:length(nbins))
{
  f_r_n = sapply(X = (0:nbins[i])/nbins[i], FUN = function(x){return((fnx(x,h = 1/nbins[i]))^2)})
  int_f_2 = sum(f_r_n)/nbins[i]
  v = NULL
  for(j in 1:length(dat))
  {
    v[j] = Histo(x = dat[j],data = dat[-j],x0 = mean(dat),h = 1/nbins[i])
  }
  M0h[i] = int_f_2 - 2*sum(v)/n
}

plot(nbins,M0h,type = "o",ylab = bquote(M[0](n[bins])),main = bquote("CV values for different choices of "~ n[bins]),xlab = bquote(n[bins]),pch = 20)
ind = which(M0h == min(M0h))
bin_opt = nbins[ind]
bin_opt
abline(v = bin_opt,lty = 2)
hist(x,breaks = bin_opt,prob = TRUE,col = "cyan")
lines(sq,(dnorm(sq,mean = 1,sd = 0.5)+dnorm(sq,mean = 5,sd = 1))/2,col = "darkblue")

# Section 2
# Kernel Density Estimation


# Defining different kernels
# Uniform Kernel

K_naive <- function(x)
{
  if(x >= -1 & x < 1)
  {return(1/2)}
  else
  {return(0)}
}


# Beta Kernel

K_beta <- function(x)
{
  if(x >= -1/2 & x < 1/2)
  {return(6*(x+1/2)*(1/2-x))}
  else
  {return(0)}
}

# Gaussian kernel

K <- function(x){return(dnorm(x))}

# Exponential Kernel

K <- function(x){return(dexp(x))}

# Epanechnikov Kernel

K_epnch <- function(x)
{
  if(abs(x) <= 1)
  {return((3/4)*(1-x^2))}
  else
  {return(0)}
}

# Type II

K_opt <- function(x)
{
  if(abs(x) <= sqrt(5))
  {return((3/(4*sqrt(5)))*(1-x^2/5))}
  else
  {return(0)}
}

# Triangular

K_tri <- function(x)
{
  if(abs(x) <= 1)
  {return(1-abs(x))}
  else
  {return(0)}
}

# Cosine Kernel
K_cosine <- function(x)
{
  if(abs(x) <= 1)
  {return((pi/4)*cos((pi*x)/2))}
  else
  {return(0)}
}


# Estimating density at a single point

f_kernel <- function(t,h,dta,K = K_opt)
{
  n = length(dta)
  t = ((n*h)^(-1))*sum(sapply((t - dta)/h , FUN = function(s){return(K(s))}))
  return(t)
}

# Estimating density at vector of points

f_kernel_vec <- function(t_vec,h,dta,K = K_opt)
{
  val = sapply(t_vec, FUN = function(t){return(f_kernel(t,h,dta,K))})
  return(val)
}

# Subsection 1

# Demonstrating how kernel density estimation works
par(mfrow = c(1,1))
plot.new()
n = 5;h = 0.5
dat = rnorm(n)
#hist(dat,prob = TRUE,ylim = c(0,1),breaks = 6,xlim = c(-3,3),col = "red")
plot(rg,f_kernel_vec(t_vec = rg,h = h,dta = dat,K = K_gaussian),col = "red",type = "l",ylab ="Density",main = "Estimated Density Function",lwd = 2,xlab = "Observed Samples")
points(dat,y = rep(0,length(dat)),pch = 20,xlim = c(-3,3),ylim = c(0,1))

for(i in 1:length(dat))
{
  sq = seq(dat[i]-2,dat[i]+2,length.out = 10^2)
  lines(sq,sapply(sq,FUN = function(x){return(((n*h)^-1)*K_gaussian((x-dat[i])/h))}),lty = 2,col = "darkblue")
  segments(x0 = dat[i],y0 = 0,x1 = dat[i],((n*h)^-1)*K_gaussian(0),lty = 2)
}



# Demonstrating how estimation depends on choice of different kernels
x = rnorm(100)
r = range(x)
sq = seq(r[1]-2,r[2]+2,length.out = 10^4)

#Drawing true and estimated densities
par(mfrow = c(2,3))

#Drawing the kernels
plot(seq(-1.5,1.5,length.out = 30),sapply(seq(-1.5,1.5,length.out = 30), FUN = function(x){return(K_unif(x))}),type = "l",main = "Uniform Kernel",ylim = c(0,1),xlab = "",ylab = "",col = "orange")
plot(seq(-3,3,length.out = 100),K_gaussian(seq(-3,3,length.out = 100)),type = "l",main = "Gaussian Kernel",xlab = "",ylab = "",col = "blue")
plot(seq(-3,3,length.out = 100),sapply(seq(-3,3,length.out = 100), FUN = function(x){return(K_opt(x))}),type = "l",main = "Epanenechnikov Kernel",xlab = "",ylab = "",col = "green")

#Drawing the KDEs
M = max(dnorm(sq))
plot(sq,f_kernel_vec(t_vec = sq,h = 0.422,dta = x,K = K_unif),col = "red",lwd = 1,type = "l",ylab = "Estimated Density",ylim = c(0,M+0.1))
lines(sq,dnorm(sq),type = "l",lwd = 1.5,col = "darkorchid1")

plot(sq,f_kernel_vec(t_vec = sq,h = 0.422,dta = x,K = K_gaussian),col = "red",lwd = 1,type = "l",ylab = "Estimated Density",ylim = c(0,M+0.1))
lines(sq,dnorm(sq),type = "l",col = "darkorchid1",lwd = 1.5)

plot(sq,f_kernel_vec(t_vec = sq,h = 0.422,dta = x,K = K_opt),col = "red",lwd = 1,type = "l",ylab = "Estimated Density",ylim = c(0,M+0.1))
lines(sq,dnorm(sq),type = "l",col = "darkorchid1",lwd = 1.5)

#Drawing the kernels
plot(seq(-1,1,length.out = 100),sapply(seq(-1,1,length.out = 100), FUN = function(x){return(K_beta(x))}),type = "l",main = "Beta Kernel",xlab = "",ylab = "",col = "orange")
plot(seq(-1.5,1.5,length.out = 100),sapply(seq(-1.5,1.5,length.out = 100), FUN = function(x){return(K_tri(x))}),type = "l",main = "Triangular Kernel",xlab = "",ylab = "",col = "green")
plot(seq(-1.5,1.5,length.out = 100),sapply(seq(-1.5,1.5,length.out = 100), FUN = function(x){return(K_cosine(x))}),type = "l",main = "Cosine Kernel",xlab = "",ylab = "",col = "blue")


#Drawing the KDEs
plot(sq,f_kernel_vec(t_vec = sq,h = 0.422,dta = x,K = K_beta),col = "red",lwd = 1,type = "l",ylab = "Estimated Density",ylim = c(0,M+0.1))
lines(sq,dnorm(sq),type = "l",lwd = 1.5,col = "darkorchid1")
plot(sq,f_kernel_vec(t_vec = sq,h = 0.422,dta = x,K = K_tri),col = "red",lwd = 1,type = "l",ylab = "Estimated Density",ylim = c(0,M+0.1))
lines(sq,dnorm(sq),type = "l",lwd = 1.5,col = "darkorchid1")
plot(sq,f_kernel_vec(t_vec = sq,h =0.422,dta = x,K = K_cosine),col = "red",lwd = 1,type = "l",ylab = "Estimated Density",ylim = c(0,M+0.1))
lines(sq,dnorm(sq),type = "l",lwd = 1.5,col = "darkorchid1")

####

# How estimated density looks if we use kernels which are asymmetric 
par(mfrow = c(1,1))
plot(sq,f_kernel_vec(t_vec = sq,h = 0.7,dta = x,K = K_expo),col = "red",lwd = 1,type = "l",ylab = "Density",ylim = c(0,M+0.1),main = "Estimated Density using Exponential Kernel")
lines(sq,dnorm(sq),type = "l",lwd = 1.5,col = "darkorchid1")
#Since the exponential kernel is positive valued for the values of x that
#lie on the left side of the observations hence in the first part its underestimated
#and on the right side its overestimated due to presence of many observations that are contributing


# Subsection 2
# Choosing the optimal bandwidth for density estimation
# Subjective Choice
library(ProbBayes)
dat = buffalo_jan$JAN
dat = dat[-1]
r = range(dat)
sq = seq(r[1]-0.1,r[2]+0.1,length.out = 10^4)
hist(dat,breaks = 15,probability = TRUE)
points(dat,rep(0,length(dat)))
# For h=5
lines(sq,f_kernel_vec(t_vec = sq,h = 5,dta = dat,K = K_opt),col = "red",lwd = 1,type = "l",ylab = "Density",main = "Estimated Density For ")

# For h=12
plot(sq,f_kernel_vec(t_vec = sq,h = 12,dta = dat,K = K_opt),col = "red",lwd = 1,type = "l",ylab = "Density",main = "Estimated Density For ")

# A much more interesting example
dat1 = rnorm(n = 500,mean = 0,sd = 1)
dat2 = rnorm(n = 100,mean = -1,sd = 1/10)
dat3 = rnorm(n = 100,mean = -1/2,sd = 1/10)
dat4 = rnorm(n = 100,mean = 0,sd = 1/10)
dat5 = rnorm(n = 100,mean = 1/2,sd = 1/10)
dat6 = rnorm(n = 100,mean = 1,sd = 1/10)
dat = c(dat1,dat2,dat3,dat4,dat5,dat6)
r = range(dat)
sq = seq(-3,3,length.out = 10^4)

den <- function(x)
{
  return(sapply(x,FUN = function(x){return((1/2)*dnorm(x) + (dnorm(x,-1,0.1)+dnorm(x,-0.5,0.1)+dnorm(x,0,0.1)+dnorm(x,0.5,0.1)+dnorm(x,1,0.1))/10)}))
}
plot(sq,den(sq),col = "green",main = "Estimated Density",type = "l",lwd = 1,ylab = "Density",lty = 2)
# For h=1
lines(sq,f_kernel_vec(t_vec = sq,h = 1,dta = dat,K = K_opt),lwd = 1,type = "l",col = "coral")
legend("topright",legend = c("h = 1"),fill = c("coral"))
# For h=0.1
plot(sq,den(sq),col = "green",main = "Estimated Density",type = "l",lwd = 1,ylab = "Density",lty = 2)
lines(sq,f_kernel_vec(t_vec = sq,h = 0.1,dta = dat,K = K_opt),col = "blueviolet",lwd = 1,type = "l")
legend("topright",legend = c("h = 0.1"),fill = c("blueviolet"))
# For h=0.05
plot(sq,den(sq),col = "green",main = "Estimated Density",type = "l",lwd = 1,ylab = "Density",lty = 2)
lines(sq,f_kernel_vec(t_vec = sq,h = 0.05,dta = dat,K = K_opt),col = "blue",lwd = 1,type = "l")
legend("topright",legend = c("h = 0.05"),fill = c("blue"))


# Demonstrating how important its in case of estimating exponential type 
# densities
n = 10^3
x = rexp(n)
r = range(x)

CV_h=function(n,h,K=K_opt,data)
{
  n = length(data)
  p=0
  for(i in 1:n)
  {
    p=p+log(f_kernel(t = data[i],K = K,h = h,dta = data[-i]))
  }
  return(p/n)
}

h1 = seq(.02,1,length.out=20)
data = x
n = length(x)
cv = sapply(h1,function(x){return(CV_h(n = n,h = x,K = K_opt,data = data))})
plot(h1,cv,type="l",xlab = "Binwidth h",ylab = "CV(h)")
h_opt = h1[which(cv==max(cv))]
abline(v = h_opt,lty = 2)

sq1 = seq(0,r[2]+2,length.out = 10^2)
fx = dexp(sq1)
plot(sq1,fx,type = "l",ylim = c(0,0.6),xlim = c(-2,r[2]+2),col = "violet",lwd = 2,main = bquote("Estimated density using optimal h"),ylab = "density",xlab = "x")
abline(v = 0,lty = 2)
sq2 = seq(-2,r[2]+2,length.out = 10^3)
den_est = f_kernel_vec(t_vec = sq2,h = h_opt,dta = x,K = K_opt)
lines(sq2,den_est,col = "red",lwd = 1,lty = 2)

#h = 0.03
den_est = f_kernel_vec(t_vec = sq2,h = 0.03,dta = x,K = K_opt)
lines(sq2,den_est,col = "blue",lwd = 0.5,lty = 2)

#h = 0.5
den_est = f_kernel_vec(t_vec = sq2,h = 0.5,dta = x,K = K_opt)
lines(sq2,den_est,col = "violet",lwd = 0.5,lty = 2)

#h = 1
den_est = f_kernel_vec(t_vec = sq2,h = 1,dta = x,K = K_opt)
lines(sq2,den_est,col = "violet",lwd = 1,lty = 2)

legend("topright",fill = c("red","violet"),legend = c("h = h_opt","h = 1"))
#What if we didn't take the optimal value of h

# Reference to a standard distribution

# Another natural approach is to consider a standard family of distribution
# and then evaluate the exact expression for h_opt

dat = rnorm(500,mean = 2,sd = 2.1)
r = range(dat)
sq = seq(r[1]-0.1,r[2]+0.1,length.out = 10^4)
sg = sd(dat)
n = length(dat)
h_opt = 1.06*sg*(n)^(-1/5)
h_opt
M = max(dnorm(sq,mean = 2,sd = 2.1))

#Density estimate for optimal choice of h
plot(sq,dnorm(sq,mean = 2,sd = 2.1),type = "l",lwd = 1,ylim = c(0,M+0.01),lty = 2,main = "Estimated Normal Density using Optimal Bandwidth",ylab = "Density")
lines(sq,f_kernel_vec(t_vec = sq,h = h_opt,dta = dat,K = K_opt),lwd = 1,type = "l",col = "coral")
plot(sq,dnorm(sq,mean = 2,sd = 2.1),type = "l",lwd = 1,ylim = c(0,M+0.01),lty = 2,main = "Estimated Normal Density using h = 0.2",ylab = "Density")
lines(sq,f_kernel_vec(t_vec = sq,h = 0.2,dta = dat,K = K_opt),lwd = 1,type = "l",col = "blue")
plot(sq,dnorm(sq,mean = 2,sd = 2.1),type = "l",lwd = 1,ylim = c(0,M+0.01),lty = 2,main = "Estimated Normal Density using h = 1",ylab = "Density")
lines(sq,f_kernel_vec(t_vec = sq,h = 1,dta = dat,K = K_opt),lwd = 1,type = "l",col = "pink")

# Cross Validation
CV_h=function(n,h,K=K_unif,data)
{
  n = length(data)
  p=0
  for(i in 1:n)
  {
    p=p+log(f_kernel(t = data[i],K = K,h = h,dta = data[-i]))
  }
  return(p/n)
}

h1=seq(.3,2,length.out=1000)
n=50
data=rnorm(n)
cv=sapply(h1,function(x){return(CV_h(n = n,h = x,K = K_opt,data = data))})
plot(h1,cv,type="l",xlab = "Binwidth h",ylab = "CV(h)")
h_opt=h1[which(cv==max(cv))]
abline(v = h_opt,lty = 2)

# Density Plot
r = range(data)
sq = seq(r[1]-0.1,r[2]+0.1,length.out = 10^4)
M = max(dnorm(sq,mean = 0,sd = 1))
#Density estimate for optimal choice of h
plot(sq,dnorm(sq,mean = 0,sd = 1),type = "l",lwd = 2,ylim = c(0,M+0.1))
lines(sq,f_kernel_vec(t_vec = sq,h = h_opt,dta = data,K = K_opt),lwd = 1,type = "l",col = "coral")
points(data,rep(0,n),pch = 20,col = rgb(0,0,1,alpha = 0.1))

## Put the optimal choice of binwidth thing for histogram here as it will then
## cover both the cross validation methods

# Subsection 3
# Large Sample Properties 

# Supremum between fn(x) & f(x)
z2=seq(0,1,length.out = 100)
z1=dbeta(z2,3,4)

Vn = function(n,K=K_opt,h,opt = TRUE)
{
  z=rbeta(n,3,4)
  y=f_kernel_vec(t_vec = z2,h = n^(-1/5),dta = z,K = K)
  return(max(abs(z1-y)))
}

n=seq(1,10^3,length.out = 100)
V=sapply(n,function(t){Vn(n = t,K = K_opt,opt = TRUE)})
plot(n,V,type = "o",pch = 20,ylab = bquote(V[n](x) ~ "= sup" ~ "|" ~ f[n](x)-f(x) ~"|"),main = bquote("Values of " ~ V[n](x) ~ "for different values of n"),ylim = c(-0.1,2))
abline(h = 0,lty = 2)

V

CV_h=function(n,h,K=K_unif,data)
{
  n = length(data)
  p=0
  for(i in 1:n)
  {
    p=p+log(f_kernel(t = data[i],K = K,h = h,dta = data[-i]))
  }
  return(p/n)
}

h1=seq(.001,0.3,length.out=20)
data=rbeta(n,3,4)
n=length(data)
cv=sapply(h1,function(x){return(CV_h(n = n,h = x,K = K_opt,data = data))})
plot(h1,cv,type="l",xlab = "Binwidth h",ylab = "CV(h)")
h_opt=h1[which(cv==max(cv))]
abline(v = h_opt,lty = 2)


Vn_draw = function(n,K=K_opt,h,...)
{
  z2=seq(0,1,by=.001)
  f_x=dbeta(z2,3,4)
  z=rbeta(n,3,4)
  fn_x=f_kernel_vec(t_vec = z2,h = h,dta = z,K = K)
  # Estimated pdf
  plot(z2,f_x,type = "l",ylim = c(0,max(c(f_x,fn_x))))
  # Actual pdf
  lines(z2,fn_x,lty=2,col="red")
  diff = abs(fn_x-f_x)
  ind = which(diff == max(diff))
  segments(x0 = z2[ind],y0 = f_x[ind],x1 = z2[ind],y1 = fn_x[ind])
}

Vn_draw(n = ,K = K_opt,h = h_opt,alpha = 0.1)



# Large Sample Distribution
n = 10^3
sim = 10^4
sg = 1
h_opt = 1.06*sg*(n)^(-1/5)


f1=NULL
for(i in 1:sim)
{
  data2=rnorm(n)
  k=(f_kernel(t = t,K = K_expo,h = h_opt,dta = data2)-m)/sqrt(v)
  f1=c(f1,k)
}
hist(f1,prob = TRUE)
sq = seq(-3,3,length.out = 100)
lines(sq,dnorm(sq))

# When the kernel taken, doesn't satisfy the required conditions :-
n = 10^3
sim = 10^4
sg = 1
h_opt = 1.06*sg*(n)^(-1/5)


f1=NULL
for(i in 1:sim)
{
  data2=rnorm(n)
  k=(f_kernel(t = t,K = K_expo,h = h_opt,dta = data2)-m)/sqrt(v)
  f1=c(f1,k)
}
hist(f1,prob = TRUE)
sq = seq(-3,3,length.out = 100)
lines(sq,dnorm(sq))
# Some interesting observations to make here :D

# How the estimated density varies near the actual density
# For beta density
n = 10^3
x = rbeta(n = n,shape1 = 2,shape2 = 4)
r = range(x)
sq = seq(0,1,length.out = 10^2)
fx = dbeta(sq,shape1 = 2,shape2 = 4)
p = max(fx)
plot(sq,fx,type = "l",col = "black",lwd = 1.5,ylim = c(0,p))
sim = 10

for(i in 1:sim)
{
  x = rbeta(n = n,shape1 = 2,shape2 = 4)
  r = range(x)
  den_est = f_kernel_vec(t_vec = sq,h = 0.1,dta = x,K = K_opt)
  lines(sq,den_est,col = "red",lwd = 0.5,lty = 2)
}

