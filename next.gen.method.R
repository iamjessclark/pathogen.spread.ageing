# Next generation method R0 analytically solved in Mathematica 
# Reviewer comments Am Nat 
# Started Febraury 2020

#### Packages ####

require(RColorBrewer)
require(fields)
require(ggfortify)
require(ggplot2)
require(lme4)
require(car)
require(MASS)
require(fitdistrplus)
require(reshape2)
require(survival)
require(plyr)
require(lmerTest)
require(survminer)
require(aod)
require(effects)
require(xtable)
require(cowplot)
library("gridExtra")
require(ggpubr)
require(tidyverse)
require(flexsurv)
require(eha)
require(interplot)
require(sjPlot)
require(nlme)

#### functions ####


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


#### Parameters ####

# all rates are per day
r <- 5 # reproduction 
k <- r # max pop density
m <- 0.5 # maturation rate (average age at maturation is 10 days)
v <- 1.5 # virulence 
c <- 0.1 
B <- 0.9 # baseline value

low_d <- 0.1 # low mortality
med_d <- 0.5 # medium mortality
high_d <- 1 # high mortality 
from=-1 # plots
to = 1 # plots 

effect_length = 501 # length out for sequences 

#### === === === === === === === === === === === === === ===  #### 
#### === === === ===     Main Effects    === === === === ===  #### 
#### === === === === === === === === === === === === === ===  #### 

#### Dr0Dx ####

x <- seq(from = from, to = to, length.out = effect_length)
y <- 0
M <- 0 # maternal effect on susceptibility
A <- 0 # ageing effect on susceptibility

# Low mortality 
d = low_d

R0xlow = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# Medium mortality 
d = med_d
R0xmed = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# High Mortality 
d = high_d
R0xhigh = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### Dr0Dy ####

x <- 0
y <- seq(from = from, to = to, length.out = effect_length)
M <- 0 # maternal effect on susceptibility
A <- 0 # ageing effect on susceptibility

# Low mortality 
d = low_d

R0ylow = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# Medium mortality 
d = med_d
R0ymed = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# High Mortality 
d = high_d
R0yhigh = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)


#### Dr0DM ####

x <- 0
y <- 0
M <- seq(from = from, to = to, length.out = effect_length) # maternal effect on susceptibility
A <- 0 # ageing effect on susceptibility

# Low mortality 
d = low_d

R0Mlow = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# Medium mortality 
d = med_d
R0Mmed = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# High Mortality 
d = high_d
R0Mhigh = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### Dr0DA ####

x <- 0
y <- 0
M <- 0
A <- seq(from = from, to = to, length.out = effect_length) 

# Low mortality 
d = low_d

R0Alow = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# Medium mortality 
d = med_d
R0Amed = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

# High Mortality 
d = high_d
R0Ahigh = (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### Density ####

#### dnudx ####

y <- 0
x <- seq(from = -1, to = 1, length.out = effect_length)
dnudx <- matrix(nrow = length(x), ncol = 4)
dnudx2 <- matrix(nrow = length(x), ncol = 4)
dnudx3 <- matrix(nrow = length(x), ncol = 4)

d <- low_d
dnudx[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudx[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudx[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudx[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudx[dnudx<0] <- 0

d <- med_d
dnudx2[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudx2[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudx2[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudx2[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudx2[dnudx2<0] <- 0

d <- high_d
dnudx3[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudx3[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudx3[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudx3[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudx3[dnudx3<0] <- 0

#### dnudy ####

x <- 0
y <- seq(from = -1, to = 1, length.out = effect_length)

dnudy <- matrix(nrow = length(y), ncol = 4)
dnudy2 <- matrix(nrow = length(y), ncol = 4)
dnudy3 <- matrix(nrow = length(y), ncol = 4)

d <- low_d
dnudy[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudy[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudy[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudy[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudy[dnudy<0] <- 0


d <- med_d
dnudy2[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudy2[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudy2[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudy2[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudy2[dnudy2<0] <- 0

d <- high_d
dnudy3[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudy3[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudy3[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudy3[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudy3[dnudy3<0] <- 0

#### Main effects plot ####

axis.text = 1.5

x <- seq(from = from, to = to, length.out = effect_length)
y <- seq(from = from, to = to, length.out = effect_length)
A <- seq(from = from, to = to, length.out = effect_length)
M <- seq(from = from, to = to, length.out = effect_length)

dev.new(height = 15, width = 6)
par(mar=c(4,5,3,3))
par(mfrow=c(5,2), mgp=c(3,1,0))

plot(NA,NA, xlim = c(-1,1), ylim = c(0,2.75), xlab="Fecundity Senescence", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
lines(R0xlow~x, col="#8B0000", lwd=2) # low dark red
lines(R0xmed~x, col="#EE0000", lwd=2) # low dark red
lines(R0xhigh~x, col="#E9967A", lwd=2) # low dark red

plot(NA,NA, xlim = c(-1,1), ylim = c(0,2.75), xlab="Mortality Senescence", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
lines(R0ylow~y, col="#8B0000", lwd=2) # low dark #EE0000
lines(R0ymed~y, col="#EE0000", lwd=2) # low dark red
lines(R0yhigh~y, col="#E9967A", lwd=2) # low dark red

plot(NA,NA, xlim = c(-1,1), ylim = c(0,5), xlab=" Transgenerational Defense", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
lines(R0Mlow~M, col="#8B0000", lwd=2) # low dark #EE0000
lines(R0Mmed~M, col="#EE0000", lwd=2) # low dark red
lines(R0Mhigh~M, col="#E9967A", lwd=2) # low dark red

plot(NA,NA, xlim = c(-1,1), ylim = c(0,5), xlab="Age Specific Defense", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
lines(R0Alow~A, col="#8B0000", lwd=2) # low dark red
lines(R0Amed~A, col="#EE0000", lwd=2) # low dark red
lines(R0Ahigh~A, col="#E9967A", lwd=2) # low dark red

plot(NA,NA,xlim=c(-1,1),ylim=c(0,5), xlab = " Reproduction ", ylab="Total density", cex.lab=axis.text, cex.axis=axis.text)
lines(dnudx[,1]~x,lwd=2, lty = 1,col="firebrick1")

lines(dnudx[,2]~x,lwd=2, lty = 1,col="dodgerblue1")#you

lines(dnudx[,3]~x,lwd=2, lty=1,col="goldenrod1")

lines(dnudx[,4]~x,lwd=2, lty=1,col="forestgreen")

lines(dnudx[,1]+dnudx[,2]+dnudx[,3]+dnudx[,4]~x, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,5), xlab = "Mortality ", ylab="Total density",cex.lab=axis.text, cex.axis=axis.text)
lines(dnudy[,1]~y,lwd=2, lty = 1,col="firebrick1")

lines((dnudy[,2]-0.02)~y,lwd=3, lty = 1,col="dodgerblue1")#you

lines(dnudy[,3]~y,lwd=2, lty=1,col="goldenrod1")

lines(dnudy[,4]~y,lwd=2, lty=1,col="forestgreen")

lines(dnudy[,1]+dnudy[,2]+dnudy[,3]+dnudy[,4]~y, lwd=2, lty=1,col="blueviolet")


plot(NA,NA,xlim=c(-1,1),ylim=c(0,5), xlab = " Reproduction ", ylab="Total density", cex.lab=axis.text, cex.axis=axis.text)
lines((dnudx2[,1]-0.02)~x,lwd=3, lty = 1,col="firebrick1")

lines((dnudx2[,2]-0.02)~x,lwd=2, lty = 1,col="dodgerblue1")#you

lines(dnudx2[,3]~x,lwd=2, lty=1,col="goldenrod1")

lines(dnudx2[,4]~x,lwd=2, lty=1,col="forestgreen")

lines(dnudx2[,1]+dnudx2[,2]+dnudx2[,3]+dnudx2[,4]~x, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,5), xlab = " Mortality ", ylab="Total density",cex.lab=axis.text, cex.axis=axis.text)
lines(dnudy2[,1]~y,lwd=2, lty = 1,col="firebrick1")

lines((dnudy2[,2]-0.02)~y,lwd=3, lty = 1,col="dodgerblue1")#you

lines(dnudy2[,3]~y,lwd=2, lty=1,col="goldenrod1")

lines(dnudy2[,4]~y,lwd=2, lty=1,col="forestgreen")

lines(dnudy2[,1]+dnudy2[,2]+dnudy2[,3]+dnudy2[,4]~y, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,5), xlab = " Reproduction ", ylab="Total density", cex.lab=axis.text, cex.axis=axis.text)
lines(dnudx3[,1]~x,lwd=2, lty = 1,col="firebrick1")

lines(dnudx3[,2]~x,lwd=2, lty = 1,col="dodgerblue1")#you

lines(dnudx3[,3]~x,lwd=2, lty=1,col="goldenrod1")

lines(dnudx3[,4]~x,lwd=2, lty=1,col="forestgreen")

lines(dnudx3[,1]+dnudx3[,2]+dnudx3[,3]+dnudx3[,4]~x, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,5), xlab = " Mortality ", ylab="Total density",cex.lab=axis.text, cex.axis=axis.text)
lines(dnudy3[,1]~y,lwd=2, lty = 1,col="firebrick1")

lines((dnudy3[,2]-0.02)~y,lwd=3, lty = 1,col="dodgerblue1")#you

lines(dnudy3[,3]~y,lwd=2, lty=1,col="goldenrod1")

lines(dnudy3[,4]~y,lwd=2, lty=1,col="forestgreen")

lines(dnudy3[,1]+dnudy3[,2]+dnudy3[,3]+dnudy3[,4]~y, lwd=2, lty=1,col="blueviolet")


#add_legend("bottom", legend=c(c("Uyy","Uyo","Uoy","Uoo","Total")), lty=1, lwd=2,
           #col=c("firebrick1","dodgerblue1","goldenrod1","forestgreen", "blueviolet"),
           #horiz=TRUE, bty='n', cex=1.2)

#add_legend("top", legend=c(expression(paste(delta,"= 0.1")), expression(paste(delta,"= 0.5")), expression(paste(delta,"= 1.5"))), lty=1, lwd=2,
           #col=c("#8B0000", "#EE0000", "#E9967A"),
           #horiz=TRUE, bty='n', cex=1.2)

dev.copy(pdf, "main.effects.NGM.26.2.20.pdf", height = 17, width = 7)
dev.off()
graphics.off()

#### === === === === === === === === === === === === === ===  #### 
#### === === === ===     Interactions    === === === === ===  #### 
#### === === === === === === === === === === === === === ===  #### 


#### A by M ####

A = matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length)
M = t(matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length))
x=0
y=0

d = low_d
r01 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = med_d
r02 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = high_d
r03 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### x by y ####

x = matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length)
y = t(matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length))
A=0
M=0

d = low_d
r04 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = med_d
r05 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = high_d
r06 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

r04[which(is.nan(r04))] <- NA
r05[which(is.nan(r05))] <- NA
r05[which(is.nan(r05))] <- NA

#### x by A ####

x = matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length)
A = t(matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length))
y=0
M=0

d = low_d
r07 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = med_d
r08 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = high_d
r09 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### x by M ####

x = matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length)
M = t(matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length))
y=0
A=0

d = low_d
r010 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = med_d
r011 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = high_d
r012 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### y by A ####

y = matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length)
A = t(matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length))
x=0
M=0

d = low_d
r013 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = med_d
r014 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = high_d
r015 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### y by M ####

y = matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length)
M = t(matrix(rep(seq(from = -1, to = 1, length.out = effect_length), each=effect_length),nrow = effect_length))
x=0
A=0

d = low_d
r016 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = med_d
r017 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)
d = high_d
r018 <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

#### Interaction Plot ####

#### Z axis limits #### 

min(c(r01,r02,r03,r04,r05,r06,r07,r08,r09,r010,r011,r012,r013,r014,r015,r016,r017,r018),na.rm=T)
max(c(r01,r02,r03,r04,r05,r06,r07,r08,r09,r010,r011,r012,r013,r014,r015,r016,r017,r018),na.rm=T)
zmin <- 0
zmax <- ceiling(max(c(r01,r02,r03,r04,r05,r06,r07,r08,r09,r010,r011,r012,r013,r014,r015,r016,r017,r018),na.rm=T))

cols<-(colorRampPalette(brewer.pal(9,"YlOrRd"),bias=1.5,interpolate="linear")(501))

nlevels<-10 # number of countours

dev.new(height = 10, width = 8)

par(mfrow = c(6,3))
par(mar=c(4,6,2,5),mgp=c(3,1,0))
par(oma=c(1,1,1,1))

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r01)), xlab= expression("Trans. Gen. Resis"~italic("(M)")), ylab= expression("Resistance"~italic("(A)")),col=cols,zlim=c(min(r01),max(r01)),useRaster=T, cex.lab = 1.5)
#mtext( 'A', side=90, line=656, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r01)),zlim=c(min((r01),na.rm=T),max((r01),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r01),na.rm=T),max((r01),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r02)), xlab=expression("Trans. Gen. Resis"~italic("(M)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min(r02),max(r02)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r02)),zlim=c(min((r02),na.rm=T),max((r02),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(0,1), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r03)), xlab=expression("Trans. Gen. Resis"~italic("(M)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min(r03),max(r03)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r03)),zlim=c(min((r03),na.rm=T),max((r03),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r03),na.rm=T),max((r03),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r04)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Mortality"~italic("(y)")),col=cols,zlim=c(min((r04),na.rm=T),max((r04),na.rm=T)),useRaster=T, cex.lab = 1.5)
#mtext( 'B', side=90, line=541, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r04)),zlim=c(min((r04),na.rm=T),max((r04),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r04),na.rm=T),max((r04),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r05)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Mortality"~italic("(y)")),col=cols,zlim=c(min((r05),na.rm=T),max((r05),na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r05)),zlim=c(min((r05),na.rm=T),max((r05),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r05),na.rm=T),max((r05),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r06)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Mortality"~italic("(y)")),col=cols,zlim=c(min((r06),na.rm=T),max((r06),na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r06)),zlim=c(min((r06),na.rm=T),max((r06),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r06),na.rm=T),max((r06),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r07)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min((r07),na.rm=T),max((r07),na.rm=T)),useRaster=T, cex.lab = 1.5)
#mtext( 'C', side=90, line=422, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r07)),zlim=c(min((r07),na.rm=T),max((r07),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r07),na.rm=T),max((r07),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r08)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min(r08),max(r08)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r08)),zlim=c(min((r08),na.rm=T),max((r08),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(0,10), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r09)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min((r09), na.rm=T),max((r09), na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r09)),zlim=c(min((r09),na.rm=T),max((r09),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r09),na.rm=T),max((r09),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r010)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Trans. Gen. Resis"~italic("(M)")),col=cols,zlim=c(min(r010),max(r010)),useRaster=T, cex.lab = 1.5)
#mtext( 'D', side=90, line=310, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r010)),zlim=c(min((r010),na.rm=T),max((r010),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r010),na.rm=T),max((r010),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r011)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Trans. Gen. Resis"~italic("(M)")),col=cols,zlim=c(min(r011),max(r011)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r011)),zlim=c(min((r011),na.rm=T),max((r011),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r011),na.rm=T),max((r011),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r012)), xlab=expression("Reproduction"~italic("(x)")), ylab=expression("Trans. Gen. Resis"~italic("(M)")),col=cols,zlim=c(min((r012), na.rm=T),max((r012), na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r012)),zlim=c(min((r012),na.rm=T),max((r012),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r012),na.rm=T),max((r012),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r013)), xlab=expression("Mortality"~italic("(y)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min(r013),max(r013)),useRaster=T, cex.lab = 1.5)
#mtext( 'E', side=90, line=190, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r013)),zlim=c(min((r013),na.rm=T),max((r013),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r013),na.rm=T),max((r013),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r014)), xlab=expression("Mortality"~italic("(y)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min(r014),max(r014)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r014)),zlim=c(min((r014),na.rm=T),max((r014),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r014),na.rm=T),max((r014),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r015)), xlab=expression("Mortality"~italic("(y)")), ylab=expression("Resistance"~italic("(A)")),col=cols,zlim=c(min(r015),max(r015)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r015)),zlim=c(min((r015),na.rm=T),max((r015),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r015),na.rm=T),max((r015),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r016)), xlab=expression("Mortality"~italic("(y)")), ylab=expression("Trans. Gen. Resis"~italic("(M)")),col=cols,zlim=c(min(r016),max(r016)),useRaster=T, cex.lab = 1.5)
#mtext( 'F', side=90, line=78, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r016)),zlim=c(min((r016),na.rm=T),max((r016),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r016),na.rm=T),max((r016),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r017)), xlab=expression("Mortality"~italic("(y)")), ylab=expression("Trans. Gen. Resis"~italic("(M)")),col=cols,zlim=c(min(r017),max(r017)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r017)),zlim=c(min((r017),na.rm=T),max((r017),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r017),na.rm=T),max((r017),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r018)), xlab=expression("Mortality"~italic("(y)")), ylab=expression("Trans. Gen. Resis"~italic("(M)")),col=cols,zlim=c(min(r018),max(r018)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = effect_length), seq(-1, 1, length.out = effect_length), z=t((r018)),zlim=c(min((r018),na.rm=T),max((r018),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r018),na.rm=T),max((r018),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

dev.copy(pdf, "interactions.figure.5.3.20.pdf", height = 18, width = 12)
dev.off()
graphics.off()

#### === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === === ===  ####

#### slowing senescence ####

d <- 0.1 # Always using low mortality 

from = 0.5
to = -1

x <- seq(from = from, to = to, length.out = effect_length)
y <- 0.5
M <- 0.5
A <- 0.5

dr0x.slow.sen <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

y <- seq(from = from, to = to, length.out = effect_length)
x <- 0.5
M<- 0.5
A<- 0.5

dr0y.slow.sen <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)


M <- seq(from = from, to = to, length.out = effect_length)
x <- 0.5
y<- 0.5
A<- 0.5

dr0M.slow.sen <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)


A <- seq(from = from, to = to, length.out = effect_length)
x <- 0.5
M<- 0.5
y<- 0.5

dr0A.slow.sen <- (B*k*(((-1 + A)*d - (1 + A)*m)*(c + d + m + v) + d*((-1 + A)*c + 2*(-1 + A)*d + 2*A*m + (-1 + A)*v)*y + (-1 + A)*d^2*y^2)*(m*(1 + M)*(-1 + x) + d*(-1 + M)*(1 + x)*(1 + y))*(-m*r*(-1 + x) - d*(m - r*(1 + x))*(1 + y) + d^2*(-1 + y^2)))/(r*(c + d + m + v - d*y)*(d + m + d*y)*(c + d + v + d*y)*(m - m*x +d*(1 + x)*(1 + y))^2)

from = 0
to = 1.5

x <- seq(from = 0, to = 1.5, length.out = effect_length)

dev.new(height = 4, width = 4)
par(mar=c(5,4,1,1))
par(mfrow=c(1,1), mgp=c(3,1,0))

plot(NA, NA, xlim = c(from,to), ylim = c(0,7), xlab = "Senescence", ylab = "Transmission Potential",cex.lab=1, cex.axis=1, bty="n")
#title(outer=TRUE,adj=0.23,main="A",cex=1, col="black",font=1,line=-1.85)
lines(dr0x.slow.sen ~ x,lwd=2, lty = 1,col="#EBCC2A")
lines(dr0y.slow.sen ~ x,lwd=2, lty = 1,col="#3B9AB2")
lines(dr0A.slow.sen ~ x,lwd=2, lty = 1,col="#b2df8a")
lines(dr0M.slow.sen ~ x,lwd=2, lty = 1,col="#e78ac3")


add_legend("top", legend=c("Reprodution", "Mortality", "Susceptibility", "Trans-gen Susceptibility"), lty=1, lwd=3,
           col=c("#EBCC2A", "#3B9AB2", "#b2df8a", "#e78ac3"),
           horiz=TRUE, bty='n', cex=0.7)

dev.copy(pdf, "slow.sen.5.3.20.pdf", height = 5.5, width = 6)
dev.off()
graphics.off()

#### SM pop age structure and density ####

y <- 0.5
x <- seq(from =0.5, to = 1.5, length.out = effect_length)
dnudx <- matrix(nrow = length(x), ncol = 4)

d <- low_d
dnudx[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudx[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudx[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudx[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudx[dnudx<0] <- 0

x <- 0.5
y <- seq(from = 0.5, to = 1.5, length.out = effect_length)

dnudy <- matrix(nrow = length(y), ncol = 4)

d <- low_d
dnudy[,1] <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dnudy[,2] <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dnudy[,3] <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dnudy[,4] <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dnudy[dnudy<0] <- 0


axis.text <- 1
x <- seq(from =0.5, to = 1.5, length.out = effect_length)
y <- seq(from = 0.5, to = 1.5, length.out = effect_length)

dev.new(height = 11, width = 6)
par(mar=c(5,4,1,1))
par(mfrow=c(2,1), mgp=c(3,1,0))


plot(NA,NA,xlim=c(0.5,1.5),ylim=c(0,5), xlab = " Reproduction ", ylab="Total density", cex.lab=axis.text, cex.axis=axis.text)
lines(dnudx[,1]~x,lwd=2, lty = 1,col="firebrick1")

lines(dnudx[,2]~x,lwd=2, lty = 1,col="dodgerblue1")#you

lines(dnudx[,3]~x,lwd=2, lty=1,col="goldenrod1")

lines(dnudx[,4]~x,lwd=2, lty=1,col="forestgreen")

lines(dnudx[,1]+dnudx[,2]+dnudx[,3]+dnudx[,4]~x, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(0.5,1.5),ylim=c(0,5), xlab = "Mortality ", ylab="Total density",cex.lab=axis.text, cex.axis=axis.text)
lines(dnudy[,1]~y,lwd=2, lty = 1,col="firebrick1")

lines((dnudy[,2]-0.02)~y,lwd=3, lty = 1,col="dodgerblue1")#you

lines(dnudy[,3]~y,lwd=2, lty=1,col="goldenrod1")

lines(dnudy[,4]~y,lwd=2, lty=1,col="forestgreen")

lines(dnudy[,1]+dnudy[,2]+dnudy[,3]+dnudy[,4]~y, lwd=2, lty=1,col="blueviolet")

dev.copy(pdf, "sm.density.5.3.20.pdf", height = 11, width = 6)
dev.off()
graphics.off()
