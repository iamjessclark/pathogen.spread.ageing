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
require(binom)

# function for adding legend 

add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#### Experimental Data ####

#use mat.effect.data.SA for unexposed individuals

#use mat.effect.data for all individuals

#use infection for only exposed individuals

#early.repro has only the reprouctive rate for early reproduction. This was a new df made with limited data and melted to be long. 

#late.repro.df is the same as above but for late rate repro. 

#weekly.repro is the melted data set with only unexposed indivudals with their total reproduction per week as the value

#daily.repro df is for the breakdown of babies per day where day is maternal age.

#ici.repro df is used for interclutch interval analysis and number of offspring/ clutch analysis (ici.repro2 was the df created when melting mat.effect.data.SA it is not used for anything else)

#### Mortality Data old mums ####

#### Data ####

mortality <- read.csv("Controls_infected.csv")

mortality.count <- mortality[-which(is.na(mortality$days.alive)==TRUE),]

mortality.count <- subset(mortality.count, F0.food == "H")

mortality.count %>%
  group_by(F1.food, E) %>%
  summarise(nrows= length(F1.food))

mortality$treatment <- as.factor(paste(mortality$F1.food, mortality$E, sep =""))

mortality.unex <- subset(mortality, E == "U")

str(mortality$treatment)

mortality$E <-relevel(mortality$E, ref='U') # to change the reference level for the model 

fitweib <- flexsurvreg(formula = Surv(days.alive) ~ F1.food + E, data = mortality, dist="weibull")
fitweib
fitweib$res.t # these parameters are on the log-scale, use fitweib$res if you want them not on the log scale

# The results from the weibreg model are presented in the same parameterization as flexsurvreg, except 
# that the shape and scale parameters are log-transformed, if you want the out put of the flexsurvreg model to match use the res.t output as it is
# transformed to the log scale
# and (unless the argument param = "lifeExp" is supplied) the covariate effects have the opposite sign.

# wiebull hazards

hazards.weib <- summary(fitweib, type = "hazard", fn = NULL, t = NULL, start = 0, ci = TRUE, B = 1000, cl = 0.95, tidy = FALSE)
surv.weib <- summary(fitweib, type = "survival", fn = NULL, t = NULL, start = 0, ci = TRUE, B = 1000, cl = 0.95, tidy = FALSE)

hazards.HE.weib <- hazards.weib[["F1.food=H,E=E"]]
HE <- "HE"
hazards.HE.weib$treatment <- as.factor(HE)

surv.HE.weib <- surv.weib[["F1.food=H,E=E"]]
HE <- "HE"
surv.HE.weib$treatment <- as.factor(HE)

hazards.HU.weib <- hazards.weib[["F1.food=H,E=U"]]
HU <- "HU"
hazards.HU.weib$treatment <- as.factor(HU)

surv.HU.weib <- surv.weib[["F1.food=H,E=U"]]
HU <- "HU"
surv.HU.weib$treatment <- as.factor(HU)

hazards.LU.weib <- hazards.weib[["F1.food=L,E=U"]]
LU <- "LU"
hazards.LU.weib$treatment <- as.factor(LU)

surv.LU.weib <- surv.weib[["F1.food=L,E=U"]]
LU <- "LU"
surv.LU.weib$treatment <- as.factor(LU)

hazards.LE.weib <- hazards.weib[["F1.food=L,E=E"]]
LE <- "LE"
hazards.LE.weib$treatment <- as.factor(LE)

surv.LE.weib <- surv.weib[["F1.food=L,E=E"]]
LE <- "LE"
surv.LE.weib$treatment <- as.factor(LE)


# plot hazards

dev.new(height = 6, width = 6)
par(mfrow=c(2,1))
par(mar=c(6,6,4,1),mgp=c(3,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(14,125), ylim = c(0,2.5), xlab="Days Alive", ylab="Hazard", cex.axis=1, cex.lab=1, bty="l")
lines(hazards.HE.weib$est~hazards.HE.weib$time, col="#8B0000", lwd=2, lty = 2) # high dark red
lines(hazards.HU.weib$est~hazards.HE.weib$time, col="#8B0000", lwd=2) # high dark red
lines(hazards.LE.weib$est~hazards.LE.weib$time, col="#E9967A", lwd=2, lty = 2) # low salmon
lines(hazards.LU.weib$est~hazards.LU.weib$time, col="#E9967A", lwd=2) # low salmon
par(fig = c(0.1,0.6, 0.5, 0.95), new = T, mgp = c(2,1,0))  
plot(NA,NA, xlim = c(14,125), ylim = c(0,1), xlab="Days Alive", ylab="Survival", cex.axis=1, cex.lab=1, bty="l")
lines(surv.HE.weib$est~surv.HE.weib$time, col="#8B0000", lwd=2, lty = 2) # high dark red
lines(surv.HU.weib$est~surv.HE.weib$time, col="#8B0000", lwd=2) # high dark red
lines(surv.LE.weib$est~surv.LE.weib$time, col="#E9967A", lwd=2, lty = 2) # low salmon
lines(surv.LU.weib$est~surv.LU.weib$time, col="#E9967A", lwd=2) # low salmon

dev.copy(pdf, "weib.mortality.pdf", height = 6, width = 6)
dev.off()
graphics.off()

#-----------------------------------------------------------

####Load Offspring Data####

mat.effect.data<-read.csv("mat.effect.sen.data.csv")
mat.effect.data <- mat.effect.data[mat.effect.data$infected !="DBD",]

mat.effect.data%>%
  group_by(Maternal.Food, clutch, exposed) %>%
  summarise(no_rows = length(Maternal.Food))
#-----------------------------------------------------------

#### Infection ####

mat.effect.data.infection<-subset(mat.effect.data, exposed == "E")
mat.effect.data.infection<-mat.effect.data.infection[mat.effect.data.infection$infected !="DBD",]
mat.effect.data.infection$Treatment<-paste(mat.effect.data.infection$Maternal.Food, mat.effect.data.infection$clutch, sep="")

# make a table to figure out proportion infected, then plot the table

infection.count<-xtabs(~infected + Treatment, mat.effect.data.infection)
infection.count

# add together the columns into a total and add the total into the df of counts

Total<-colSums(infection.count)
Total
infection.count<-rbind(infection.count, Total)
infection.count

Prop.Infect<-infection.count[2,]/infection.count[5,]#dividing row 2 by row 5 (infected/total sample)
Prop.Infect
infection.count<-rbind(infection.count, Prop.Infect)
infection.count

infection.count<-melt(infection.count, id.var = c ("LFA","LFB","LFC","NFA","NFB","NFC"), measure.vars=c("Total", "Prop.Infect"))
infection.count
infection.count<-infection.count[infection.count$Var1 !="DBD",]
infection.count
infection.count<-infection.count[infection.count$Var1 !="N",]
infection.count
infection.count<-infection.count[infection.count$Var1 !="NE",]
infection.count

# make data wide again

infection.count<-reshape(infection.count, idvar = "Var2", timevar = "Var1", direction = "wide")
infection.count
Maternal.Food <- "Maternal.Food"
infection.count[ , "Maternal.Food"] <- NA
infection.count[1:3,5]<-"LF"
infection.count[4:6,5]<-"NF"

Clutch <- "Clutch"
infection.count[ , "Clutch"] <- NA
infection.count[1,6]<-"1"
infection.count[2,6]<-"2"
infection.count[3,6]<-"Old"
infection.count[4,6]<-"1"
infection.count[5,6]<-"2"
infection.count[6,6]<-"Old"

infection.count$mat.food<-as.factor(infection.count$Maternal.Food)
infection.count$Clutch<-as.factor(infection.count$Clutch)
colnames(infection.count)[which(names(infection.count)== "Var2")] <- "Treatment"

# 

lfa.conf.int <- binom.confint(x=29, n=48,  methods = "wilson")
lfb.conf.int <- binom.confint(x=10, n=24,  methods = "wilson")
lfc.conf.int <- binom.confint(x=10, n=37,  methods = "wilson")

nfa.conf.int <- binom.confint(x=35, n=50,  methods = "wilson")
nfb.conf.int <- binom.confint(x=27, n=56,  methods = "wilson")
nfc.conf.int <- binom.confint(x=12, n=50,  methods = "wilson")

conf.int.inf <- data.frame(infection.count$Maternal.Food, infection.count$Clutch)
colnames(conf.int.inf) <- c("Maternal.Food", "Clutch")
conf.int.inf$upper <- NA
conf.int.inf$lower <- NA

conf.int.inf[1,3] <- 0.72
conf.int.inf[1,4] <- 0.46
conf.int.inf[2,3] <- 0.61
conf.int.inf[2,4] <- 0.24
conf.int.inf[3,3] <- 0.43
conf.int.inf[3,4] <- 0.15
conf.int.inf[4,3] <- 0.82
conf.int.inf[4,4] <- 0.55
conf.int.inf[5,3] <- 0.62
conf.int.inf[5,4] <- 0.35
conf.int.inf[6,3] <- 0.38
conf.int.inf[6,4] <- 0.13

# with error bars as suggested by the reviewer but it doesn't make sense. 
infection.bar <- ggplot(infection.count, aes(x = Clutch, y = value.Prop.Infect, fill=Maternal.Food))+ 
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=conf.int.inf$lower, ymax=conf.int.inf$upper),width=.2,position=position_dodge(.9))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=16))+
  scale_fill_manual(values=c("#FFA07A","#DC143C"),
                    labels = c("DR", "Normal"))+ 
  ylab("Proportion Infected")+xlab("Clutch")+
  geom_text(aes(label=value.I), position=position_dodge(width=0.9), vjust=-0.25)+
  ggsave("confint_susceptibility.pdf")

infection.bar 

# without error bars

infection.bar <- ggplot(infection.count, aes(x = Clutch, y = value.Prop.Infect, fill=Maternal.Food))+ 
  geom_bar(stat="identity", position=position_dodge())+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  scale_fill_manual(values=c("#FFA07A","#DC143C"),
                    labels = c("DR", "Normal"))+ 
  ylab("Proportion Infected")+xlab("Clutch")+
  geom_text(aes(label=value.I), position=position_dodge(width=0.9), vjust=-0.25)

infection.bar 
# This will save a 400x400 file at 100 ppi
ggsave("infection.bar.png", width=8, height=6, dpi=300)

####Analysis####

mat.effect.data.infection<-subset(mat.effect.data, exposed == "E")
mat.effect.data.infection<-mat.effect.data.infection[mat.effect.data.infection$infected !="DBD",]

#logistic regression
infection<-glm(infected~mat.food + clutch + clutch*mat.food, family=binomial(link='logit'), data=mat.effect.data.infection)
summary(infection)
#interaction not significant

#without interaction
infection2<-glm(infected~mat.food + clutch, family=binomial(link='logit'), data=mat.effect.data.infection)
summary(infection2)


#anova to assess deviance 
infection.dev<-anova(infection2, test = "Chisq")
infection.dev

wald.test(b = coef(infection2), Sigma =vcov(infection2), Terms = 3:4)#obtains a chi squared value for the specific term in the model - clutch
#Wald test:
#Chi-squared test:
#X2 = 28.2, df = 2, P(> X2) = 7.7e-07


wald.test(b = coef(infection2), Sigma =vcov(infection2), Terms = 2)#just two because not the intercept - food
#X2 = 0.57, df = 1, P(> X2) = 0.45

#### Subset unexposed ####

mat.effect.data.SA<-subset(mat.effect.data, exposed == "U")
mat.effect.data.SA$Treatment<-as.factor(paste(mat.effect.data.SA$Maternal.Food, mat.effect.data.SA$clutch, sep=""))
str(mat.effect.data.SA$Treatment)

#### Suvival Analysis ####

# Gompertz fits best

fitg.f1 <- flexsurvreg(formula = Surv(days.alive) ~ Maternal.Food + clutch, data = mat.effect.data.SA, dist="Gompertz")
fitg.f1
fitg.f1$res.t

fitexp.f1 <- flexsurvreg(formula = Surv(days.alive) ~ Maternal.Food + clutch, data = mat.effect.data.SA, dist="Exponential")
fitexp.f1

fitweib.f1 <- flexsurvreg(formula = Surv(days.alive) ~ Maternal.Food + clutch, data = mat.effect.data.SA, dist="weibull")
fitweib.f1
fitweib.f1$res.t # these parameters are on the log-scale, use fitweib$res if you want them not on the log scale

# gompertz hazards

hazards.GOMP.f1 <- summary(fitg.f1, type = "hazard", fn = NULL, t = NULL, start = 0, ci = TRUE, B = 1000, cl = 0.95, tidy = FALSE)
surv.GOMP.f1 <- summary(fitg.f1, type = "survival", fn = NULL, t = NULL, start = 0, ci = TRUE, B = 1000, cl = 0.95, tidy = FALSE)

hazards.NFA.GOMP.f1 <- hazards.GOMP.f1[["Maternal.Food=NF,clutch=A"]]
NFA <- "NFA"
hazards.NFA.GOMP.f1$treatment <- as.factor(NFA)

surv.NFA.GOMP.f1 <- surv.GOMP.f1[["Maternal.Food=NF,clutch=A"]]
NFA <- "NFA"
surv.NFA.GOMP.f1$treatment <- as.factor(NFA)

hazards.NFB.GOMP.f1 <- hazards.GOMP.f1[["Maternal.Food=NF,clutch=B"]]
NFB <- "NFB"
hazards.NFB.GOMP.f1$treatment <- as.factor(NFB)

surv.NFB.GOMP.f1 <- surv.GOMP.f1[["Maternal.Food=NF,clutch=B"]]
NFB <- "NFB"
surv.NFB.GOMP.f1$treatment <- as.factor(NFB)

hazards.NFC.GOMP.f1 <- hazards.GOMP.f1[["Maternal.Food=NF,clutch=C"]]
NFC <- "NFC"
hazards.NFC.GOMP.f1$treatment <- as.factor(NFC)

surv.NFC.GOMP.f1 <- surv.GOMP.f1[["Maternal.Food=NF,clutch=C"]]
NFC <- "NFC"
surv.NFC.GOMP.f1$treatment <- as.factor(NFC)

hazards.LFA.GOMP.f1 <- hazards.GOMP.f1[["Maternal.Food=LF,clutch=A"]]
LFA <- "LFA"
hazards.LFA.GOMP.f1$treatment <- as.factor(LFA)

surv.LFA.GOMP.f1 <- surv.GOMP.f1[["Maternal.Food=LF,clutch=A"]]
LFA <- "LFA"
surv.LFA.GOMP.f1$treatment <- as.factor(LFA)

hazards.LFB.GOMP.f1 <- hazards.GOMP.f1[["Maternal.Food=LF,clutch=B"]]
LFB <- "LFB"
hazards.LFB.GOMP.f1$treatment <- as.factor(LFB)

surv.LFB.GOMP.f1 <- surv.GOMP.f1[["Maternal.Food=LF,clutch=B"]]
LFB <- "LFB"
surv.LFB.GOMP.f1$treatment <- as.factor(LFB)

hazards.LFC.GOMP.f1 <- hazards.GOMP.f1[["Maternal.Food=LF,clutch=C"]]
LFC <- "LFC"
hazards.LFC.GOMP.f1$treatment <- as.factor(LFC)

surv.LFC.GOMP.f1 <- surv.GOMP.f1[["Maternal.Food=LF,clutch=C"]]
LFC <- "LFC"
surv.LFC.GOMP.f1$treatment <- as.factor(LFC)

# plot hazards

axis.text <- 1

dev.new(height = 6, width = 6)
par(mfrow=c(1,1))
par(mar=c(3,4,2,1),mgp=c(2,1,0))
par(fig = c(0,1,0,1))

plot(NA,NA, xlim = c(14,135), ylim = c(0,0.5), xlab="Days Alive", ylab="Hazard", cex.axis=axis.text, cex.lab=axis.text)
lines(hazards.NFA.GOMP.f1$est~hazards.NFA.GOMP.f1$time, col="#ffcb91", lwd=2) 
lines(hazards.NFB.GOMP.f1$est~hazards.NFB.GOMP.f1$time, col="#e75758", lwd=2)
lines(hazards.NFC.GOMP.f1$est~hazards.NFC.GOMP.f1$time, col="#8b0000", lwd=2) 
lines(hazards.LFA.GOMP.f1$est~hazards.LFA.GOMP.f1$time, col="#00bfff", lwd=2) 
lines(hazards.LFB.GOMP.f1$est~hazards.LFB.GOMP.f1$time, col="#2b75cb", lwd=2) 
lines(hazards.LFC.GOMP.f1$est~hazards.LFC.GOMP.f1$time, col="#1c2c99", lwd=2) 


par(fig = c(0.1,0.6, 0.5, 0.95), new = T)  

plot(NA,NA, xlim = c(14,135), ylim = c(0,1), xlab="Days Alive", ylab="Survival", cex.axis=axis.text, cex.lab=axis.text)
lines(surv.NFA.GOMP.f1$est~surv.NFA.GOMP.f1$time, col="#ffcb91", lwd=2) 
lines(surv.NFB.GOMP.f1$est~surv.NFB.GOMP.f1$time, col="#e75758", lwd=2)
lines(surv.NFC.GOMP.f1$est~surv.NFC.GOMP.f1$time, col="#8b0000", lwd=2) 
lines(surv.LFA.GOMP.f1$est~surv.LFA.GOMP.f1$time, col="#00bfff", lwd=2) 
lines(surv.LFB.GOMP.f1$est~surv.LFB.GOMP.f1$time, col="#2b75cb", lwd=2) 
lines(surv.LFC.GOMP.f1$est~surv.LFC.GOMP.f1$time, col="#1c2c99", lwd=2) 

add_legend("top", legend=c(c("LFA","LFB","LFC","NFA","NFB","NFC")), lty=1, lwd=2,
           col=c("#00bfff","#2b75cb","#1c2c99","#ffcb91", "#e75758","#8b0000" ),
           horiz=TRUE, bty='n', cex=0.8)

dev.copy(pdf, "gomp.mortality.f1.pdf", height = 6, width = 6)
dev.off()
graphics.off()


#### Body Size ####

mat.effect.data.body<-mat.effect.data[-which(is.na(mat.effect.data$body.size)==TRUE),]
mat.effect.data.body$Treatment<-paste(mat.effect.data.body$Maternal.Food, mat.effect.data.body$clutch, sep="")

####Analysis####

#ANOVA

body.sizem1<-lm(aov(body.size ~ Maternal.Food + clutch + Maternal.Food*clutch, data=mat.effect.data.body))
body.sizem2<-aov(body.size ~ clutch + Maternal.Food + clutch*Maternal.Food, data=mat.effect.data.body)
summary(body.sizem2)
summary(body.sizem1)
sjt.lm(body.sizem1)


#no difference between m1 or m2

TukeyHSD(body.sizem2)

#make a column of residuals in the mat.effect.data df

mat.effect.data.body$bs.residuals <-residuals(body.sizem1)

#make a histogram of the residuals to assess distriution 
hist(mat.effect.data.body$bs.residuals)

#use plot function to get the Q-Q plot
plot(body.sizem1)
#make a standard error for body size
body.size.SE<-as.data.frame(as.list(aggregate(body.size ~ Treatment, data = mat.effect.data.body, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))
#use the standard error in the plot

food <- "Maternal.Food"
body.size.SE[ , food] <- NA
body.size.SE[1:3,4]<-"low"
body.size.SE[4:6,4]<-"normal"
body.size.SE$mat.food<-as.factor(body.size.SE$Maternal.Food)

str(mat.effect.data.body$Maternal.Food)
str(body.size.SE$Maternal.Food)
#### dot plot ####
body.size.SE<-as.data.frame(as.list(aggregate(body.size ~ clutch + Maternal.Food, data = mat.effect.data.body, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))
body.size.SE$Treatment <- paste(body.size.SE$clutch, body.size.SE$Maternal.food, sep="")

mat.effect.data.body$clutch <- revalue(mat.effect.data.body$clutch, c("A"="1","B"="2","C"="Old"))
body.size.SE$clutch <- revalue(body.size.SE$clutch, c("A"="1","B"="2","C"="Old"))

bs.point <- ggplot(mat.effect.data.body, aes(x = factor(clutch), y = body.size, colour = Maternal.Food)) +
  geom_jitter(aes(color = Maternal.Food), position = position_jitter(width = 0.05), alpha = 0.9, shape = 19)+
  geom_errorbar(data=body.size.SE, aes(x = clutch, ymin=body.size.mean-1.96*body.size.standerror,ymax=body.size.mean+1.96*body.size.standerror), width = 0.1,inherit.aes = FALSE, position=position_dodge(w=0.75))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=14))+
  scale_colour_manual(values=c("#FFA07A","#DC143C"),
                      labels = c("DR", "Normal")) +
  ylab("Body Size (mm)") + xlab("Clutch")+
  scale_y_continuous(limits = c(0,NA))

bs.point
# This will save a 400x400 file at 100 ppi
ggsave("bs.point.png", width=8, height=6, dpi=300)

#### Reproduction Data Manipulation Offspring ####

per.clutch<-melt(mat.effect.data.SA, id.vars = c("uniqueID","Maternal.Food", "clutch" ), measure.vars = c("c1n","c2n","c3n","c4n","c5n","c6n","c7n","c8n","c9n","c10n","c11n","c12n","c13n","c14n","c15n","c16n","c17n","c18n","c19n","c20n") )
per.clutch$Treatment<-as.factor(paste(per.clutch$Maternal.Food, per.clutch$clutch, sep =""))
levels(per.clutch$variable)<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")

colnames(per.clutch)[which(names(per.clutch) == "clutch")] <- "mat.age"
colnames(per.clutch)[which(names(per.clutch) == "variable")] <- "clutch"
colnames(per.clutch)[which(names(per.clutch) == "value")] <- "offspring"
per.clutch$uniqueID<-as.factor(per.clutch$uniqueID)

per.clutch$clutch<-as.integer(per.clutch$clutch)
str(per.clutch)

#make a new column with clutch scaled
#use this scaled clutch number in the polynomial analysis
#this should help with collinearity
per.clutch$clutchscale <- scale(per.clutch$clutch)

#raise the scaled clutch number to the power and then scale that power.
per.clutch$square<-per.clutch[,"clutchscale"]^2

#now scale this squared term
per.clutch$square.scale<-scale(per.clutch$square)

#do the same for a cube term
per.clutch$cube<-per.clutch[,"clutchscale"]^3
per.clutch$cube.scale<-scale(per.clutch$cube)

#### Plot ####
offspring.repro.SE<-as.data.frame(as.list(aggregate(offspring ~ Maternal.Food + clutch + mat.age, data = per.clutch, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))

offspring.pc <- ggplot(data = per.clutch, aes(x = clutch, y = offspring, colour = Maternal.Food, shape = mat.age)) +
  geom_point(position =position_dodge(w=0),stat="summary", fun.y="mean", size = 4) +
  geom_errorbar(data=offspring.repro.SE, aes(color= Maternal.Food, x = clutch, ymin=offspring.mean-1.96*offspring.standerror,ymax=offspring.mean+1.96*offspring.standerror),width = 0.1,inherit.aes = FALSE, position=position_dodge(w=0))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text( size = 14),
        axis.title.y = element_text( size = 14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))+
  scale_colour_manual(values=c("#FFA07A","#DC143C"),
                      labels = c("DR", "Normal")) +
  scale_shape_manual(values=c(19,17,11),name = "Maternal Age", labels = c("1st Clutch", "2nd Clutch", "Old"))+
  scale_x_continuous(name="Clutch", limits=c(0, 20), breaks=c(1:20))+
  scale_y_continuous(name="Clutch", limits=c(0, 18), breaks=c(1:20))
offspring.pc


# This will save a 400x400 file at 100 ppi
ggsave("per.clutch.plot.5.3.20.png", width=10, height=6, dpi=300)

####Analysis####

hist(per.clutch$offspring)

####polynomial per clutch####
clutch.offspring.m1<-lmer(offspring ~ Maternal.Food + mat.age + clutchscale + square.scale + cube.scale + Maternal.Food*clutchscale + Maternal.Food*cube.scale + Maternal.Food * square.scale + mat.age*clutchscale + mat.age* square.scale + mat.age* cube.scale+ (1|uniqueID),  data = per.clutch)
summary(clutch.offspring.m1)
qqnorm(resid(clutch.offspring.m1))
qqline(resid(clutch.offspring.m1))
anova(clutch.offspring.m1, test = "Chisq")

#remove the least significant highest order term
#remove Maternal.Food*cube
clutch.offspring.m1.1<-lmer(offspring ~ Maternal.Food + mat.age + clutchscale + square.scale + cube.scale + Maternal.Food*clutchscale + Maternal.Food * square.scale + mat.age*clutchscale + mat.age* square.scale + mat.age* cube.scale+ (1|uniqueID),  data = per.clutch)
summary(clutch.offspring.m1.1)
qqnorm(resid(clutch.offspring.m1.1))
qqline(resid(clutch.offspring.m1.1))
anova(clutch.offspring.m1.1, test = "Chisq")

#remove Maternal.Food*square
clutch.offspring.m1.2<-lmer(offspring ~ Maternal.Food + mat.age + clutchscale + square.scale + cube.scale + Maternal.Food*clutchscale +  mat.age*clutchscale + mat.age* square.scale + mat.age* cube.scale+ (1|uniqueID),  data = per.clutch)
summary(clutch.offspring.m1.2)
qqnorm(resid(clutch.offspring.m1.2))
qqline(resid(clutch.offspring.m1.2))
anova(clutch.offspring.m1.2, test = "Chisq")

#remove mat.age*cube
clutch.offspring.m1.3<-lmer(offspring ~ Maternal.Food + mat.age + clutchscale + square.scale + cube.scale + Maternal.Food*clutchscale +  mat.age*clutchscale + mat.age* square.scale + (1|uniqueID),  data = per.clutch)
summary(clutch.offspring.m1.3)
qqnorm(resid(clutch.offspring.m1.3))
qqline(resid(clutch.offspring.m1.3))
anova(clutch.offspring.m1.3, test = "Chisq")

#remove mat.age*square
clutch.offspring.m1.4<-lmer(offspring ~ Maternal.Food + mat.age + clutchscale + square.scale + cube.scale + Maternal.Food*clutchscale +  mat.age*clutchscale + (1|uniqueID),  data = per.clutch)
summary(clutch.offspring.m1.4)
qqnorm(resid(clutch.offspring.m1.4))
qqline(resid(clutch.offspring.m1.4))
anova(clutch.offspring.m1.4, test = "Chisq")


#remove the squared and cubed terms this is not quadratic

clutch.offspring.m1.5<-lmer(offspring ~  Maternal.Food - 1 + mat.age + clutchscale + Maternal.Food*clutchscale +  mat.age*clutchscale + (1|uniqueID),  data = per.clutch)
summary(clutch.offspring.m1.5)
qqnorm(resid(clutch.offspring.m1.5))
qqline(resid(clutch.offspring.m1.5))
anova(clutch.offspring.m1.5, test = "Chisq")

# do this again allowing for the intercept

clutch.offspring.m1.5<-lmer(offspring ~  Maternal.Food + mat.age + clutchscale + Maternal.Food*clutchscale +  mat.age*clutchscale + (1|uniqueID),  data = per.clutch)
summary(clutch.offspring.m1.5)
qqnorm(resid(clutch.offspring.m1.5))
qqline(resid(clutch.offspring.m1.5))
anova(clutch.offspring.m1.5, test = "Chisq")

####old mums####

####data manipulation####
old.mums<-read.csv("mat.sen.old.mums.csv")
old.mums.per.clutch<-melt(old.mums, id.vars = c("uniqueID","Food" ), measure.vars = c("c1n","c2n","c3n","c4n","c5n","c6n","c7n","c8n","c9n","c10n","c11n","c12n") )
colnames(old.mums.per.clutch)[which(names(old.mums.per.clutch) == "value")] <- "offspring"
colnames(old.mums.per.clutch)[which(names(old.mums.per.clutch) == "variable")] <- "Clutch"
levels(old.mums.per.clutch$Clutch)<-c("1","2","3","4","5","6","7","8","9","10","11","12")
old.mums.per.clutch$Clutch<-as.integer(old.mums.per.clutch$Clutch)
old.mums.per.clutch$uniqueID<-as.factor(old.mums.per.clutch$uniqueID)
str(old.mums.per.clutch)

hist(old.mums.per.clutch$offspring)

old.mums.per.clutch$clutchscale <- scale(old.mums.per.clutch$Clutch)
#raise the scaled clutch number to the power and then scale that power.
old.mums.per.clutch$square<-old.mums.per.clutch[,"clutchscale"]^2
#now scale this squared term
old.mums.per.clutch$square.scale<-scale(old.mums.per.clutch$square)
#do the same for a cube term
old.mums.per.clutch$cube<-old.mums.per.clutch[,"clutchscale"]^3
old.mums.per.clutch$cube.scale<-scale(old.mums.per.clutch$cube)

hist(old.mums.per.clutch$offspring)

old.mums.per.clutch$key <- "Data"

#### per clutch poisson #### 
hist(old.mums.per.clutch$offspring)
old.mums.offspring.m1<-glmer(offspring ~ Food + clutchscale + (1|uniqueID) + Food*clutchscale,  family = poisson, data = old.mums.per.clutch)
summary(old.mums.offspring.m1)
anova(old.mums.offspring.m1, test="Chisq")
qqnorm(resid(old.mums.offspring.m1))
qqline(resid(old.mums.offspring.m1))
plot(resid(old.mums.offspring.m1))

#### per clutch polynomial ####

old.mums.offspring.m2<-lmer(offspring ~ Food + clutchscale + square.scale + cube.scale + (1|uniqueID) + Food*clutchscale + Food*square.scale + Food*cube.scale , data = old.mums.per.clutch)
summary(old.mums.offspring.m2)
anova(old.mums.offspring.m2, test="Chisq")
qqnorm(resid(old.mums.offspring.m2))
qqline(resid(old.mums.offspring.m2))


#sequentially remove the non-significant interactions
#begin with cube*food and cube

old.mums.offspring.m3<-lmer(offspring ~ Food + clutchscale + square.scale + (1|uniqueID) + Food*clutchscale*square.scale , data = old.mums.per.clutch)
summary(old.mums.offspring.m3)
old.mums.anova<-anova(old.mums.offspring.m3, test="Chisq")
qqnorm(resid(old.mums.offspring.m3))
qqline(resid(old.mums.offspring.m3))
anova(old.mums.offspring.m3, test = "Chisq")

interplot(old.mums.offspring.m3, var1 = "Food", var2 = "clutchscale")
interplot(old.mums.offspring.m3, var1 = "Food", var2 = "square.scale", point = T)

anova(old.mums.offspring.m1, old.mums.offspring.m3)

plot(resid(old.mums.offspring.m3))

fit.old.mums.per.clutch <- data.frame(predict(old.mums.offspring.m3))

fit.old.mums.per.clutch$key <- "Model"

old.mums.per.clutch <- old.mums.per.clutch[-which(is.na(old.mums.per.clutch$offspring)==TRUE),]
old.mums.per.clutch.edit <- cbind(old.mums.per.clutch[,1:4], old.mums.per.clutch[,10])
colnames(old.mums.per.clutch.edit) <- c("uniqueID", "Food", "Clutch", "offspring", "key")

data <- data.frame(cbind(old.mums.per.clutch$uniqueID, old.mums.per.clutch$Food, old.mums.per.clutch$Clutch, fit.old.mums.per.clutch))

colnames(data) <- c("uniqueID", "Food", "Clutch", "offspring", "key")

data <- rbind(old.mums.per.clutch.edit, data)

library(plyr)
data$Food <- revalue(data$Food, c("LF"="DR", "NF"="Normal"))

data.data <- subset(data, key=="Data")
data.model <- subset(data, key=="Model")

mum.repro.SE<-as.data.frame(as.list(aggregate(offspring ~ Food + Clutch + key, data = data, FUN=function(x) c(mean = mean(x), standerror = sd(x)/sqrt(length(x))))))

ompc <- ggplot(data = data, aes(x = Clutch, y = offspring, colour = Food, shape = key)) +
  geom_point(position =position_dodge(w=0),stat="summary", fun.y="mean", size = 4) +
  geom_errorbar(data=mum.repro.SE, aes(color= Food, x = Clutch, ymin=offspring.mean-1.96*offspring.standerror,ymax=offspring.mean+1.96*offspring.standerror),width = 0.1,inherit.aes = FALSE, position=position_dodge(w=0))+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text( size = 14),
        axis.title.y = element_text( size = 14),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))+
  scale_colour_manual(values=c("#FFA07A","#DC143C"),
                      labels = c("DR", "Normal")) +
  scale_shape_manual(values=c(19,2),name = "Key", labels = c("Data", "Model"))+
  scale_x_continuous(name="Clutch", limits=c(0, 13), breaks=c(1:12))

ompc

ggsave("ompc.png", width=8, height=6, dpi=300)
#### Experiment panel figure ####

panel <-  

# not poly - larger AIC go with poly 
ARmodel.offspring <- lme(offspring ~ Food*clutchscale, 
                         random = ~ 1|uniqueID, 
                         correlation = corAR1(),
                         data = old.mums.per.clutch)

summary(ARmodel.offspring)
plot(resid(ARmodel.offspring))
anova(ARmodel.offspring, test = "Chisq")

# poly 
ARmodel.offspring.poly <- lme(offspring ~ Food*clutchscale + Food*square.scale, 
                         random = ~ 1|uniqueID, 
                         correlation = corAR1(),
                         data = old.mums.per.clutch)

summary(ARmodel.offspring.poly)

qqnorm(resid(ARmodel.offspring.poly))
qqline(resid(ARmodel.offspring.poly))
plot(resid(ARmodel.offspring.poly))

plot(ACF(ARmodel.offspring.poly, resType = "normalized", alpha=0.5))

# no AR 
NOARmodel.offspring.poly <- lme(offspring ~ Food*clutchscale + Food*square.scale, 
                              random = ~ 1|uniqueID,
                              data = old.mums.per.clutch)

summary(NOARmodel.offspring.poly)

plot(resid(NOARmodel.offspring.poly))

plot(ACF(NOARmodel.offspring.poly, resType = "normalized", alpha=0.5))


anova(ARmodel.offspring.poly, NOARmodel.offspring.poly)

anova(ARmodel.offspring, ARmodel.offspring.poly)
#### === === === === === === === === === === === === === ===  #### 
#### === === === ===       Epi Model     === === === === ===  #### 
#### === === === === === === === === === === === === === ===  #### 

require(fields)

#### === === === === === === === === === === === === === ===  #### 
#### === === === ===     Main Effects    === === === === ===  #### 
#### === === === === === === === === === === === === === ===  #### 

# all rates are per day
r <- 2.5
k <- r # max pop density
m <- 0.5 # maturation rate (average age at maturation is 10 days)
a <- 1.5
c <- 0.1 
B_baseline <- 0.9

low_d <- 0.1
med_d <- 0.5
high_d <- 1
from=-1
to = 1

mat_length = 501

#### Dr0Dx ####

x <- seq(from = from, to = to, length.out = mat_length)
y <- 0
mat<-0 # maternal effect on susceptibility
age<-0 # ageing effect on susceptibility

dr0x.main.low <- matrix( nrow = length(x),ncol = 4)
dr0x.main.med <- matrix( nrow = length(x),ncol = 4)
dr0x.main.high <- matrix( nrow = length(x),ncol = 4)

#### Low Mortality ####

d <- low_d

dr0x.main.low[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0x.main.low[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0x.main.low[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0x.main.low[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou
dr0x.main.low[dr0x.main.low<0]<-0

total <- rowSums(dr0x.main.low)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(x)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0x.yyu.low <- (dr0x.main.low[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0x.oyu.low <- ((dr0x.main.low[,3]/total)+(dr0x.main.low[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0x.you.low <- (dr0x.main.low[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0x.oou.low <- ((dr0x.main.low[,4]/total)+(dr0x.main.low[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0x.low.new <- rowSums((dr0x.main.low*weights)*(r0x.yyu.low + r0x.oyu.low + r0x.you.low + r0x.oou.low))

#### Medium Mortality ####

d <- med_d

dr0x.main.med[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0x.main.med[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0x.main.med[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0x.main.med[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0x.main.med[dr0x.main.med<0]<-0

total <- rowSums(dr0x.main.med)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(x)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0x.yyu.med <- (dr0x.main.med[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0x.oyu.med <- ((dr0x.main.med[,3]/total)+(dr0x.main.med[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0x.you.med <- (dr0x.main.med[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0x.oou.med <- ((dr0x.main.med[,4]/total)+(dr0x.main.med[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0x.med.new <- rowSums((dr0x.main.med*weights)*(r0x.yyu.med + r0x.oyu.med + r0x.you.med + r0x.oou.med))

#### High Mortality ####

d <- high_d

dr0x.main.high[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0x.main.high[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0x.main.high[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0x.main.high[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0x.main.high[dr0x.main.high<0]<-0

total <- rowSums(dr0x.main.high)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(x)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0x.yyu.high <- (dr0x.main.high[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0x.oyu.high <- ((dr0x.main.high[,3]/total)+(dr0x.main.high[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0x.you.high <- (dr0x.main.high[,2]/total)*(1/(m+a+c+d*(1-y)))
r0x.oou.high <- ((dr0x.main.high[,4]/total)+(dr0x.main.high[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0x.high.new <- rowSums((dr0x.main.high*weights)*(r0x.yyu.high + r0x.oyu.high + r0x.you.high + r0x.oou.high))

#### Dr0Dy ####

y <- seq(from = from, to = to, length.out = mat_length)
x <- 0
mat<-0 # maternal effect on susceptibility
age<-0 # ageing effect on susceptibility

dr0y.main.low <- matrix( nrow = length(y),ncol = 4)
dr0y.main.med <- matrix( nrow = length(y),ncol = 4)
dr0y.main.high <- matrix( nrow = length(y),ncol = 4)

#### Low Mortality ####

d <- low_d

dr0y.main.low[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0y.main.low[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0y.main.low[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0y.main.low[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0y.main.low[dr0y.main.low<0]<-0

total <- rowSums(dr0y.main.low)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(y)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0y.yyu.low <- (dr0y.main.low[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oyu.low <- ((dr0y.main.low[,3]/total)+(dr0y.main.low[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0y.you.low <- (dr0y.main.low[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oou.low <- ((dr0y.main.low[,4]/total)+(dr0y.main.low[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0y.low.new <- rowSums((dr0y.main.low*weights)*(r0y.yyu.low + r0y.oyu.low + r0y.you.low + r0y.oou.low))

#### Medium Mortality ####

d <- med_d

dr0y.main.med[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0y.main.med[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0y.main.med[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0y.main.med[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0y.main.med[dr0y.main.med<0]<-0

total <- rowSums(dr0y.main.med)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(y)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0y.yyu.med <- (dr0y.main.med[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oyu.med <- ((dr0y.main.med[,3]/total)+(dr0y.main.med[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0y.you.med <- (dr0y.main.med[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oou.med <- ((dr0y.main.med[,4]/total)+(dr0y.main.med[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0y.med.new <- rowSums((dr0y.main.med*weights)*(r0y.yyu.med + r0y.oyu.med + r0y.you.med + r0y.oou.med))

# High Mortality

d <- high_d

dr0y.main.high[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0y.main.high[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0y.main.high[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0y.main.high[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0y.main.high[dr0y.main.high<0]<-0

total <- rowSums(dr0y.main.high)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(y)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0y.yyu.high <- (dr0y.main.high[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oyu.high <- ((dr0y.main.high[,3]/total)+(dr0y.main.high[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0y.you.high <- (dr0y.main.high[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oou.high <- ((dr0y.main.high[,4]/total)+(dr0y.main.high[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0y.high.new <- rowSums((dr0y.main.high*weights)*(r0y.yyu.high + r0y.oyu.high + r0y.you.high + r0y.oou.high))

#### dr0DA ####

y <- 0
x <- 0
age<- seq(from = from, to = to, length.out = mat_length)  
mat<- 0

dr0A.main.low <- matrix( nrow = length(age),ncol = 4)
dr0A.main.med <- matrix( nrow = length(age),ncol = 4)
dr0A.main.high <- matrix( nrow = length(age),ncol = 4)

# Low Mortality 
d <- low_d

dr0A.main.low[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0A.main.low[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0A.main.low[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0A.main.low[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0A.main.low[dr0A.main.low<0]<-0

total <- rowSums(dr0A.main.low)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(age)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0A.yyu.low <- (dr0A.main.low[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oyu.low <- ((dr0A.main.low[,3]/total)+(dr0A.main.low[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0A.you.low <- (dr0A.main.low[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oou.low <- ((dr0A.main.low[,4]/total)+(dr0A.main.low[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0A.low.new <- rowSums((dr0A.main.low*weights)*(r0A.yyu.low + r0A.oyu.low + r0A.you.low + r0A.oou.low))

# Medium Mortality

d <- med_d

dr0A.main.med[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0A.main.med[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0A.main.med[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0A.main.med[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0A.main.med[dr0A.main.med<0]<-0

total <- rowSums(dr0A.main.med)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(age)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0A.medmort <- rowSums(dr0A.main.med*weights)/(a+d)

r0A.yyu.med <- (dr0A.main.med[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oyu.med <- ((dr0A.main.med[,3]/total)+(dr0A.main.med[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0A.you.med <- (dr0A.main.med[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oou.med <- ((dr0A.main.med[,4]/total)+(dr0A.main.med[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0A.med.new <- rowSums((dr0A.main.med*weights)*(r0A.yyu.med + r0A.oyu.med + r0A.you.med + r0A.oou.med))

# High Mortality

d <- high_d

dr0A.main.high[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0A.main.high[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0A.main.high[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0A.main.high[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0A.main.high[dr0A.main.high<0]<-0

total <- rowSums(dr0A.main.high)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(age)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0A.yyu.high <- (dr0A.main.high[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oyu.high <- ((dr0A.main.high[,3]/total)+(dr0A.main.high[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0A.you.high <- (dr0A.main.high[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oou.high <- ((dr0A.main.high[,4]/total)+(dr0A.main.high[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0A.high.new <- rowSums((dr0A.main.high*weights)*(r0A.yyu.high + r0A.oyu.high + r0A.you.high + r0A.oou.high))

#### Dr0DM ####

y <- 0
x <- 0
age<- 0
mat<- seq(from = from, to = to, length.out = mat_length)

dr0M.main.low <- matrix( nrow = length(mat),ncol = 4)
dr0M.main.med <- matrix( nrow = length(mat),ncol = 4)
dr0M.main.high <- matrix( nrow = length(mat),ncol = 4)

# Low Mortality 
d <- low_d

dr0M.main.low[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0M.main.low[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0M.main.low[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0M.main.low[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0M.main.low[dr0M.main.low<0]<-0

total <- rowSums(dr0M.main.low)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(mat)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0M.yyu.low <- (dr0M.main.low[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oyu.low <- ((dr0M.main.low[,3]/total)+(dr0M.main.low[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0M.you.low <- (dr0M.main.low[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oou.low <- ((dr0M.main.low[,4]/total)+(dr0M.main.low[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0M.low.new <- rowSums((dr0M.main.low*weights)*(r0M.yyu.low + r0M.oyu.low + r0M.you.low + r0M.oou.low))

# Medium Mortality

d <- med_d

dr0M.main.med[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0M.main.med[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0M.main.med[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0M.main.med[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0M.main.med[dr0M.main.med<0]<-0

total <- rowSums(dr0M.main.med)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(mat)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0M.yyu.med <- (dr0M.main.med[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oyu.med <- ((dr0M.main.med[,3]/total)+(dr0M.main.med[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0M.you.med <- (dr0M.main.med[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oou.med <- ((dr0M.main.med[,4]/total)+(dr0M.main.med[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0M.med.new <- rowSums((dr0M.main.med*weights)*(r0M.yyu.med + r0M.oyu.med + r0M.you.med + r0M.oou.med))

# High Mortality

d <- high_d

dr0M.main.high[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0M.main.high[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0M.main.high[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0M.main.high[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0M.main.high[dr0M.main.high<0]<-0

total <- rowSums(dr0M.main.high)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(mat)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0M.yyu.high <- (dr0M.main.high[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oyu.high <- ((dr0M.main.high[,3]/total)+(dr0M.main.high[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0M.you.high <- (dr0M.main.high[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oou.high <- ((dr0M.main.high[,4]/total)+(dr0M.main.high[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0M.high.new <- rowSums((dr0M.main.high*weights)*(r0M.yyu.high + r0M.oyu.high + r0M.you.high + r0M.oou.high))

#### Density ####

#### dnudx ####

y <- 0
x <- seq(from = -1, to = 1, length.out = mat_length)
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
y <- seq(from = -1, to = 1, length.out = mat_length)

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

#### Main effects and Density panel plot in one ####

x <- seq(from = from, to = to, length.out = mat_length)
y <- seq(from = from, to = to, length.out = mat_length)
A <- seq(from = from, to = to, length.out = mat_length)
M <- seq(from = from, to = to, length.out = mat_length)

axis.text <- 2

dev.new(height = 15, width = 6)
par(mar=c(4,5,3,3))
par(mfrow=c(5,2), mgp=c(3,1,0))

plot(NA,NA, xlim = c(-1,1), ylim = c(0,2.5), xlab="Resistance", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
#title(outer=TRUE,adj=0.115,main="A",cex=1, col="black",font=1,line=-3.8)
lines(r0A.low.new~A, col="#8B0000", lwd=2) # low dark red
lines(r0A.med.new~A, col="#EE0000", lwd=2) 
lines(r0A.high.new~A, col="#E9967A", lwd=2) # high salmon

plot(NA,NA, xlim = c(-1,1), ylim = c(0,2.5), xlab="Transgenerational Resistance", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
#title(outer=TRUE,adj=0.625,main="B",cex=1, col="black",font=1,line=-3.8)
lines(r0M.low.new~M, col="#8B0000", lwd=2) # low dark red
lines(r0M.med.new~M, col="#EE0000", lwd=2)
lines(r0M.high.new~M, col="#E9967A", lwd=2) # high salmon

plot(NA,NA, xlim = c(-1,1), ylim = c(0,1.5), xlab="Reproduction", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
#title(outer=TRUE,adj=0.115,main="C",cex=1, col="black",font=1,line=-21.7)
lines(r0x.low.new~x, col="#8B0000", lwd=2) # low dark red
lines(r0x.med.new~x, col="#EE0000", lwd=2)
lines(r0x.high.new~x, col="#E9967A", lwd=2) # high salmon

plot(NA,NA, xlim = c(-1,1), ylim = c(0,1.5), xlab="Mortality", ylab="Transmission Potential", cex.axis=axis.text, cex.lab=axis.text)
#title(outer=TRUE,adj=0.625,main="D",cex=1, col="black",font=1,line=-21.7)
lines(r0y.low.new~y, col="#8B0000", lwd=2) # low dark red
lines(r0y.med.new~y, col="#EE0000", lwd=2)
lines(r0y.high.new~y, col="#E9967A", lwd=2) # high salmon

plot(NA,NA,xlim=c(-1,1),ylim=c(0,2.5), xlab = " Reproduction ", ylab="Total density", cex.lab=axis.text, cex.axis=axis.text)
lines(dnudx[,1]~x,lwd=2, lty = 1,col="firebrick1")

lines(dnudx[,2]~x,lwd=2, lty = 1,col="dodgerblue1")#you

lines(dnudx[,3]~x,lwd=2, lty=1,col="goldenrod1")

lines(dnudx[,4]~x,lwd=2, lty=1,col="forestgreen")

lines(dnudx[,1]+dnudx[,2]+dnudx[,3]+dnudx[,4]~x, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,2.5), xlab = "Mortality ", ylab="Total density",cex.lab=axis.text, cex.axis=axis.text)
lines(dnudy[,1]~y,lwd=2, lty = 1,col="firebrick1")

lines((dnudy[,2]-0.02)~y,lwd=3, lty = 1,col="dodgerblue1")#you

lines(dnudy[,3]~y,lwd=2, lty=1,col="goldenrod1")

lines(dnudy[,4]~y,lwd=2, lty=1,col="forestgreen")

lines(dnudy[,1]+dnudy[,2]+dnudy[,3]+dnudy[,4]~y, lwd=2, lty=1,col="blueviolet")


plot(NA,NA,xlim=c(-1,1),ylim=c(0,2.5), xlab = " Reproduction ", ylab="Total density", cex.lab=axis.text, cex.axis=axis.text)
lines((dnudx2[,1]-0.02)~x,lwd=3, lty = 1,col="firebrick1")

lines((dnudx2[,2]-0.02)~x,lwd=2, lty = 1,col="dodgerblue1")#you

lines(dnudx2[,3]~x,lwd=2, lty=1,col="goldenrod1")

lines(dnudx2[,4]~x,lwd=2, lty=1,col="forestgreen")

lines(dnudx2[,1]+dnudx2[,2]+dnudx2[,3]+dnudx2[,4]~x, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,2.5), xlab = " Mortality ", ylab="Total density",cex.lab=axis.text, cex.axis=axis.text)
lines(dnudy2[,1]~y,lwd=2, lty = 1,col="firebrick1")

lines((dnudy2[,2]-0.02)~y,lwd=3, lty = 1,col="dodgerblue1")#you

lines(dnudy2[,3]~y,lwd=2, lty=1,col="goldenrod1")

lines(dnudy2[,4]~y,lwd=2, lty=1,col="forestgreen")

lines(dnudy2[,1]+dnudy2[,2]+dnudy2[,3]+dnudy2[,4]~y, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,2.5), xlab = " Reproduction ", ylab="Total density", cex.lab=axis.text, cex.axis=axis.text)
lines(dnudx3[,1]~x,lwd=2, lty = 1,col="firebrick1")

lines(dnudx3[,2]~x,lwd=2, lty = 1,col="dodgerblue1")#you

lines(dnudx3[,3]~x,lwd=2, lty=1,col="goldenrod1")

lines(dnudx3[,4]~x,lwd=2, lty=1,col="forestgreen")

lines(dnudx3[,1]+dnudx3[,2]+dnudx3[,3]+dnudx3[,4]~x, lwd=2, lty=1,col="blueviolet")

plot(NA,NA,xlim=c(-1,1),ylim=c(0,2.5), xlab = " Mortality ", ylab="Total density",cex.lab=axis.text, cex.axis=axis.text)
lines(dnudy3[,1]~y,lwd=2, lty = 1,col="firebrick1")

lines((dnudy3[,2]-0.02)~y,lwd=3, lty = 1,col="dodgerblue1")#you

lines(dnudy3[,3]~y,lwd=2, lty=1,col="goldenrod1")

lines(dnudy3[,4]~y,lwd=2, lty=1,col="forestgreen")

lines(dnudy3[,1]+dnudy3[,2]+dnudy3[,3]+dnudy3[,4]~y, lwd=2, lty=1,col="blueviolet")


add_legend("bottom", legend=c(c("Uyy","Uyo","Uoy","Uoo","Total")), lty=1, lwd=2,
           col=c("firebrick1","dodgerblue1","goldenrod1","forestgreen", "blueviolet"),
           horiz=TRUE, bty='n', cex=1.2)

add_legend("top", legend=c(expression(paste(delta,"= 0.1")), expression(paste(delta,"= 0.5")), expression(paste(delta,"= 1.5"))), lty=1, lwd=2,
          col=c("#8B0000", "#EE0000", "#E9967A"),
          horiz=TRUE, bty='n', cex=1.2)

dev.copy(pdf, "main.effects.finalsizes.pdf", height = 15, width = 7)
dev.off()
graphics.off()

#### === === === === === === === === === === === === === ===  #### 
#### === === === ===     Interactions    === === === === ===  #### 
#### === === === === === === === === === === === === === ===  #### 

r <- 2.5
k <- r # max pop density
m <- 0.5 # maturation rate (average age at maturation is 10 days)
a <- 1.5
c <- 0.1
B_baseline <- 0.9

low_d <- 0.1 
med_d <- 0.5
high_d <- 1
from=-1
to = 1

mat_length = 501

library(RColorBrewer)

#### === === === === === === === === === === === ===  #### 
#### === === === === === === A by M  === === === ===  #### 
#### === === === === === === === === === === === ===  #### 

#### Low Mortality ####

d <- low_d

age <- matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length)
y <- 0
mat <- t(matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length))
x <- 0

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r01 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### Medium Mortality ####

d <- med_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r02 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

min(r02)
max(r02)

r02[]<-min(r02)

#### High Mortality ###

d <- high_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r03 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### === === === === === === === === === === === ===  #### 
#### === === === === === === x by y  === === === ===  #### 
#### === === === === === === === === === === === ===  #### 

#### Low Mortality ####

d <- low_d

x <- matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length)
mat <- 0
y <- t(matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length))
age <- 0

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow = mat_length, ncol = mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r04 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### Medium Mortality ####

d <- med_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r05 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### High Mortality ###

d <- high_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r06 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### === === === === === === === === === === === ===  #### 
#### === === === === === === x by A  === === === ===  #### 
#### === === === === === === === === === === === ===  #### 

#### Low Mortality ####

d <- low_d

x <- matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length)
y <- 0
mat <- 0
age <- t(matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length))

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r07 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### Medium Mortality ####

d <- med_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r08 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

min(r08)
max(r08)

r08[]<-min(r08)

#### High Mortality ###

d <- high_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r09 <- ((yyu*BYy) + (you*BYo) + (oyu*BOy) + (oou*BOo))*(yyum+oyum+youm+ooum)

#### === === === === === === === === === === === ===  #### 
#### === === === === === === x by M  === === === ===  #### 
#### === === === === === === === === === === === ===  #### 

#### Low Mortality ####

d <- low_d

x <- matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length)
y <- 0
mat <- t(matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length))
age <- 0

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r010 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### Medium Mortality ####

d <- med_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r011 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### High Mortality ###

d <- high_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r012 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### === === === === === === === === === === === ===  #### 
#### === === === === === === y by A  === === === ===  #### 
#### === === === === === === === === === === === ===  #### 

#### Low Mortality ####

d <- low_d
y <- matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length)
x <- 0
age <- t(matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length))
mat <- 0

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r013 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### Medium Mortality ####

d <- med_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r014 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### High Mortality ###

d <- high_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r015 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### === === === === === === === === === === === ===  #### 
#### === === === === === === y by M  === === === ===  #### 
#### === === === === === === === === === === === ===  #### 

#### Low Mortality ####

d <- low_d

y <- matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length)
x <- 0
mat <- t(matrix(rep(seq(-1,1,length.out = mat_length),each=mat_length),nrow = mat_length))
age <- 0

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r016 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### Medium Mortality ####

d <- med_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r017 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

#### High Mortality ###

d <- high_d

yyu <- (d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
you <- (-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
oyu <- (d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
oou <- (k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

yyu[yyu<0] <- 0
you[you<0] <- 0
oyu[oyu<0] <- 0
oou[oou<0] <- 0

total <- yyu + you + oyu + oou

BYy <- matrix(B_baseline*(1-mat)*(1-age), nrow=mat_length, ncol=mat_length)
BYo <- matrix(B_baseline*(1+mat)*(1-age), nrow = mat_length, ncol = mat_length) # transmission to young individuals with old mothers
BOy <- matrix(B_baseline*(1-mat)*(1+age), nrow = mat_length, ncol = mat_length)# transmission to old individuals with young mothers
BOo <- matrix(B_baseline*(1+mat)*(1+age), nrow = mat_length, ncol = mat_length) # transmission to old individuals with old mothers

yyum <- (yyu/total)*(1/(m+a+c+(d*(1-y))))
oyum <- (((oyu/total)+(yyu/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))
youm <- (you/total)*(1/(m+a+c+(d*(1-y))))
ooum <- (((oou/total)+(you/total)*(m/(m+a+c+(d*(1-y))))))*(1/(a+c+(d*(1+y))))

r018 <- (yyu*BYy + you*BYo + oyu*BOy + oou*BOo)*(yyum+oyum+youm+ooum)

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
par(mar=c(3,4,1,5),mgp=c(2,1,0))
par(oma=c(1,1,1,1))

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r01)), xlab="Trans. Gen. Resis ", ylab="Resistance",col=cols,zlim=c(min(r01),max(r01)),useRaster=T, cex.lab = 1.5)
mtext( 'A', side=90, line=656, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r01)),zlim=c(min((r01),na.rm=T),max((r01),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r01),na.rm=T),max((r01),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r02)), xlab="Trans. Gen. Resis ", ylab="Resistance",col=cols,zlim=c(min(r02),max(r02)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r02)),zlim=c(min((r02),na.rm=T),max((r02),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(0,1), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r03)), xlab="Trans. Gen. Resis ", ylab="Resistance",col=cols,zlim=c(min(r03),max(r03)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r03)),zlim=c(min((r03),na.rm=T),max((r03),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r03),na.rm=T),max((r018),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r04)), xlab="Reproduction", ylab="Mortality",col=cols,zlim=c(min((r04),na.rm=T),max((r04),na.rm=T)),useRaster=T, cex.lab = 1.5)
mtext( 'B', side=90, line=541, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r04)),zlim=c(min((r04),na.rm=T),max((r04),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r04),na.rm=T),max((r04),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r05)), xlab="Reproduction", ylab="Mortality",col=cols,zlim=c(min((r05),na.rm=T),max((r05),na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r05)),zlim=c(min((r05),na.rm=T),max((r05),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r05),na.rm=T),max((r05),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r06)), xlab="Reproduction", ylab="Mortality ",col=cols,zlim=c(min((r06),na.rm=T),max((r06),na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r06)),zlim=c(min((r06),na.rm=T),max((r06),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r06),na.rm=T),max((r06),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r07)), xlab="Reproduction", ylab="Resistance",col=cols,zlim=c(min((r07),na.rm=T),max((r07),na.rm=T)),useRaster=T, cex.lab = 1.5)
mtext( 'C', side=90, line=422, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r07)),zlim=c(min((r07),na.rm=T),max((r07),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r07),na.rm=T),max((r07),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r08)), xlab="Reproduction", ylab="Resistance",col=cols,zlim=c(min(r08),max(r08)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r08)),zlim=c(min((r08),na.rm=T),max((r08),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(0,1), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r09)), xlab="Reproduction", ylab="Resistance",col=cols,zlim=c(min((r09), na.rm=T),max((r09), na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r09)),zlim=c(min((r09),na.rm=T),max((r09),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r09),na.rm=T),max((r09),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols, bias=2)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r010)), xlab="Reproduction", ylab="Trans. Gen. Resis",col=cols,zlim=c(min(r010),max(r010)),useRaster=T, cex.lab = 1.5)
mtext( 'D', side=90, line=310, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r010)),zlim=c(min((r010),na.rm=T),max((r010),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r010),na.rm=T),max((r010),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r011)), xlab="Reproduction", ylab="Trans. Gen. Resis",col=cols,zlim=c(min(r011),max(r011)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r011)),zlim=c(min((r011),na.rm=T),max((r011),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r011),na.rm=T),max((r011),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r012)), xlab="Reproduction", ylab="Trans. Gen. Resis",col=cols,zlim=c(min((r012), na.rm=T),max((r012), na.rm=T)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r012)),zlim=c(min((r012),na.rm=T),max((r012),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r012),na.rm=T),max((r012),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)


image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r013)), xlab="Mortality", ylab="Resistance",col=cols,zlim=c(min(r013),max(r013)),useRaster=T, cex.lab = 1.5)
mtext( 'E', side=90, line=190, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r013)),zlim=c(min((r013),na.rm=T),max((r013),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r013),na.rm=T),max((r013),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r014)), xlab="Mortality", ylab="Resistance",col=cols,zlim=c(min(r014),max(r014)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r014)),zlim=c(min((r014),na.rm=T),max((r014),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r014),na.rm=T),max((r014),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r015)), xlab="Mortality ", ylab="Resistance",col=cols,zlim=c(min(r015),max(r015)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r015)),zlim=c(min((r015),na.rm=T),max((r015),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r015),na.rm=T),max((r015),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r016)), xlab="Mortality ", ylab="Trans. Gen. Resis",col=cols,zlim=c(min(r016),max(r016)),useRaster=T, cex.lab = 1.5)
mtext( 'F', side=90, line=78, at=grconvertX(60,'npc','nic'), outer=TRUE )
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r016)),zlim=c(min((r016),na.rm=T),max((r016),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r016),na.rm=T),max((r016),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r017)), xlab="Mortality ", ylab="Trans. Gen. Resis",col=cols,zlim=c(min(r017),max(r017)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r017)),zlim=c(min((r017),na.rm=T),max((r017),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r017),na.rm=T),max((r017),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

image(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r018)), xlab="Mortality ", ylab="Trans. Gen. Resis",col=cols,zlim=c(min(r018),max(r018)),useRaster=T, cex.lab = 1.5)
contour(seq(-1,1, length.out = mat_length), seq(-1, 1, length.out = mat_length), z=t((r018)),zlim=c(min((r018),na.rm=T),max((r018),na.rm=T)),add=T,nlevels=nlevels,drawlabels=F,lwd=0.2)
image.plot(zlim=c(min((r018),na.rm=T),max((r018),na.rm=T)), nlevel=501,legend.only=TRUE, horizontal=F,col=cols)

dev.copy(pdf, "fig. 3 text sizes.pdf", height = 10, width = 8)
dev.off()
graphics.off()

#### === === === === === === === === === === === ===  #### 
#### === === === === === === Humans  === === === ===  #### 
#### === === === === === === === === === === === ===  ####

# The effects on pathogen transmission, of reducing human senescence  

r <- 2.5 # Reproduction 
k <- r # Maximum Population Density 
m <- 0.5 # Rate of Maturation 
a <- 1.5 # Virulence 
B_baseline <- 0.9 # Transmission 
mat_length = 501 # Vector Length 
from = 0.5
to = -1
d <- 0.1 # Always using low mortality 
c <- 0.1

# Reproductive Senescence

x <- seq(from = from, to = to, length.out = mat_length)
y <- 0.5
mat<- 0.5
age<- 0.5

dr0x.human <- matrix( nrow = length(x),ncol = 4)

dr0x.human[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0x.human[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0x.human[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0x.human[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0x.human[dr0x.human<0]<-0

total <- rowSums(dr0x.human)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(x)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0x.yyu.human <- (dr0x.human[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0x.oyu.human <- ((dr0x.human[,3]/total)+(dr0x.human[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0x.you.human <- (dr0x.human[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0x.oou.human <- ((dr0x.human[,4]/total)+(dr0x.human[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0x.human <- rowSums((dr0x.human*weights)*(r0x.yyu.human + r0x.oyu.human + r0x.you.human + r0x.oou.human))

# Mortality Senescence 

y <- seq(from = from, to = to, length.out = mat_length)
x <- 0.5
mat<-0.5 # maternal effect on susceptibility
age<-0.5 # ageing effect on susceptibility

dr0y.human <- matrix( nrow = length(y),ncol = 4)

dr0y.human[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0y.human[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0y.human[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0y.human[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0y.human[dr0y.human<0]<-0

total <- rowSums(dr0y.human)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(y)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0y.yyu.human <- (dr0y.human[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oyu.human <- ((dr0y.human[,3]/total)+(dr0y.human[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0y.you.human <- (dr0y.human[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0y.oou.human <- ((dr0y.human[,4]/total)+(dr0y.human[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0y.human <- rowSums((dr0y.human*weights)*(r0y.yyu.human + r0y.oyu.human + r0y.you.human + r0y.oou.human))


# Age Specific Resistance

y <- 0.5
x <- 0.5
age<- seq(from = from, to = to, length.out = mat_length)  
mat<- 0.5

dr0A.human <- matrix( nrow = length(age),ncol = 4)

dr0A.human[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0A.human[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0A.human[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0A.human[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0A.human[dr0A.human<0]<-0

total <- rowSums(dr0A.human)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(age)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0A.yyu.human <- (dr0A.human[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oyu.human <- ((dr0A.human[,3]/total)+(dr0A.human[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0A.you.human <- (dr0A.human[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0A.oou.human <- ((dr0A.human[,4]/total)+(dr0A.human[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0A.human <- rowSums((dr0A.human*weights)*(r0A.yyu.human + r0A.oyu.human + r0A.you.human + r0A.oou.human))


# Transgenerational Resistance

y <- 0.5
x <- 0.5
age<- 0.5
mat<- seq(from = from, to = to, length.out = mat_length)

dr0M.human <- matrix( nrow = length(mat),ncol = 4)

dr0M.human[,1]<-(d^2*k*(1+x)*(1+y)^2*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Yyu
dr0M.human[,2]<-(-((d*k*m*(-1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)))#You
dr0M.human[,3]<-(d*k*m*(1+x)*(1+y)*(-m*r*(-1+x)-d*(m-r*(1+x))*(1+y)+d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oyu
dr0M.human[,4]<-(k*m^2*(-1+x)*(m*r*(-1+x)+d*(m-r*(1+x))*(1+y)-d^2*(-1+y^2)))/(r*(d+m+d*y)*(m-m*x+d*(1+x)*(1+y))^2)#Oou

dr0M.human[dr0M.human<0]<-0

total <- rowSums(dr0M.human)

BYy <- B_baseline*(1-mat)*(1-age) # transmission to young individuals with young mothers
BYo <- B_baseline*(1+mat)*(1-age) # transmission to young individuals with old mothers
BOy <- B_baseline*(1-mat)*(1+age) # transmission to old individuals with young mothers
BOo <- B_baseline*(1+mat)*(1+age) # transmission to old individuals with old mothers

weights<-matrix(ncol=4,nrow=length(mat)) # weighting for transmission potential, we weight by susceptibility
weights[,1]<-BYy # weighting for young individuals from young mums
weights[,2]<-BYo # weighting for young individuals from old mums
weights[,3]<-BOy # weighting for old individuals from young mums
weights[,4]<-BOo # weighting for old individuals from old mums

r0M.yyu.human <- (dr0M.human[,1]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oyu.human <- ((dr0M.human[,3]/total)+(dr0M.human[,1]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))
r0M.you.human <- (dr0M.human[,2]/total)*(1/(m+a+c+(d*(1-y))))
r0M.oou.human <- ((dr0M.human[,4]/total)+(dr0M.human[,2]/total)*(m/(m+a+c+(d*(1-y)))))*(1/(a+c+(d*(1+y))))

r0M.human <- rowSums((dr0M.human*weights)*(r0M.yyu.human + r0M.oyu.human + r0M.you.human + r0M.oou.human))

#### Plot ####

#### Human Transmission Plot ####

from = 0
to = 1.5

x <- seq(from = 0, to = 1.5, length.out = mat_length)

dev.new(height = 4, width = 4)
par(mar=c(5,4,1,1))
par(mfrow=c(1,1), mgp=c(3,1,0))

plot(NA, NA, xlim = c(from,to), ylim = c(0,3.5), xlab = "Senescence", ylab = "Transmission Potential",cex.lab=1, cex.axis=1, bty="n")
#title(outer=TRUE,adj=0.23,main="A",cex=1, col="black",font=1,line=-1.85)
lines(r0x.human ~ x,lwd=2, lty = 1,col="#EBCC2A")
lines(r0y.human ~ x,lwd=2, lty = 1,col="#3B9AB2")
lines(r0A.human ~ x,lwd=2, lty = 1,col="#b2df8a")
lines(r0M.human ~ x,lwd=2, lty = 1,col="#e78ac3")


add_legend("top", legend=c("Reprodution", "Mortality", "Susceptibility", "Trans-gen Susceptibility"), lty=1, lwd=3,
           col=c("#EBCC2A", "#3B9AB2", "#b2df8a", "#e78ac3"),
           horiz=TRUE, bty='n', cex=0.7)

dev.copy(pdf, "slow.sen.3.5.20.pdf", height = 5.5, width = 6)
dev.off()
graphics.off()


