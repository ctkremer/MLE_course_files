
fit.timeFert<-mle2(log.biomass~dnorm(mean=b0+b1*centered.year,sd=exp(s0)),
                       parameters=list(b0~fertilized),
                       start=list(b0=4,b1=-0.05, s0=0.5),data=c1)

fit.time2Fert<-mle2(log.biomass~dnorm(mean=b0,sd=exp(s0)),
                        parameters=list(b0~fertilized+centered.year+I(centered.year^2)),
                        start=list(b0=4,s0=0.5),data=c1)

fit.time.var<-mle2(log.biomass~dnorm(mean=b0+b1*centered.year,sd=exp(s0+s1*centered.year)),start=list(b0=4,b1= -0.05, s0=0.5,s1=0),data=c1)

fit.timeFert.var<-mle2(log.biomass~dnorm(mean=b0,sd=exp(s0+s1*centered.year)),
                        parameters=list(b0~fertilized+centered.year),
                        start=list(b0=4,s0=0.5,s1=0),data=c1)

head(c1)
fit.time2Fert.var<-mle2(log.biomass~dnorm(mean=b0,sd=exp(s0+s1*centered.year)),
                   parameters=list(b0~fertilized+centered.year+I(centered.year^2)),
                   start=list(b0=4,s0=0.5,s1=0),data=c1)
summary(fit.time.var)
summary(fit.time2Fert.var)


fit.time2Fert.varDist<-mle2(log.biomass~dnorm(mean=b0,sd=exp(s0)),
                        parameters=list(b0~fertilized+centered.year+I(centered.year^2),
                                        s0~disturbed+centered.year),
                        start=list(b0=4,s0=0.5),data=c1)
summary(fit.time2Fert.varDist)


AICtab(fit.time.var,fit.timeFert.var,fit.timeFert,fit.time2Fert.var,fit.time2Fert,fit.time2Fert.varFert)

pd<-predict(fit.time2Fert.var)
pd<-data.frame(pd,centered.year=c1$centered.year,fertilized=c1$fertilized)
pd<-unique(pd)
head(pd)

ggplot(c1,aes(x=centered.year,y=log.biomass))+
  geom_point(aes(colour=fertilized))+
  geom_line(data=pd,aes(y=pd,colour=fertilized))

curve(exp(-0.44026043+-0.0360528*x)^2,-12,12,ylim=c(0,1.2),xlab='centered.year',ylab='var(log.biomass), residual')


ggplot(c1,aes(x=centered.year,y=log.biomass))+
  geom_point(aes(colour=fertilized))+
  facet_wrap(~disturbed,nrow=2)

plot(log.biomass~centered.year,data=c1,col=c1$fert.code+1,pch=c1$dist.code+1)


