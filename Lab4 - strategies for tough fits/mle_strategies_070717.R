#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Lab 4: Strategies for difficult fits with MLE
# ELME MLE course
# Colin T. Kremer

# updated 5/2015, 7/2017

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Dealing with difficult fits ####

# So what do you do when you're trying to use MLE, maybe for a complicated 
# or noisy data set, and you don't know what parameter guesses to try?

# Options:

# 1) Educated guesses/method of moments
# 2) Trial and error + plotting
# 3) Brute force (via grid.mle)

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Data setup ####

data<-read.csv("/Users/colin/Teaching/ELME/ELME 2017/labs/Lab4 - strategies for tough fits/thermal_tolerance.csv")

# This data set contains the exponential growth rate of a culture of 
# phytoplankton measured at different temperatures, sometimes called a 
# thermal tolerance curve or thermal reaction norm. These curves usually 
# are unimodal, and left skewed - definitely not linear! Sounds like a 
# great opportunity for using MLE.

# here's what this one looks like
plot(grths~temps,data=data)

# one of the many nonlinear functions proposed for the shape of these 
# curves comes from Norberg et al. 2004,
nbcurve<-function(x,opt,w,a,b){
  res<-a*exp(b*x)*(1-((x-opt)/(w/2))^2)
  res
}

# It leads to shapes like this:
curve(nbcurve(x,15,22,0.02,0.12),0,30,ylim=c(0,0.24),col='red')

# Let's try to fit this nonlinear function to our growth rate data,
# using a variety of different approaches.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### 1) Educated guessing ####

# For this approach to work, you have to really think about the data you 
# have and the deterministic equation you're fitting.

# Look at each parameter of our function in turn, and make logical guesses:

# 'opt', controls the position of the maximum of the function nbcurve()
# - So, maybe pick the temperature corresponding to the highest growth rate 
#   found in the data we're trying to fit, 15 C:
data$temps[data$grths==max(data$grths)]

# 'w', controls the width of the function nbcurve
# - So, maybe pick the range of temperatures in this growth rate data set:
diff(range(data$temps))

# 'a', corresponds to the growth rate at temperature = 0
# - So, what is the average growth rate observed at 0 C?
#   This is about 0.04
mean(data$grths[data$temps==0])

# 'b', controls how quickly growth rates increase from low temperatures up until
# the optimum temperature (or 'opt' parameter) is reached. This is a little harder
# to guess, but we can do a little math. From the function, we  have:
#     growth(x) = a*exp(b*x)*(1-((x-opt)/(w/2))^2)
# If we focus on growth rate at the optimum temperature, then x = 'opt', and our
# equation simplifies to:
#     growth(opt) = a*exp(b*opt)
# We already came up with guesses for 'a' and 'opt', and we can estimate growth(opt)
# as well - just caclulate the mean growth rate from our data set for a temperature
# of 15 C:
    mean(data$grths[data$temps==15])
# which yields y = 0.207. Doing a little algebra, we can now solve for 'b':
#   a*exp(b*opt)=y    
#   exp(b*opt)=y/a
#   b*opt = log(y/a)
#   b = log(y/a)/opt
#   b = log(0.207/0.04)/15 = 0.11     
log(0.207/0.04)/15    
        
# Or, more messily:
log(mean(data$grths[data$temps==15])/0.04)/15

### Now we have guesses for all of the parameters of nbcurve(), which we can use 
# in our mle fit

# use these guesses in an mle...
fit<-mle2(grths~dnorm(mean=nbcurve(temps,o,w,a,b),sd=s),start=list(w=27,o=15,a=0.04,b=0.11,s=10),data=data)
summary(fit)

# look at the resulting fit
cfs<-as.list(coef(fit))
plot(grths~temps,data=data)
curve(nbcurve(x,cfs$o,cfs$w,cfs$a,cfs$b),0,30,ylim=c(0,0.5),add=T,col='red')

# How well did we do?  What changed?


####__Bonus Exercise__####

# Why might a & b have changed so much, comparatively?
# Hint: Try plotting the likelihood surface for this model as a function of a and b,
#       while holding opt and w fixed at their mle estimates.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### 2) Trial and error/plotting ####

library(manipulate)

# interactive plotting
?manipulate

plot(grths~temps,data=data)

# make a little plotting function that shows the data and a potential curve
testplot<-function(o,w,a,b){
  plot(grths~temps,data=data)
  curve(nbcurve(x,o,w,a,b),0,30,ylim=c(0,0.5),add=T,col='red')
}
testplot(12,32,0.11,0.05)

# Check out the interactive plot
manipulate(testplot(o,w,a,b),w=slider(0,50,initial=10),o=slider(0,40,initial=12),a=slider(-0.5,0.5,initial=0.05),b=slider(0,0.5,initial=0.01))

# Try out different values of the parameters until the curve looks like a good
# fit to the data 'by eye'.

# Then use these guesses in an mle...
fit<-mle2(grths~dnorm(mean=nbcurve(temps,o,w,a,b),sd=s),start=list(w=29,o=15,a=0.204,b=0.01,s=10),data=data)
summary(fit)

# Look at the resulting fit
cfs<-as.list(coef(fit))
plot(grths~temps,data=data)
curve(nbcurve(x,cfs$o,cfs$w,cfs$a,cfs$b),0,30,ylim=c(0,0.5),add=T,col='red')


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### 3) Brute force ####

# Load a tool I've built called grid.mle
source("/Users/colin/Teaching/ELME/ELME 2017/labs/Lab4 - strategies for tough fits/grid.mle2 053114.R")

# We have 5 parameters to estimate:
# w, o, a, b, and s 

# We can make either a single guess, or a set of guesses for each variable individually

# Parameters for which we want to specify a set of guesses get listed in a single list
# (I usually call it grids)
grids<-list(w=seq(5,30,5),o=seq(10,25,5))

# All parameters (including those in grids) must also be specified in a second list,
# which I usually call 'start'. Parameters entered in grids are set to NA in this list,
# whereas parameters we want to specify single estimates to get numerical guesses:
start<-list(w=NA,o=NA,a=0.11,b=0.05,s=10)

# Now run grid.mle2. It will repeatedly call mle2 using all of the possible combinations
# of starting parameters that we have supplied, and return a long list of model fits.
fit<-grid.mle2(minuslogl=grths~dnorm(mean=nbcurve(temps,o,w,a,b),sd=s),grids=grids,start=start,data=data)


### The result of a grid.mle2 fit contains 3 objects,
names(fit)

# a) res.mat - a matrix of the parameter estimates from the mle fit using each 
#     combination of starting guesses
head(fit$res.mat)

# b) res.mod - a list of all of the model fits for each combo of starting guesses
fit$res.mod[1:3]

# c) res.best - the best model out of the whole set, determined by AIC
summary(fit$res.best)


# We can do everything with fit$res.best that we can do with a usual mle2 fit:
coef(fit$res.best)
pf<-profile(fit$res.best)
plot(pf)
confint(pf)

# look at the resulting fit
cfs<-as.list(coef(fit$res.best))
plot(grths~temps,data=data)
curve(nbcurve(x,cfs$o,cfs$w,cfs$a,cfs$b),0,30,ylim=c(0,0.5),add=T,col='red')



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Notes ####

# These approaches are complementary, and come with different pros and cons.
# - fitting using grid.mle2 is automated and doesn't require a lot of thinking,
#   which can be great, but also can become very computationally demanding very quickly.
# - coming up with logical guesses takes more thought, but can yield robust, 
#   quick to calculate results.
# - the manipulate() approach is visually appealing, and requires less thinking about math,
#   but some models become very difficult to visualize, and writing code to display them
#   can be time consuming. Finally, the human eye can be both remarkably good and pretty
#   terrible at making good guesses, perhaps especially when parameter covariance is high.





#### Plotting confidence bands around nonlinear MLE fits:



# use these guesses in an mle...
fit<-mle2(grths~dnorm(mean=nbcurve(temps,o,w,a,b),sd=s),start=list(w=27,o=15,a=0.04,b=0.11,s=10),data=data)
summary(fit)

# look at the resulting fit
cfs<-as.list(coef(fit))
plot(grths~temps,data=data)
curve(nbcurve(x,cfs$o,cfs$w,cfs$a,cfs$b),0,30,ylim=c(0,0.5),add=T,col='red')

# Follow the delta method (see Bolker book, pg. 255):
library(emdbook)

# confidence band for pre-drought pattern (see Bolker book, delta method, pg 255)
xs<-seq(min(data$temps),max(data$temps),0.1)
dvs<-deltavar(fun=nbcurve(xs,o,w,a,b),meanval=coef(fit),Sigma=vcov(fit))
sdapprox<-sqrt(dvs)
mlmean<-nbcurve(xs,cfs$o,cfs$w,cfs$a,cfs$b)
pred.vals<-data.frame(temps=xs,ml.mean=mlmean,ml.se=sdapprox,ml.ci=1.96*sdapprox)

library(ggplot2)

ggplot(pred.vals,aes(x=temps,y=ml.mean))+
  geom_line()+
  geom_ribbon(aes(ymin=ml.mean-ml.ci,ymax=ml.mean+ml.ci),alpha=0.2)+
  geom_point(data=data,aes(y=grths))+
  coord_cartesian(ylim=c(0,0.4))+
  theme_bw()

