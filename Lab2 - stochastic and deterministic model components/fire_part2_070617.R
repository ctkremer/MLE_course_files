#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Lab 2: Fire Part 2
# ELME MLE course
# Colin T. Kremer

# Reading:
#	Grace and Keeley 2006, (Grace2006.pdf on wiki)

# Data:
#	Keeley_rawdata_ELME.csv (on wiki)

# Topics covered:
#	- fitting deterministic models to species richness data from Grace and Keeley 2006
#	- combined various stochastic and deterministic components to model species richness
#	- saw some different ways of visualizing the results
#	- delved more deeply into model comparison techniques
#		- AICc and BIC approaches
#	- learned about estimating confidence intervals based on the shape of the likelihood surface, including:
#		- likelihood slices
#		- likelihood profiles
#		- Fisher Information
#	- looked briefly at how our estimates of the standard errors of each coefficient are used to estimate p-values signifying whether or not each coefficient is significantly different from 0 (ie, does it have an effect.)  This result is only as good as the validity with which we can apply the Fisher Information approach.

#----------------------------------------------------------------#
#----------------------------------------------------------------#

library(bbmle)	# load up our mle tools.
source("/Users/colin/Teaching/ELME/ELME 2017/labs/mle.tools_070717.R")

# set your working directory
setwd("/Users/colin/Teaching/ELME/ELME 2017/labs/Lab2 - stochastic and deterministic model components/")

# load the data
data<-read.csv("Keeley_rawdata_ELME.csv")
# data file can be found on the class wiki under lab activities

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Recall Part 1 ####

# Previously we investigated models that attempted to use various 
# probability distributions to describe variation in observed species 
# richness and fire severity:
#	 - variation in species richness was best described by a normal 
#     distribution (or possibly negative binomial distribution)
#	 - variation in fire severity was best described as normal or 
#     gamma distributed.

# At least in the case of species richness, however, the normal/negative
# biniomial distributions, despite being the 'best', weren't necessarily 
# very satisfactory in terms of how well they actually fit the data.  
# At the end of the lab, we thought about some other factors that might
# be influencing variation in species richness.  In this lab, we'll 
# investigate some of these hypotheses.

# Remember - distribution of species richness values across study
plot(density(data$rich,bw=2),xlab='Species richness',lwd=2,ylim=c(0,0.06),main="")

# So far, we've treated all of this variation as if it were coming from
# some unobserved, stochastic process.  Now we'll turn to trying to 
# explain trends in this variation in a more mechanistic way, combining
# deterministic models with our stochastic model of error/variation in 
# species richness.


# Generate a big plot with all possible pairwise scatterplot comparisons
# of the variables in data frame called 'data'. Usually this is too 
# overwhelming, but sometimes it provides a quick and dirty first glance
# at possible bi-variate relationships in the data
pairs(data)

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Graphical exploration ####

# Grace and Keeley 2006 found a relationship between species richness
# and the distance of sites from the pacific coast (a proxy for many 
# environmental gradients, including moisture and elevation)

# Let's explore this pattern visually:
plot(rich~distance,data=data)


#~~~ Digression - making a better plot ~~~#

# add axes labels
plot(rich~distance,xlab='Distance from coast',ylab='Species richness',data=data)

# add plot title
plot(rich~distance,main='Diversity response to fire',xlab='Distance from coast',ylab='Species richness',data=data)

# change the x and y range of the plot
plot(rich~distance,main='Diversity response to fire',xlab='Distance from coast',ylab='Species richness',xlim=c(0,100),ylim=c(0,90),data=data)

#~~~ End Digression ~~~#


# Capture and plot local trends in data:
lw.fit<-lowess(data$rich~data$distance)
lines(lw.fit,lwd=2,col='red')

# lowess() fits a smooth but wiggly function to the data, helpful for 
# detecting local trends. Methods aren't important for us, just provides
# a visual guide to motivate our next steps.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Adding a deterministic covariate to our model of species richness ####

# Let's return to thinking about our poisson model, even though it 
# seemed too simplistic in our earlier analyses - it's a good place to start.

# The poisson distribution has a single parameter, lambda, which 
# describes both its mean and its variance. Previously we used MLE to
# estimate the most likely value of lambda, give our data:

fit.pois<-mle2(rich~dpois(lambda=l),start=list(l=mean(data$rich)),data=data)
summary(fit.pois)
# and found that the MLE estimate of lambda = 49.23333

# This model implicitly assumes that the mean of the poisson distribution
# is insensitive to the distance of the plots from the pacific coast.
# We could write this as a linear function in the following way:

# lambda = l + l1*distance, where l1 = 0.

# What we'd like to try now is to fit a model where lambda changes with
# distance. In a linear model, lambda would either increase or decrease 
# linearly with increasing distance from the coast. Our model would be 
# identical, except now l1 can take on any value (and provides the slope
# of the linear relationship between mean species richness and distance 
# from the coast.)

# lambda = l + l1*distance, where l1 can be any value between -Inf and +Inf.

# This is the first model we'll try, and is essentially the same as 
# performing linear regression, with l being the intercept and l1 being
# the slope.


### Fit a covariate model

# Method 1: explicit negative log likelihood calculator
pois.distance.NLL<-function(l0,l1){
	lambda.values<-l0+l1*data$distance
	-sum(dpois(data$rich,lambda=lambda.values,log=T))
}

fit.pois.dist.NLL<-mle2(pois.distance.NLL,start=list(l0=30,l1=0.5),data=data)
summary(fit.pois.dist.NLL)

# Method 2: provide a statistical model formula in call to mle2

fit.pois.dist<-mle2(rich~dpois(lambda=l0+l1*distance),start=list(l0=30,l1=0.5),data=data)
summary(fit.pois.dist)



#~~~ Digression: ~~~#

# What does this error message mean/where does it come from?

#	Warning messages:
#	1: In dpois(x, lambda, log) : NaNs produced

# Thinking back to the properties of the poisson distribution, 
# one of the constraints is that the mean, lambda, must always be 
# greater than zero.

# So, if the optimizer underneath mle2 examines values for l0 
# and l1 such that the estimate of lambda is less than zero, 
# we'll get a warning message (for example, if we provide a starting
# guess for l0 set to -30):

mle2(rich~dpois(lambda=l0+l1*distance),start=list(l0=-30,l1=0.5),
     data=data,eval.only = T)

mle2(rich~dpois(lambda=l0+l1*distance),start=list(l0=30,l1=0.5),data=data)

# Note that option 'eval.only=T' forces mle2 to return only the 
# likelihood at the initial values, rather than actually optimizing.

# There are some tricks for preventing this kind of issue that 
# we'll run into later. For this example at least, the optimizer is 
# robust enough that even though it complains, it eventually settles
# on valid parameter estimates.

#~~~ End digression ~~~#



# Show the results of our fit graphically
plot(rich~distance,data=data)

# Add our regression line:
?abline()
abline(coef(fit.pois.dist))

# Because we're trying to explain the same y data (species richness), 
# we can use AIC model comparison to compare our basic poisson model
# (with constant lambda) with our shiny new model in which lambda 
# changes linearly with distance:

# AIC model comparison
AICtab(fit.pois,fit.pois.dist)

# Which model does better?

# What do these results mean, biologically?


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Developing a quadratic model ####

# Adding information about the distance of sites from the coast
# seems to greatly improve our ability to model variation in 
# species richness. So far we've examined a linear relationship 
# (about as simple as it gets).  But there's no particular 
# reason to believe that this relationship should be only linear
# - and there are tons of other functional forms we could try.

# In general, using MLE, we can check out pretty much any model of the form

# lambda = f(covariates)

# where f is some function (linear or nonlinear) of any covariate(s)
# that we are interested in, including, but not limited to, distance
# from the coast.

# Let's look at our local regression guide again:
lines(lw.fit,lwd=2,col='red')

# This suggests that there may be some important variation in the 
# relationship between species richness and distance from the coast
# not captured by our linear trend - species richness values are 
# on average lower than predicted by our linear model both close to
# the coast and far away, while at intermediate distances, species 
# richness is on average higher than we predict.


####__Exercise__####

# Try to fill in the code below.

# Let's try out a quadratic regression (second order polynomial) next:
# lambda = l0 + l1*distance + l2*distance^2


# Method 1: write NLL calculator

pois.dist2.NLL<-function(){
	
}

fit.pois.dist2<-mle2(pois.dist2.NLL,start=list(_____?______),data=data)


# Method 2: statistical model inside of mle2()

fit.pois.dist2<-mle2(____?_____,start=list(___?____),data=data)



# model comparison
AICtab(fit.pois,fit.pois.dist,fit.pois.dist2)

# visualize these results
plot(rich~distance,data=data)

# abline() won't work any more - only for linear regression, so we have to plot this 
# another way:

# the easiest way involves using the ?curve() function:
cfs.dist2<-coef(fit.pois.dist2)
curve(cfs.dist2[[1]]+cfs.dist2[[2]]*x+cfs.dist2[[3]]*x^2,0,100,col='blue',add=T)

# compare our quadratic fit to the linear fit and local trend:
abline(fit.pois.dist)
lines(lw.fit,col='red')


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Visualization of variance modelled by quadratic regression ####

# Think about these results from our original question: How do we explain variation in 
# species richness?

# Each different site located a different distance from the coast will have variable species
# richness described by a poisson distribution with a unique mean (driven by our hypothesized 
# deterministic, quadratic relationship). Collectively, the total variation in species 
# richness will be described by:
#	1) variation in mean species richness (poisson means) driven by distance from the coast 
#    (deterministic), and
#	2) remaining (stochastic) variation from a poisson distribution with the specified mean.


# Let's try to visualize what's happened in an intuitive way.

# **NOTE** This graphic is a little complex; don't get caught up on the details of
# how it's made (unless you want to). The important thing is to make the image, and 
# think about what it means.

# Remember, the model we've constructed accounts for variation in species richness using
# a poisson distribution, whose mean depends on (is conditonal on) distance from the coast.
# We can visualize how the poisson distribution, conditional on distance from the coast,
# changes over the range of our regression using the following code. This is my favorite
# way of thinking about what regressions are doing.
plot(rich~distance,data=data,ylim=c(0,100),xlab='Distance',ylab='Species richness')
cfs.dist2<-coef(fit.pois.dist2)
curve(cfs.dist2[[1]]+cfs.dist2[[2]]*x+cfs.dist2[[3]]*x^2,0,100,col='blue',add=T)
pred.x<-seq(5,85,20)
pred.y<-predict(fit.pois.dist2,newdata=list(distance=pred.x))
for(i in 1:length(pred.x)){
	rng<-30
	y2<-seq(round(pred.y[i])-rng,round(pred.y[i])+rng,1)
	xshift<-100*dpois(y2,lambda=pred.y[i])
	x2<-pred.x[i]-xshift
	lines(x2,y2,lwd=2)
	lines(c(pred.x[i],pred.x[i]),c(pred.y[i]-rng,pred.y[i]+rng),lty=3,lwd=0.75)
	points(pred.x[i],pred.y[i],pch=19)
}

# Here the blue line shows our quadratic fit (predicting species richness) and the 
# open circles show observed data points. Black points are examples of sites that
# might exist, and the black curves illustrate the stochastic variation (from a 
# poisson distribution) around the mean species richness we would predict at that 
# site using our quadratic trend.




### To compare what we've done in this script to the visual framework we were using
# in the previous lab:

# You could imagine adding up all of the densities contributed by each site's poisson
# distribution, across the range of species richnesses. This aggregated distribution
# (when appropriately normalized) would represent our new predicted distribution of
# species richness, and can be viewed in the same way we looked at predicting variation
# in species richness previously. (ie, we can look at the frequency/density of 
# different species richness values, without explicitly visually showing distance from
# the coast.)

# If you feel like some probability terminology, what I described above is basically the
# 'marginal distribution' of species richness. It describes variation in species richness
# integrated across the entire range of distances from the coast. We're collapsing down 
# a dimension, in effect.

### Let's visualize this:

# Generate predictions, add up contributions from each individual site's possion 
# distribution, and divide by the sample size (each poisson distribution integrates
# to 1, so 90 of them would integrate to 90 - and we want to produce a new 
# distribution that integrates to 1)

# *** NOTE: the following code will ONLY work if you used mle2 Method #2 above to
# create the model fit.pois.dist2. The predict() function needs to see the formula
# used to specify a statistical model in order to work. When we use Method #1 to fit
# mle2 models, the formula gets obscured within our NLL calculator.
# To see the difference, you can use the summary() function on fit.pois.dist2 after
# using Method 1 and 2. Look at the "Call" section; for Method 1, we see:
#   "minuslogl = pois.dist2.NLL"
# and for Method 2
#   "minuslogl = rich ~ dpois(lambda = l0 + l1 * distance + l2 * distance^2)"
# The predict() function isn't smart enough to know how to work with pois.dist2.NLL,
# but it can make use of the statistical formula provided for Method 2.

xs<-seq(0,100,1)
ys<-predict(fit.pois.dist2)
dense<-c()
for(i in 1:length(xs)){
	dense<-append(dense,sum(sapply(ys,function(y) dpois(xs[i],lambda=y)))/90)
}

# Plot it:
plot(1,type="n",xlim=c(0,100),ylim=c(0,0.06),xlab="Species richness",ylab="Density",main='')
for(i in 1:length(ys)){
	curve(dpois(x,lambda=ys[i])/3,0,100,lwd=0.5,add=T)
}
lines(xs,dense,col='red',lwd=2,lty=1)
# The red line shows this new, marginal distribution we predicted, describing variation in 
# species richness along with all of the individual poisson distributions (conditional 
# distributions, whose means depend on distance from coast) that combine to form the marginal
# distribution in red.


## Let's compare this distribution to those that we've found in previous models 
# without covariates:

# Recall our poisson fit to species richness data w/out considering covariates. 
# (from fire_part1.R script)
fit.pois.2<-mle2(rich~dpois(lambda=l),start=list(l=mean(data$rich)),data=data)
l.mle<-coef(fit.pois.2)[[1]]

# Now for comparison, we can visualize the fits of our original poisson model and
# the marginal distribution of our new poisson model incorporating distance:
plot(density(data$rich,bw=4),xlab='Species richness',lwd=2,ylim=c(0,0.06),main="")
curve(dpois(x,lambda=l.mle),0,100,col='red',lwd=2,add=T)
lines(xs,dense,col='red',lwd=2,lty=3)
legend(0,0.062,legend=c('data distribution','poisson distribution','distance^2 + poisson'),lty=c(1,1,3),col=c('black','red','red'),lwd=2,cex=0.5)

# This new marginal distribution appears to do a much better job of capturing the observed variation 
# in species richness, in comparison to the original poisson. It's multi-modal, like the 
# data distribution, and much more platykurtic (it has 'shoulders').




#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Compare with earlier models ####

# Run model fits from previous lab:
fit.pois<-mle2(rich~dpois(lambda=l),start=list(l=mean(data$rich)),data=data)
fit.norm<-mle2(rich~dnorm(mean=mu,sd=s),start=list(mu=mean(data$rich),s=sd(data$rich)),data=data)
fit.nbinom<-mle2(rich~dnbinom(size=s,mu=m),start=list(m=mean(data$rich),s=10),data=data)

# Save coefficients
norm.cfs<-coef(fit.norm)
nbinom.cfs<-coef(fit.nbinom)
pois.cfs<-coef(fit.pois)

# Visualize various fits to the species richness distribution
plot(density(data$rich,bw=4),xlab='Species richness',lwd=2,ylim=c(0,0.06),main="")
curve(dnorm(x,mean=norm.cfs[[1]],sd=norm.cfs[[2]]),0,100,col='blue',lwd=2,add=T)
curve(dpois(x,lambda=pois.cfs[[1]]),0,100,col='red',lwd=2,add=T)
curve(dnbinom(x,size=nbinom.cfs[[2]],mu=nbinom.cfs[[1]]),0,100,col='purple',lwd=2,add=T)
lines(xs,dense,col='red',lwd=2,lty=3)
legend(0,0.062,legend=c('data distribution','poisson distribution','poisson + distance^2','normal distribution','negative binomial distribution'),lty=c(1,1,3,1,1),col=c('black','red','red','blue','purple'),lwd=2,cex=0.5)

# Remember that the poisson and negative binomial distributions are actually discrete, 
# so really we shouldn't draw them as smooth curves:

# Visualize fits recognizing discrete nature of the distributions
plot(table(data$rich),xaxt="n",ylab='Frequency',xlab='Species richness',ylim=c(0,8),xlim=c(0,100))
axis(1,at=seq(0,100,5))
points(seq(0,100,1),90*dpois(seq(0,100,1),lambda=pois.cfs[[1]]),col='red',pch=19)
points(xs,90*dense,col='red',pch=1)
points(seq(0,100,1),90*dnbinom(seq(0,100,1),size=nbinom.cfs[[2]],mu=nbinom.cfs[[1]]),col='purple',pch=18)
curve(90*dnorm(x,mean=norm.cfs[[1]],sd=norm.cfs[[2]]),0,100,col='blue',lwd=2,add=T)
legend(0,8,legend=c('data distribution','poisson distribution','poisson+distance^2','negative binomial distribution','normal distribution'),lty=c(1,0,0,0,1),col=c('black','red','red','purple','blue'),pch=c(-1,19,1,18,-1),lwd=c(2,1,1,1,2),cex=0.5)


# At this point I want to stress again that these graphics are intended to be teaching aids - 
# they give us a rough visual sense of why one model might be better than another at 
# describing variation in our dependent variable (species richness). To be rigorous,
# we need to make use of actual quantitative comparisons of models and model fits.


# So, we've just been doing some visual comparisons of these model fits, but let's 
# get quantitative again.
AICtab(fit.norm,fit.nbinom,fit.pois,fit.pois.dist,fit.pois.dist2,weights=T)

# What do these results suggest to you?

# What could you do to advance the modeling that we've done, in light of what we've 
# just found out?



####__Exercise__####

# Consider linear and quadratic deterministic relationships between distance and species 
# richness in the context of stochastic (error) distributions that are normal, negative binomial, etc.


# What's the best, biologically reasonable model that you can come up with?





#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### ~More on model comparison and evaluation~ ####


#### Variations on AIC ####

# We've been focusing on using Akaike Information Criteria, but many other flavors of 
# 'information criteria' exist, two of which are worth pointing out at this point:

### AICc - AIC corrected for sample size
# http://en.wikipedia.org/wiki/Akaike_information_criterion#AICc
# Chapter 6, pg. 210 in Bolker.

# AICc = -2*ln(likelihood) + 2*k + (2*k*(k+1))/(n-k-1) = AIC + (2*k*(k+1))/(n-k-1) 
# where k is the number of parameters and n is the number of data points

# AICc introduces a greater penalty for the number of parameters in a model 
# (relative to the number of data points)

# the bbmle package comes with a function that lets us calculate AICc
AICc(fit.norm,nobs=length(data$rich))

### BIC - Bayesian Information Criteria
#http://en.wikipedia.org/wiki/Bayesian_information_criterion

# Derived from a bayesian approach, rather than from information theory, this 
# criteria also bases its value on the number of data points a model is being fit to:

# BIC = -2*ln(likelihood) + k*ln(n)
# where k is the number of parameters and n is the number of data points

# in this sense, BIC, like AICc, can be more conservative than AIC in the presence of 
# small amounts of data.

# calculating BIC:
BIC(fit.norm)

# We can make tables comparing models based on AIC, AICc, or BIC:

AICctab(fit.norm,fit.nbinom,fit.pois,fit.pois.dist,fit.pois.dist2,weights=T,nobs=length(data$rich))

BICtab(fit.norm,fit.nbinom,fit.pois,fit.pois.dist,fit.pois.dist2,weights=T,nobs=length(data$rich))

AICtab(fit.norm,fit.nbinom,fit.pois,fit.pois.dist,fit.pois.dist2,weights=T)

# Note that because of the formulas for these information criteria, models having 
# the same number of parameters won't demonstrate changes in their AICc/BIC/AIC values 
# with respect to each other, while the differences between models having different 
# numbers of parameters will shift.

# Ideally, we have sufficient amounts of data and are able to draw the same conclusions
# based on AIC and AICc tables and no further thought is required. When the conclusions
# of these approaches differ, we'd need to think more carefully about the conclusions we make.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Confidence intervals ####

# we'd like to be able to estimate confidence intervals around our maximum likelihood estimates of parameters.
# this helps us describe how certain or uncertain we are of our estimates, and to interpret our results.
# There are (at least) three different ways to calculate confidence intervals for maximum 
# likelihood estimates, which we'll cover in lecture and in the following portions of this lab:
#   1) Likelihood slices
#   2) Likelihood profiles
#   3) Fisher information approach

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Likelihood slices ####

# So far we've used maximum likelihood estimation to estimate the precise value of a single
# parameter (lambda). But this kind of 'point estimate' isn't terribly useful on its own, 
# as it seems extremely precise, and leaves us with no sense of just how certain or uncertain
# we are in the value of the parameter that we selected. 

# one way of capturing or quantifying our uncertainty in our point estimate for lambda 
# is by deriving confidence intervals for our estimate. For univariate models like this one 
# (having a single parameter), this process is pretty straightforward.

# Here's our poisson negative log likelihood function from earlier:
poissonNLL<-function(lambda){
	-sum(dpois(data$rich,lambda,log=T))	
}

# We can use this function to take a look at our likelihood surface with respect to different
# values of lambda:

# Negative log likelihood curve
curve(sapply(x,function(y) poissonNLL(y)),45,55,ylab='Negative log likelihood',xlab='lambda')

# our maximum likelihood estimate for lambda is:
coef(fit.pois)

# the negative log likelihood at this point estimate for lambda is:
-logLik(fit.pois)

# locate this point on our likelihood curve:
points(coef(fit.pois),-logLik(fit.pois),pch=19,col='red')


# Now let's think about what happens to the value of the negative log likelihood as 
# we increase (or decrease) lambda, moving away from the maximum likelihood estimate 
# we obtained.

# As lambda changes, the negative log likelihood increases to either side of our best
# estimate of lambda. In some sense, the 'badness of fit' increases as we move away 
# from the point where lambda=49.233.

# Turns out there's a lot of information contained in the shape of this likelihood 
# curve that can inform our efforts to describe how confident we are that lambda=49.233.

# If you recall from lecture, when we've moved far enough away from lambda=49.233 
# (our point estimate) that the negative log likelihood increases to be 1.92 units 
# higher than the negative log likelihood of our point estimate, we've reached one 
# extreme of a 95% confidence interval.

# For reference, we can add a line to our graphic showing where this critical value falls:
abline(-logLik(fit.pois)[1]+1.92,0,lty=3)

# We could calculate this critical value (based on the chi-squared distribution):
qchisq(0.95,df=1)/2
crit.value <- -logLik(fit.pois)[1]+qchisq(0.95,df=1)/2


### Now we'd like to figure out just what exactly are the lambda values having negative 
# log likelihoods equal to that of this critical value.

# 2 approaches:

# Method 1: numerical root finding
?uniroot()

lower.root<-uniroot(function(x){poissonNLL(x)-crit.value},lower=47,upper=48)
upper.root<-uniroot(function(x){poissonNLL(x)-crit.value},lower=50,upper=51)

# So our estimated confidence interval is:
c(lower.root$root,upper.root$root)


# Method 2: approximation
# (follows Bolker's EMD book, pg 193-194)
?approx()

# generate a list of negative log likelihoods
xs<-seq(46,52,0.1)
ys<-sapply(xs, poissonNLL)

# lower bound of the confidence interval
prof.lower<-ys[1:which.min(ys)]
prof.lvec<-xs[1:which.min(ys)]
lower.ci<-approx(prof.lower,prof.lvec,xout=crit.value)

# upper bound of CI
prof.upper<-ys[which.min(ys):length(ys)]
prof.lvec<-xs[which.min(ys):length(ys)]
upper.ci<-approx(prof.upper,prof.lvec,xout=crit.value)

# resulting confidence interval:
c(lower.ci$y,upper.ci$y)

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Likelihood slices in multi-parameter models ####

# We can calculate these likelihood slices for multi-parameter models as well, using our 
# negative log likelihood functions:

pois.dist.NLL<-function(l0,l1){
	-sum(dpois(data$rich,lambda=l0+l1*data$distance,log=T))
}

# Maximum likelihood point estimates for our parameters:
pois.dist<-mle2(pois.dist.NLL,start=list(l0=35.6,l1=0.884),data=data,method="Nelder-Mead")
cfs<-coef(pois.dist)
cfs

# Critical value for 95% confidence intervals/Likelihood Ratio test
critical.value <- -logLik(pois.dist)[1]+qchisq(0.95,df=1)/2

# Calculate our negative log likelihood slices by varying each parameter in turn, while 
# holding the other 2 fixed.

# Visualize each of these likelihood slices, one for each parameter.
par(mfrow=c(1,2))
curve(sapply(x, function(y) pois.dist.NLL(y,cfs[[2]])),39,43,ylab="NLL",xlab="l0",main="Likelihood slice")
abline(critical.value,0,lty=3)
curve(sapply(x, function(y) pois.dist.NLL(cfs[[1]],y)),0.2,0.3,ylab="NLL",xlab="l1",main="Likelihood slice")
abline(critical.value,0,lty=3)
par(mfrow=c(1,1))

####__Exercise__####

# Calculate the corresponding estimates of confidence intervals for 
# l0 and l1, using either of the above approaches.



#----------------------------------------------------------------#
#----------------------------------------------------------------#


#### Likelihood profiles ####

# There's a draw-back to this approach to finding a confidence interval.
# When we work with more than 1 parameter, parameter estimates can co-vary, and 
# likelihood surfaces get more complicated.  Ultimately, likelihood slices are a biased 
# way of estimating confidence intervals in these situations. Slices assume that 
# varying parameters in orthogonal combinations traces the curvature of the likelihood 
# (negative log likelihood) surface taking on the maximum (minimum) values.

# Let's try to visualize this surface.
pnts<-expand.grid(x=seq(38,45,0.01),y=seq(0.18,0.34,0.005))
NLL.values<-rep(NA,dim(pnts)[[1]])
for(i in 1:dim(pnts)[[1]]){
	NLL.values[i]<-pois.dist.NLL(pnts[i,1],pnts[i,2])
}
pnts$NLL<-NLL.values

# ***NOTE*** 
# We're busting out some fancy graphics; again, don't get lost in the code
# that accomplishes this (unless you want to). Mostly I want you to look
# at the resulting visual and think about what it means.

# load an extra library with tools for graphing
library(ggplot2)

p1<-ggplot(mapping=aes(x,y,z=NLL),data=pnts)+
  geom_tile(aes(fill=NLL),data=pnts)+
  stat_contour(bins=20,size=0.25,data=pnts)+
  stat_contour(breaks=c(critical.value),colour='red',size=0.6,data=pnts)+
  geom_segment(aes(x=cfs[[1]],y=min(pnts$y),xend=cfs[[1]],yend=max(pnts$y)),colour='white',linetype=2,size=1)+
  geom_segment(aes(x=min(pnts$x),y=cfs[[2]],xend=max(pnts$x),yend=cfs[[2]]),colour='white',linetype=2,size=1)+
  geom_point(aes(x=cfs[[1]],y=cfs[[2]]),colour='red')+
  scale_x_continuous(name="l0",expand = c(0,0))+
  scale_y_continuous(name="l1",expand = c(0,0))
p1  

# This graph plots the negative log likelihood surface, and shows the location of our 
# MLE estimates for l0 and l1 (red dot). The red ellipse shows the contour of the NLL 
# surface corresponding to our critical value. Can you see how to read off 
# likelihood-slice confidence intervals from this plot?



### Concept of likelihood profiles:

# As discussed in lecture, we can account for this potential covariance between 
# parameters by calculating likelihood profiles instead of likelihood slices.

# Instead of varying 1 parameter and holding all of the others constant at their 
# MLE estimate, we switch to varying 1 parameter and re-fitting the model (ie, 
# all of the other parameters can adjust their value in order to minimize the 
# corresponding negative log likelihood). This finds the 'most likely' combination
# of other parameters, given the fixed value of our focal parameter, and the data.

# Here's an example of calculating a profile by hand:
# 	- let's run it for l0 for starters.

l0.values<-seq(38,45,0.2)	# select a range of l0 values to investigate

i<-1
l0.prof.NLL<-rep(NA,length(l0.values))	# vector to save our negative log likelihood values in
l1.save<-rep(NA,length(l0.values))		# for saving l1 estimates
for(i in 1:length(l0.values)){	# loop through each fixed value of parameter l0
	# re-fit the model using mle2, specifying that parameter l0 should be fixed at a specific value (ie, fixed=list(l0=l0.values[i]) in this case)
	tmp.fit<-mle2(rich~dpois(lambda=l0+l1*distance),start=list(l0=35.6,l1=0.884),fixed=list(l0=l0.values[i]),data=data,method="Nelder-Mead")

	# save our findings
	l0.prof.NLL[i]<- -logLik(tmp.fit)[1]
	l1.save[i]<-coef(tmp.fit)[[2]]
}
l0.prof.NLL	# here are all of the corresponding negative log likelihood values, given the adjustments made to the value of all parameters except the fixed value of l0

# Visualize the differences:
par(mfrow=c(1,1))	
curve(sapply(x, function(y) pois.dist.NLL(y,cfs[[2]])),38,45,ylab="NLL",xlab="l0",main="l0")	
abline(critical.value,0,lty=3)
lines(l0.prof.NLL~l0.values,col='red')
points(cfs[[1]],-logLik(pois.dist),pch=3)

# Get confidence intervals from the slice:
slice.lower.root<-uniroot(function(x){pois.dist.NLL(x,cfs[[2]])-critical.value},lower=35,upper=cfs[[1]])$root
slice.upper.root<-uniroot(function(x){pois.dist.NLL(x,cfs[[2]])-critical.value},lower=cfs[[1]],upper=51)$root
c(slice.lower.root,slice.upper.root)

# Compare with confidence intervals from the profile:

# Set up a profile function
prof.l0.func<-function(x){
	tmp<-mle2(rich~dpois(lambda=l0+l1*distance),start=list(l0=35.6,l1=0.884),fixed=list(l0=x),data=data,method="BFGS")
	res<- -logLik(tmp)
	return(res[[1]])
}
prof.l0.func(40)

# Calculate confidence intervals using this profile function
prof.lower.root<-uniroot(function(x){prof.l0.func(x)-critical.value},lower=35,upper=cfs[[1]])$root
prof.upper.root<-uniroot(function(x){prof.l0.func(x)-critical.value},lower=cfs[[1]],upper=51)$root
c(prof.lower.root,prof.upper.root)

# Save the l1 estimates at these locations (for graphing later)
tmp.lw<-mle2(rich~dpois(lambda=l0+l1*distance),start=list(l0=35.6,l1=0.884),fixed=list(l0=prof.lower.root),data=data,method="BFGS")
tmp.up<-mle2(rich~dpois(lambda=l0+l1*distance),start=list(l0=35.6,l1=0.884),fixed=list(l0=prof.upper.root),data=data,method="BFGS")
l1.vals.at.ci<-c(coef(tmp.lw)[[2]],coef(tmp.up)[[2]])

# Show location of these confidence interval estimates on plot to verify them:
curve(sapply(x, function(y) pois.dist.NLL(y,cfs[[2]])),38,45,ylab="NLL",xlab="l0",main="l0")	
abline(critical.value,0,lty=3)
lines(l0.prof.NLL~l0.values,col='red')
points(cfs[[1]],-logLik(pois.dist),pch=3)
points(c(critical.value,critical.value)~c(prof.lower.root,prof.upper.root),col='red',pch=8) 
points(c(critical.value,critical.value)~c(slice.lower.root,slice.upper.root),pch=2)


### We can visualize this as the 2D likelihood surface as well, to get a 
# better understanding of what's happening:

# Set up data for plotting
tmp.df<-data.frame(l0.values,l1.save)
tmp.l0.slice.ci<-data.frame(x=c(slice.lower.root,slice.upper.root),y=c(cfs[[2]],cfs[[2]]))
tmp.l0.prof.ci<-data.frame(x=c(prof.lower.root,prof.upper.root),y=l1.vals.at.ci)

p2<-ggplot(mapping=aes(x,y,z=NLL),data=pnts)+
  geom_tile(aes(fill=NLL),data=pnts)+stat_contour(bins=20,size=0.25,data=pnts)+
	stat_contour(breaks=c(critical.value),colour='red',size=0.6,data=pnts)+
	geom_segment(aes(x=cfs[[1]],y=min(pnts$y),xend=cfs[[1]],yend=max(pnts$y)),
	             colour='white',linetype=2)+
	geom_segment(aes(x=min(pnts$x),y=cfs[[2]],xend=max(pnts$x),yend=cfs[[2]]),
	             colour='white',linetype=2)+
	scale_x_continuous(name="l0",expand=c(0,0))+
	scale_y_continuous(name="l1",expand=c(0,0))+
	geom_line(data=tmp.df,aes(x=l0.values,y=l1.save,z=NULL),colour='white')+
	geom_point(aes(x=cfs[[1]],y=cfs[[2]]),colour='red')+
	geom_point(data=tmp.l0.slice.ci,aes(x,y,z=NULL),colour='red',shape=8,size=3)+
	geom_point(data=tmp.l0.prof.ci,aes(x,y,z=NULL),colour="white",shape=2,size=3)
p2


#----------------------------------------------------------------#
#----------------------------------------------------------------#


####__Bonus Exercise__####

# Try repeating this exercise to find estimates of confidence
# intervals for the second parameter, l1, using likelihood profiling



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Likelihood profiling: the easy way... ####

# Whew - that was a lot of work!

# Fortunately, the mle2 package comes to our rescue again, and can save us 
# a great deal of code. It has built in functions that let us calculate the 
# likelihood profiles of any given mle2 model fit:
pf1<-profile(pois.dist)
pf1

# and given the profile, will calculate confidence intervals for us directly:
confint(pf1)

# we can even plot the profiles:
plot(pf1)
# here the profiles look v-shaped instead of quadratic, because the default plot shows 
# the square-root transformed likelihood profiles 

# We could plot our profiles this way, too:
plot(NULL,xlim=c(38,45),ylim=c(0,4),xlab='l0',ylab='sqrt(NLL - min(NLL))',main='likelihood profile for l0')
lines(sqrt(l0.prof.NLL-min(l0.prof.NLL))~l0.values,col='red')
abline(qchisq(0.95,df=1)/2,0,lty=3)




#~~~ Digression ~~~#

# One of the assumptions of the LRT as applied to confidence intervals requires 
# that the shape of the likelihood surface (in the direction of the profile) 
# be quadratic. When plotting the square-root of the likelihood profile, 
# quadratic shapes become linear - so if these plots look v-shaped, life is good.

#~~~ End digression ~~~#


# profile() will also work on more complicated models where it becomes very 
# difficult to visualize/plot what's happening. (although for very large models
# it often runs into trouble, as likelihood surfaces typically get bumpy and 
# badly behaved)

# Take a look at the profiles for a larger model:
fit.pois.dist2<-mle2(rich~dpois(lambda=l0+l1*distance+l2*distance*distance),start=list(l0=30,l1=0.5,l2=0),data=data,method="Nelder-Mead")

pf2<-profile(fit.pois.dist2)
confint(pf2)
plot(pf2)



#----------------------------------------------------------------#
#----------------------------------------------------------------#


#### Fisher Information ####

# As discussed in lecture, it is also possible to approximate confidence intervals that are 
# asymptotically valid by using Fisher information

# This is a function I wrote that uses the Hessian matrix calculated as part of our mle2 
# fits to estimate confidence intervals via the Fisher information approach. For
# now, just run it and look at the results; I describe how it works later on.

# function to generate confidence intervals based on Fisher Information criteria
confint.FI<-function(model){
	cfs<-coef(model)
	ses<-sqrt(diag(vcov(model)))	# standard errors
	lw<-cfs-1.96*ses
	up<-cfs+1.96*ses
	res<-cbind(lw,up)
	dimnames(res)<-list(names(cfs),c("2.5 %","97.5 %"))	
	res
}


# We can apply it directly to an mle2 model object:
confint.FI(fit.pois.dist2)

# Compare these values to our earlier estimates:
confint(pf2)

# Do the Fisher CI's seem more or less conservative than those obtained by profiling?



# Here's a step by step explanation of how the function confint.FI works 
# (and recall lecture discussion):

# save our estimated coefficient values
cfs<-coef(fit.pois.dist2)

# Now we need to get the variance-covariance matrix describing our estimates of the 
# variance or covariance of our model parameters, given the likelihood surface.

# Buried inside of each of our mle2 model objects is a hessian matrix, calculated 
# during mle2's optimization (and sometimes even used in the optimization process itself):
fit.pois.dist2@details$hessian

# we can determine the variance-covariance matrix describing our estimates of the variance 
# or covariance of our model parameters, given the likelihood surface. R does this for us 
# automatically if we ask for the variance-covariance using vcov() function (based on that Hessian matrix):
vcov(fit.pois.dist2)

# we can focus on the diagonal elements (Describing the variance we're attributing to our 
# estimates of each parameter, ultimately based on the shape of the likelihood surface):
diag(vcov(fit.pois.dist2))

# our estimates for the standard error/deviation of each of these coefficients/parameters 
# in our data is obtained by taking the square root of these variances:
ses<-sqrt(diag(vcov(fit.pois.dist2)))
ses

# if we can assume that the distribution of our MLE estimates for each parameter is 
# normally distributed (common assumption), then we can estimate the lower (upper) 
# CI values by subtracting (adding) 1.96*standard error to our maximum likelihood 
# estimate for each parameter:
lw<-cfs-1.96*ses
up<-cfs+1.96*ses
cbind(lw,up)


# 95% of the probability density of the normal distribution falls within
#     1.96 * standard deviation of the distribution

### Or, instead of doing this by hand, you can either:
# - copy and paste the confint.FI function into the R script where you're going to be 
#   doing an analysis and run confint.FI(mle2modelresult)
# - load mle.tools_060815.R before running your analyses, which sets up a bunch of useful 
#   functions such as confint.FI

confint.FI(fit.pois.dist2)


#~~~ Digression ~~~#

# Remember how we're always looking at the summaries of our mle2 fits, and we've been 
# ignoring the 'Pr(z)' or p-values attached to the individual coefficients?
summary(fit.pois.dist2)

# Now we can talk about where they come from:

# Step 1: get a z-score for each parameter estimate, given the standard normal 
# distribution's 'null hypothesis' that each parameter should be equal to 0. 
# This zscore will then help us quantify just how extremely different our coefficient
# estimate is from zero (higher or lower).


# Note:  the "z value" reported in the mle2() summaries is calculated by:
# 1) taking our coefficient estimate, and subtracting 0 from it
# 2) dividing by our estimate of the standard error of the coefficient

# eg,
coef(fit.pois.dist2)[[1]]
coef(fit.pois.dist2)[[1]]/ses[1]
zscore<-coef(fit.pois.dist2)[[1]]/ses[1]


# Step 2: from this z-score or z value, we can estimate a p-value, describing 
# the 'probability' of observing a coefficient as extreme (different from 0)
# as our estimated coefficient.

# (default is to use a 2-tailed test)

# this value is the Pr(z) or p-value estimate associated with our l0 
# coefficient in summary(fit.pois.dist2)
2*pnorm(zscore,mean=0,sd=1,lower.tail=FALSE)


# so for our l2 parameter:
2*pnorm(coef(fit.pois.dist2)[[3]]/ses[[3]],mean=0,sd=1,lower.tail=TRUE)

# compare this value to the Pr(z) value for summary(fit.pois.dist2)
summary(fit.pois.dist2)

#~~~ End digression ~~~#




#----------------------------------------------------------------#
#----------------------------------------------------------------#

####__Bonus Exercise__####

# In the previous lab you may have looked at variation in fire severity between sites.

# Building on what you found, try using covariates (like the age of a site) to help 
# explain variation in fire severity, by adding a deterministic component to your models.

# Fit models using MLE, and then use model comparison, biological sense, and the other new 
# techniques that we've covered (likelihood slices/profiles, confidence intervals) 
# to compare and evaluate these models. What can you discover?


#----------------------------------------------------------------#
#----------------------------------------------------------------#


#### Recap ####

# In this lab script, we've covered the following:
#
#	- fitting deterministic models to species richness data from Grace and Keeley 2006
#	- combined various stochastic and deterministic components to model species richness
#	- saw some different ways of visualizing the results
#	- delved more deeply into model comparison techniques
#		- AICc and BIC approaches
#	- learned about estimating confidence intervals based on the shape of the likelihood surface, including:
#		- likelihood slices
#		- likelihood profiles
#		- Fisher Information
#	- looked briefly at how our estimates of the standard errors of each coefficient are used to estimate p-values signifying whether or not each coefficient is significantly different from 0 (ie, does it have an effect.)  This result is only as good as the validity with which we can apply the Fisher Information approach.

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### "Answers" ####


# Exercise: fitting a quadratic model with a poisson distribution:

# Method 1:
pois.dist2.NLL<-function(l0,l1,l2){
  lambda.values<-l0+l1*data$distance+l2*data$distance^2
  -sum(dpois(data$rich,lambda=lambda.valuess,log=T))
}

fit.pois.dist2<-mle2(pois.dist2.NLL,start=list(l0=30,l1=0.5,l2=0),data=data)
summary(fit.pois.dist2)

# Method 2:
fit.pois.dist2<-mle2(rich~dpois(lambda=l0+l1*distance+l2*distance^2),start=list(l0=30,l1=0.5,l2=0),data=data)
summary(fit.pois.dist2)


####

# Try a variety of mle models with different distributions and deterministic components

fit.norm.dist<-mle2(rich~dnorm(mean=m0+m1*distance,sd=s),start=list(m0=30,m1=0.5,s=5),data=data)

# Note: run into the NaNs error again here - can introduce trick of wrapping sd estimate in exp()
fit.norm.dist<-mle2(rich~dnorm(mean=m0+m1*distance,sd=exp(s)),start=list(m0=30,m1=0.5,s=1),data=data)

fit.norm.dist2<-mle2(rich~dnorm(mean=m0+m1*distance+m2*distance*distance,sd=s),start=list(m0=30,m1=0.5,m2=0,s=5),data=data)

fit.nbinom.dist<-mle2(rich~dnbinom(size=s,mu=m0+m1*distance),start=list(m0=30,m1=0.5,s=5),data=data)
fit.nbinom.dist2<-mle2(rich~dnbinom(size=s,mu=m0+m1*distance+m2*distance*distance),start=list(m0=30,m1=0.5,m2=0,s=5),data=data)


AICtab(fit.norm,fit.norm.dist,fit.norm.dist2,fit.nbinom,fit.nbinom.dist,fit.nbinom.dist2,fit.pois,fit.pois.dist,fit.pois.dist2)

#                  dAIC  df
# fit.norm.dist2     0.0 4 
# fit.nbinom.dist2   2.9 4 
# fit.norm.dist     14.2 3 
# fit.nbinom.dist   14.9 3 
# fit.norm          33.6 2 
# fit.nbinom        35.6 2 
# fit.pois.dist2    87.0 3 
# fit.pois.dist    136.7 2 
# fit.pois         227.9 1 

# So, normal distribution with quadratic relationship to distance from the coast 
# is the best model by AIC comparison, closely followed by negative binomial with
# quadratic relationship

## Visualize this

# now, plot again, but normalize this result...
xs<-seq(0,100,1)
ys<-predict(fit.norm.dist2)
#plot(rich~distance,data=data)
#points(data$distance,ys,col='red')
dense<-c()
for(i in 1:length(xs)){
	dense<-append(dense,sum(sapply(ys,function(y) dnorm(xs[i],mean=y,sd=coef(fit.norm.dist2)[[4]])))/90)
}

# imagine adding up all of the densities contributed by each poisson distribution, across the range of species richnesses (and dividing by a normalizing constant, in this case the number of observations)
# that would give us the red line in the following figure.
plot(1,type="n",xlim=c(0,100),ylim=c(0,0.06),xlab="Species richness",ylab="Density",main='')
for(i in 1:length(ys)){
	curve(dnorm(x,mean=ys[i],sd=coef(fit.norm.dist2)[[4]])/3,0,100,lwd=0.5,add=T)
}
lines(xs,dense,col='red',lwd=2,lty=1)


# visualize the poisson and normal fits to our data
plot(density(data$rich,bw=4),xlab='Species richness',lwd=2,ylim=c(0,0.06),main="")
curve(dnorm(x,mean=coef(fit.norm)[[1]],sd=coef(fit.norm)[[2]]),0,100,col='blue',lwd=2,add=T)
lines(xs,dense,col='blue',lwd=2,lty=3)
legend(40,0.05,legend=c('data distribution','normal distribution','distance^2 + normal'),lty=c(1,1,3),col=c('black','blue','blue'),lwd=2)


