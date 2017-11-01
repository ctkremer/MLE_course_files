#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Lab 3: Categorical data analysis
# ELME MLE course
# Colin T. Kremer

# Reading:
#	Huberty, Gross, and Miller 1998, (Huberty1998.pdf on course site)

# Data:
#	LTER_oldfield_fertilization.csv

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Lab activity: Using categorical covariates as the deterministic portion of our models ####

# Previously we thought about species richness and distance from the pacific coast 
# (correlated to broad scale patterns in moisture/elevation and other environmental gradients)
# Distance is a continuous covariate; now let's explore ways of modeling categorical 
# covariates using MLE methods.

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Data description ####

# Old field fertilization/disturbance experiment, KBS LTER

# Former agricultural fields were abandoned from production here at KBS back in 1989.  
# Since that time, data has been collected at least yearly, quantifying the biomass of 
# all plant species that occur in these plots. The goal of the study was to investigate 
# succession in abandoned agricultural fields. Additionally, two treatments were applied 
# factorially to portions of the plots: 
# 1) Some plots received additional nitrogen (a fertilization treatment)
# 2) Some plots were tilled (to simulate disturbance)

# Together, these treatments allow KBS researchers to study how fertilization and disturbance
# interact to affect the successional trajectory of these old-field plots.

# for more on this, see papers published on this study (posted on the course site):
#	- Huberty et al. 1998.pdf
#	- Grman et al. 2010.pdf

# Today's goal:
# - investigate the effects of fertilization and disturbance on above ground plant biomass 
#   in this succession experiment.

#----------------------------------------------------------------#
#----------------------------------------------------------------#

library(bbmle)
source("/Users/colin/Teaching/ELME/ELME 2017/labs/mle.tools_070717.R")

# set working directory
setwd("/Users/colin/Teaching/ELME/ELME 2017/labs/Lab3 - categorical covars and LRTs/")

# load data
c1<-read.csv("LTER_oldfield_fert_cleaned.csv")

# Here's what the data look like:
str(c1)
head(c1)

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Data exploration ####

# Check out the treatments
boxplot(log.biomass~fertilized,data=c1)		# effect of fertilization
boxplot(log.biomass~disturbed,data=c1)		# effect of disturbance
boxplot(log.biomass~fertilized*disturbed,data=c1)	# both

# Are there temporal patterns?
boxplot(log.biomass~year,data=c1)
plot(log.biomass~year,data=c1)

# What issues might this present for analysis?


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Test for an effect of treatment using MLE ####

# As usual, before we incorporate any covariates, it's a good idea to start with a 'null' 
# model lacking any covariates, and treating the observed data as if it arose from a 
# stochastic process/through observation error.
fit.null<-mle2(log.biomass~dnorm(mean=b0,sd=s),start=list(b0=4,s=1),data=c1)
summary(fit.null)

# Or, set up slightly differently to avoid the error messages and make the optimization 
# more stable:
fit.null<-mle2(log.biomass~dnorm(mean=b0,sd=exp(s)),start=list(b0=4,s=1),data=c1)
summary(fit.null)

# We just have to make sure to interpret the parameter estimated for s appropriately. 
# It is no longer directly the MLE estimate for the standard deviation of the normal 
# distribution describing the stochastic error in our model (aka, residual standard 
# error) - but exp(s) is.


### Now investigate the effect of the fertilization treatment

# Fertilization is a discrete, categorical variable - either a plot was fertilized 
# or it wasn't. 

# Here's a visualization of variation in log.biomass broken down by fertilization 
# treatment this is the variation that we're trying to understand and model
boxplot(log.biomass~fertilized,data=c1)


# To work with this kind of data, it is helpful to generate a column (or, in the
# case of multiple categorical covariates, a matrix), containing 0's and 1's to 
# indicate whether or not a treatment has been applied. For example focusing on 
# fertilization, we could create a column containing a fertilization code, with 1 
# indicating fertilization and 0 indicating no fertilization:
c1$fert.code<-ifelse(c1$fertilized=="fert",1,0)

# As described in lecture, we can then use this coded column to set up calculations
# of our treatment effects.

# Here's a negative log likelihood calculator for this situation, with normal errors, and 
# incorporating a treatment effect for fertilization
norm.fert.NLL<-function(b0,b1,s){
	-sum(dnorm(c1$log.biomass,mean=b0+b1*c1$fert.code,sd=s,log=T))
}
# note that because of how we set up c1$fert.code, the mean of this normal 
# distribution will equal:
#	    b0+b1*1		when c1$fert.code equals 1 (ie, fertilized plot)
#	    b0+b1*0, or just b0		when c1$fert.code == 0 (unfertilized plot)

# So we can interpret our coefficients b0 and b1 in the following way:
#	  - b0 tells us the estimated mean log(biomass) in the unfertilized treatment.
#	  - b1 tells us the amount that the fertilization treatment increased (or decreased)
#     the mean log(biomass) compared to the unfertlized plots.
#	  - b0+b1 tells us the estimated mean log(biomass) in the fertilized plots.


#### Now let's run the actual MLE fits and see what we get:

# Method 1: explicit negative log likelihood calculator
fit.fert<-mle2(norm.fert.NLL,start=list(b0=4,b1=0,s=1))
summary(fit.fert)

# So our estimate for b0 is 3.802713 interpreted as the mean log(biomass) 
# for the unfertilized plots. For comparison, the value calculated directly
# from the data is 3.802812
mean(c1$log.biomass[c1$fert.code==0])

# Our estimate for b1 is 0.633651, which we interpret as the effect of 
# adding fertilizer on the mean log(biomass)
boxplot(c1$log.biomass~c1$fert.code,ylab='log(biomass)',xlab='fert.code')

# Do these results make sense?

# Model comparison
AICtab(fit.fert,fit.null)

# Does the fertilization treatment have an effect?


# Method 2: statistical formula provided to mle2
# (results should be identical...)
fit.fert<-mle2(log.biomass~dnorm(mean=b0+b1*fert.code,sd=s),start=list(b0=4,b1=0,s=1),data=c1)
summary(fit.fert)

# Compare these results to a standard, least-squares one-way ANOVA:
fit.lm<-lm(log.biomass~fert.code,data=c1)
summary(fit.lm)


#### How well do the coefficient estimates agree with each other, comparing 
# either of the maximum likelihood fits to the anova results?

# Note that parameter 's' from the mle model estimates the standard deviation
# of the normal distribution we're using to model variation in log.biomass 
# given the fertilizer treatment. It's comparable to the 'Residual standard
# error' value in summary of the lm() fit.  (the standard deviation of the 
# residuals from the lm fit)

# What are the differences?
#	- discussion of bias and mle approach
#	- recall the results of our analytical calculation of MLE estimators for 
#     the normal distribution



####__Exercise__####

# Note that we could also have written the model above using a lognormal 
# probability distribution and mle2, rather than log-transforming the 
# biomass data and using a normal distribution.  Can you set up such a model, 
# and verify that the results of the two are comparable?

# Hint: when setting up your model, be careful to pay attention to the 
# constraints of the parameters of the lognormal distribution

?dlnorm()

fit.fert.lnorm<-mle2(______?_______)

# what are the parameter values? how do they compare to those from our earlier fits?

# Question: Why is this a bad idea?
AICtab(fit.fert,fit.fert.lnorm)




### Model comparison

# Earlier, we used AIC to compare models
AICtab(fit.fert,fit.null)

# These comparisons can give us a pretty good idea of how models compare to 
# each other, including comparisons between models with covariates and 
# 'null models'.  However, they don't give us 'p-values', which can *sometimes* 
# be useful or necessary to calculate to satisfy reviewers or others preferring
# frequentist approaches.

# Fortunately, there are approaches that enable us to obtain approximate p-values 
# in certain situations, which can facilitate comparisons with traditonal p-value 
# based methods.

# Remember how we obtained confidence intervals previously?

# Via likelihood slices or profiles
fert.prof<-profile(fit.fert)
confint(fert.prof)
plot(fert.prof)

# Next we'll look at the likelihood ratio test as a way of computing p-values; 
# it shares some similarities to these approaches for finding confidence intervals.



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Likelihood Ratio Test ####


# For more on the derivation of the likelihood ratio test, see Bolker's EMD book, 
# bottom of page 191.

# In addition to being useful for estimating confidence intervals, the likelihood
# ratio test can be used to determine p-values when comparing 'nested models'.

# When we used the likelihood ratio test to derive confidence intervals, we compared
# the likelihood of our maximum likelihood best fit model to the likelihoods arising
# from fixing the value of our parameter of interest.  We searched for the fixed value
# corresponding to an increase in the negative log likelihood of 1.92 units to find 
# a 95% confidence interval.

# Now we'll use the likelihood ratio test in a similar way.  If we compare the
# following models for the mean of a normal distribution:

# 1)  mu = b0
# 2)  mu = b0 + b1*fert.code

# It's pretty easy to see that when b1 = 0, 1) and 2) are equivalent. The usual
# term for this property is that model 1) is 'nested' in model 2). 

# We can compare the model where mu = b0 + b1*fert.code and {b0,b1} are fit by 
# mle to the model where mu = b0 and {b0} is fit by mle (and b1 is implicitly 
# equal to 0) using the likelihood ratio test, just the same way that we used 
# the likelihood ratio test for confidence intervals. 

# In this particular case, the difference in the degrees of freedom between the
# two models is 1 (b1 is either fixed at 0 or allowed to vary, giving model 2
# one extra degree of freedom). We need to know the difference in the degrees of
# freedom between the models to calculate the appropriate value from the chi-squared 
# distribution. 

# In general, you can compare models that differ in more than one degree of 
# freedom, as long as the simpler model (the model with fewer degrees of freedom) 
# remains nested in the more complex model (having more degrees of freedom).  
# For example, from our previous lab exercise:

# Species richness 
#	lambda = l0
#	lambda = l0 + l1*distance + l2*distance*distance

# The first model is nested in the second model, because when l1 and l2 both equal
# zero, the second model becomes equivalent to the first. 

# Would lambda = l0 + l1*distance be nested in the second model given above?
 

# Returning to the data at hand, we can compare our model without an effect based
# on fertilization treatment (mu = b0) with our model incorporating a fertilization
# treatment effect (mu = b0 + b1*fert.code), because these models are nested.

# Using our notation from above, fit.null is the equivalent of model 1) and 
# fit.fert is equivalent to model 2)

# Save the log likelihoods of these nested models
m1.LL<-logLik(fit.null)
m2.LL<-logLik(fit.fert)

# and the degrees of freedom for each of the models
m1.df<-attr(m1.LL,"df")
m2.df<-attr(m2.LL,"df")

# calculate the likelihood ratio and difference in degrees of freedom between
# the more complex and the simpler model.
lik.ratio <- -2*m1.LL[1]+2*m2.LL[1]
df.diff<-m2.df-m1.df

# p-value estimate from likelihood ratio test.
pchisq(lik.ratio,df=df.diff,lower.tail=F)

# compare with the values from the standard frequentist/one-way ANOVA result:
summary(fit.lm)

# Maybe not terribly informative here - large amounts of data tends to make
# p-values smaller and smaller. In this case, they're so small as to be difficult to 
# compare accurately.

# This is one reason why it can be much more important to pay attention to the
# effect sizes. These results suggest that the addition of fertilizer increases
# community above ground, dry biomass by exp(0.6335) = 1.88 grams per m^2.

# Expert trick:
# you can use the function anova() to perform a likelihood ratio test on nested
# models without writing as much code - just be careful that you only provide it with properly
# nested models, as I think it isn't smart enough to complain if you send it the wrong models.
anova(fit.null,fit.fert)

# Now we'll return to looking at different kinds of deterministic models with 
# categorical covariates.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### A two-way ANOVA using MLE ####

# Now let's try to understand variation in log.biomass as a function of both
# treatments, fertilization and disturbance
boxplot(log.biomass~fertilized*disturbed,data=c1)

# We can combine information from our fert.code column and our dist.code column:

# Negative log likelihood calculator - normal errors, treatment effect of fertilization
norm.fert.dist.NLL<-function(b0,b1,b2,s){
	-sum(dnorm(c1$log.biomass,mean=b0+b1*c1$fert.code+b2*c1$dist.code,sd=s,log=T))
}

# Method 1: explicit negative log likelihood calculator
fit.fert.dist.NLL<-mle2(norm.fert.dist.NLL,start=list(b0=4,b1=0,b2=0,s=1))
summary(fit.fert.dist.NLL)

# Method 2: statistical formula provided to mle2
fit.fert.dist<-mle2(log.biomass~dnorm(mean=b0+b1*fert.code+b2*dist.code,sd=s),start=list(b0=4,b1=0,b2=0,s=1),data=c1)
summary(fit.fert.dist)

# Comparable to an ANOVA/linear model with two main effects:
fit.lm2<-lm(log.biomass~fert.code+dist.code,data=c1)
summary(fit.lm2)


### With interaction effect:

# We can also investigate a potential interaction between fertilization 
# and disturbance treatments
fit.fertxdist<-mle2(log.biomass~dnorm(mean=b0+b1*fert.code+b2*dist.code+b3*fert.code*dist.code,sd=s),start=list(b0=4,b1=0,b2=0,b3=0,s=1),data=c1)
summary(fit.fertxdist)

# Comparable to an ANOVA/linear model with two main effects and interaction:
fit.lm3<-lm(log.biomass~fert.code*dist.code,data=c1)
summary(fit.lm3)


# Expert trick:  Now that you've done a bit of work putting together your
# own coding variables and running mle fits, here's a short cut that gets 
# mle2 to do some of that work for you. It saves time, especially for models
# with many factors, but as models get more complex it doesn't always work - 
# so it's still good to know how to set these things up 'by hand'.

fit.example<-mle2(log.biomass~dnorm(mean=b0,sd=s),start=list(b0=4,s=1),
                  parameters=list(b0~fertilized*disturbed),data=c1)
summary(fit.example)	
# This generates the same results, although coefficients/parameters will 
# look different, because of a difference in how the 'contrasts' are being
# set up.


### All of these models can still be compared with AIC comparison
AICtab(fit.null,fit.fert,fit.fert.dist,fit.fertxdist)

# What do these results tell you?


### And we can continue to use Likelihood Ratio Tests to do significance testing

# This includes testing the overall significance of a model containing main effects of
# fertilizer and disturbance treatments, as well as an interaction, versus a null model

# p-value estimated via the likelihood ratio test
m1.LL<-logLik(fit.null)
m2.LL<-logLik(fit.fertxdist)
m1.df<-attr(m1.LL,"df")
m2.df<-attr(m2.LL,"df")

lik.ratio <- -2*m1.LL[1]+2*m2.LL[1]
df.diff<-m2.df-m1.df

pchisq(lik.ratio,df=df.diff,lower.tail=F)

# or 
anova(fit.null,fit.fertxdist)

# Note that if we look at the output for the standard linear model/ANOVA:
summary(fit.lm3)

# We see several p-values. There's one for the whole model (4.328e-10)
# and one for each of the individual effects in the model. The 
# likelihood ratio test we did above generates a p-value analogous to
# the whole model p-value, because it compares the full model with the 
# null model. However, we can re-create the p-value estimates for the 
# individual effects via conducting likelihood ratio tests on appropriate
# pairs of nested models.

# For example, we can look at the significance of just the interaction effect:

# p-value estimated via the likelihood ratio test, comparing model with
# an interaction to the model with only main effects (which is nested 
# in the full interaction model)
m1.LL<-logLik(fit.fert.dist)
m2.LL<-logLik(fit.fertxdist)
m1.df<-attr(m1.LL,"df")
m2.df<-attr(m2.LL,"df")

lik.ratio <- -2*m1.LL[1]+2*m2.LL[1]
df.diff<-m2.df-m1.df

pchisq(lik.ratio,df=df.diff,lower.tail=F)
summary(fit.lm3)
# do these match?


# Or, we can check the significance of main effect of disturbance 
# (by comparing against the nested model containing only fertilization
# treatment and interaction between fertilization treatment and 
# disturbance treatment)
fit.fert.fertxdist<-mle2(log.biomass~dnorm(mean=b0+b1*fert.code+b2*fert.code*dist.code,sd=s),start=list(b0=4,b1=0,b2=0,s=1),data=c1)
summary(fit.fert.fertxdist)

lik.ratio <- -2*logLik(fit.fert.fertxdist)[1]+2*logLik(fit.fertxdist)[1]
pchisq(lik.ratio,df=1,lower.tail=F)

# compare to p-value for dist.code:
summary(fit.lm3)


####__Exercise__####

# Can you compute the LRT estimate for the p-value of the main effect of fertilization?



# Note that the mle/likelihood ratio test estimates of the p-values 
# are slightly lower than those returned by the frequentist linear model

# Why might this be?


#----------------------------------------------------------------#
#----------------------------------------------------------------#

####  General linear models & ANCOVA ####

# We've now covered 2 different kinds of deterministic models: 
#	- regression type models where we determine the relationship between 
#     a continuous covariate and variation in our response variable
#	- ANOVA like models where we determine the effects of categorical 
#     covariates on variation in our response variable.

# These kinds of models are often presented as fundamentally different
# things (we even given them different names) - but they are 
# fundamentally connected. Additionally, there's no problem with mixing
# together both continuous and discrete covariates in the same model 
# (sometimes referred to as an ANCOVA).  This family of models are 
# collectively refered to as 'general linear models'.


# For fun and jollies, let's try out an mle2 analog of an ANCOVA

# Don't let the name confuse you:
# - an 'ANCOVA' or 'analysis of covariance' is really just a linear 
#   model with both categorical and continuous covariates (like a hybrid
#   between regression and ANOVA)

# Let's look at how fertilizer and time might both be affecting above 
# ground biomass in this old field succession experiment.

# We've already looked at the effect of fertilizer on biomass, but what
# do you think the trend in biomass might be over time? (if there is one!)


# Exploratory analysis
plot(log.biomass~centered.year,data=c1)

# Is that what you were expecting?  Any insight from the ecologists/plant
# enthusiasts on causes of this pattern?


# Note: we're essentially doing a regression on year and trying to estimate
# a slope coefficient. Numbers like 1998 are large numbers, so to make the 
# slope coefficient a little bit more interpretable, instead of using the 
# original year value, we're using centered year (which is just the difference
# between year and mean(year) in the data set).



####__Exercise__####

# Can you assemble an mle model using fertilization treatment
# and time (centered.year) as covariates? Start by using a linear 
# relationship with time.



# For comparison, here are results from a standard general linear model:
fit.lm4<-lm(log.biomass~fertilized*centered.year,data=c1)
summary(fit.lm4)

# plot results...
cfs<-coef(fit.lm4)

plot(log.biomass~centered.year,col=c1$fert.code+1,data=c1)
abline(cfs[[1]],cfs[[3]],col='red',lwd=2)
abline(cfs[[1]]+cfs[[2]],cfs[[3]]+cfs[[4]],lwd=2)
legend(5,7.5,legend=c('fert data','no fert data','fert trend','no fert trend'),lty=c(0,0,1,1),col=c(2,1,2,1),pch=c(1,1,-1,-1),lwd=c(1,1,2,2))
# Here we see declining linear relationships between log.biomass and 
# time, with different means between fertilization treatment (black 
# versus red).


####__Exercise__####

# Earlier, we assumed that biomass changes linearly with time. 
# How valid do you think this assumption is? Traditional approaches using 
# lm() require this assumption of linearity.  Can you put together an 
# mle model containing a quadratic relationship with time (still technically 
# 'linear' from a statistician's perspective)? How about a truly nonlinear 
# function? How do these models compare to the 'linear' model?


# We'll hopefully have a chance to talk a bit more about working with 
# time series data - there are a lot of other approaches. But, if you
# can come up with a valid, deterministic function that you think might
# explain how your response variable changes through time, there's every
# reason to try fitting it to your data. In other words, time is a 
# legitimate continuous covariate.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Beyond the mean: deterministic models for variance ####

# Thus far, every model that we've examined using a deterministic model 
# has related covariates to the means of our stochastic distributions.  
# But we're not limited to this in any way - we could just as easily 
# use covariates to explain other parameters describing our stochastic 
# distribution, such as the variance.  This can be really useful in 
# situations where the level of variance in a data set changes over time.

# Let's look at an example:

### Calculate variance patterns over time

# Let's return to a linear model of log.biomass over time:
fit.time<-mle2(log.biomass~dnorm(mean=b0+b1*centered.year,sd=s),start=list(b0=4,b1=-0.05,s=1),data=c1)
summary(fit.time)

# This model assumes that the variance around the mean relationship
# between log.biomass and time remains constant (an assumption of all
# basic general linear models - regression, ANOVAs, ANCOVAs)

# Let's check out what the variation looks like around our predicted 
# values through time to see how reasonable this assumption is:

# Unique years to obtain predictions from our model fit.time for:
u.years<-sort(unique(c1$centered.year))
new.data<-data.frame(centered.year=u.years)
new.ys<-predict(fit.time,new.data)

# Show the predictions on top of the data:
plot(log.biomass~centered.year,data=c1)
points(new.ys~new.data$centered.year,pch=19,col='blue')
abline(fit.time,col='blue')

# Variance is calculated as the sum of the squared differences between
# each data point and the mean of all of the data points, all divided 
# by (n-1), where n is the number of data points.

# So for our first year, this would look like:
first.year.biomasses<-c1$log.biomass[c1$centered.year==u.years[1]]

# We can apply the typical formula for estimating variance:
sum((first.year.biomasses-mean(first.year.biomasses))^2)/(length(first.year.biomasses)-1)

# Compare with the built-in R function:
var(first.year.biomasses)

# We're going to do a very similar calculation, but instead of using 
# mean(first.year.biomasses), we'll substitute the mean biomass estimated
# in the first year from our regression model, fit.time.

# Grab out all of the data from the first year
tmp<-c1[c1$centered.year==u.years[1],]
head(tmp)

# Calculate the local variance around our predicted mean in the first year
local.var<-sum((tmp$log.biomass-new.ys[1])^2)/(length(tmp$biomass)-1)
local.var
sqrt(local.var)	# corresponding standard deviation

# We can repeat this calculation, determining the value of the standard
# deviation of log.biomass values around the means predicted by our model
# for each year:
est.sd<-rep(NA,length(u.years))		        # data storage
for(i in 1:length(u.years)){		          # loop through years
	tmp<-c1[c1$centered.year==u.years[i],]	# grab out subset of the data
	local.var<-sum((tmp$log.biomass-new.ys[i])^2)/(length(tmp$biomass)-1)	# calculate local variance
	est.sd[i]<-sqrt(local.var)		          # calculate and save our estimated standard deviations
}

# Visualizing the results:
plot(est.sd~u.years)

# Compare the mean of these values with the value for residual standard 
# deviation from our mle fit:
mean(est.sd)
coef(fit.time)[3]

# We see that the estimated standard deviation around our model 
# predictions decreases through time.

# This can happen for several reasons:
#	- variation in the experiment decreases through time (perhaps as 
#     community assembly/succession causes the plant community to converge
#     on some equilibrium state or climax community)
#	- our model does bad job of describing the mean through time, changing 
#     our estimate of residual variation.

# What do you think happens with our data set?


# We can fit models where we describe how the residual variation changes
# through time using a deterministic relationship:
fit.time.var<-mle2(log.biomass~dnorm(mean=b0+b1*centered.year,sd=exp(s0+s1*centered.year)),start=list(b0=4,b1= -0.05, s0=0.5,s1=0),data=c1)
# here the standard deviation of our stochastic normal distribution can 
# increase or decrease through time according to an exponential function 
# (why exponential?)

# check out the results:
summary(fit.time.var)

# model comparison
AICtab(fit.time,fit.time.var)

# visualize the results
plot(est.sd~u.years)
abline(coef(fit.time)[[3]],0)	# old constant estimate
curve(exp(coef(fit.time.var)[[3]]+coef(fit.time.var)[[4]]*x),-11,12,add=T,col='red') 
# new estimate w/ exponential relationship.



####__Exercise__####

# Investigate modeling the residual variance in log.biomass
# using something other than a linear trend between biomass and time. 
# Are there still differences in variance over time?




#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Run an ANOVA with more than two levels ####

# The LTER data on plant biomass through time was useful for looking at 
# working with categorical covariates, but each of its treatments had 
# only two levels. It's quite common for data to have more than two 
# categories for any given categorical covariate, and this can introduce
# some extra complexity to setting up models.

# For the next example, we'll use data on frog predation explored in the
# EMD book (pg 263 and beyond). The relevant study looked at the effects
# of frog size, the density of frogs, and the presence of predators, on 
# frog survival.

# load the package containing the data we'll use
library(emdbook)

# grab data on frog predation
data(ReedfrogPred)

# For now, let's just consider the subset of the data where predators were present:
pred.data<-subset(ReedfrogPred,pred=="pred")
pred.data$density<-as.factor(pred.data$density)
pred.data

# Our response variable is the proportion of frogs surviving predation. 
# Modeling proportional data appropriately takes some extra consideration, 
# but for now we'll assume:
# - that the stochastic portion of our model follows a normal distribution.
# - that some of the variation in survival is attributable to a deterministic
#   relationship between frog density and survival.

# Frog density has 3 different levels, 10, 25, and 35:
pred.data$density

# Because there are three levels, the coding that we need to test for 
# the effects of different categories is more complicated (not just one 
# column of 0's and 1's).

# Set up the model matrix
col1<-ifelse(pred.data$density=='25',1,0)
col2<-ifelse(pred.data$density=='35',1,0)
mod.matrix<-cbind(col1,col2)

### There are several different ways you could set up this model, but 
# the basic idea is:

# mean = p0 + p1*col1 + p2*col2

# when density = 10, col1 and col2 both equal 0, so our estimate of the mean is mean = p0
# when density = 25, col1 = 1 and col2 = 0, so our model becomes mean = p0 + p1
# when density = 35, col1 = 0 and col2 = 1, so our model becomes mean = p0 + p2

# Interpretation: p0 is the mean proportion surviving given a density of 10. 
# p1 is difference(or change) in the mean proportion surviving between a 
#     density of 10 and a density of 25. 
# p2 is the difference (or change) in the mean proportion surviving between
#     a density of 10 and a density of 35.


### Here are 3 different ways of structuring this model:

# (1) Use each individual column (most intuitive)
fit.dense1<-mle2(propsurv~dnorm(mean=p0+p1*col1+p2*col2,sd=s),start=list(p0=0.5,p1=0,p2=0,s=0.5),data=pred.data)
summary(fit.dense1)

# (2) Use a little matrix multiplication
fit.dense2<-mle2(propsurv~dnorm(mean=p0 + mod.matrix %*% c(p1,p2),sd=s),start=list(p0=0.5,p1=0,p2=0,s=0.5),data=pred.data)
summary(fit.dense2)

# (3) Use the parameters option of mle2()
fit.dense3<-mle2(propsurv~dnorm(mean=p0,sd=s),start=list(p0=0.5,s=0.5),parameters=list(p0~pred.data$density),data=pred.data)
summary(fit.dense3)

# Results are identical for all 3 approaches.

# Making more and more complicated models that includeg categorical 
# covariates simply involves setting up the appropriate model matrices.



####__Bonus Exercise__####

# This is an interesting and fairly complex data set. Explore it more - perhaps 
# treating density as a continuous covariate instead of categorical/discrete. 
# You can check out Ben Bolker's analysis in the EMD Book, starting around page 263





#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Recap ####

# In this lab script, we've covered the following:
#
#	- fitting deterministic models with categorical covariates to LTER plant biomass data
#	- Likelihood Ratio Tests for estimating p-values
#	- One- and Two-way ANOVAs
#	- mixing categorical and continuous covariates (ANCOVAs) 
#	- modeling parameters other than the mean: beyond constant/equal variance assumptions
#	- categorical covariates with more than two levels

#----------------------------------------------------------------#
#----------------------------------------------------------------#


#### 'Answers' ####

fit.fert.lnorm<-mle2(biomass~dlnorm(meanlog=exp(b0+b1*fert.code),sdlog=exp(s)),start=list(b0=log(4),b1=0,s=1),data=c1)
summary(fit.fert.lnorm)

# no fertilizer, mean log biomass:
exp(coef(fit.fert.lnorm)[[1]])

# fertilizer, mean biomass
exp(coef(fit.fert.lnorm)[[1]])*exp(coef(fit.fert.lnorm)[[2]])

# difference:
exp(coef(fit.fert.lnorm)[[1]])*exp(coef(fit.fert.lnorm)[[2]])-exp(coef(fit.fert.lnorm)[[1]])

# standard error
exp(coef(fit.fert.lnorm)[[3]])





#----------------------------------------------------------------#
#----------------------------------------------------------------#

###  Scraps:

library(emdbook)

# shortcut for calculating likelihood slices:
?calcslice()

fit.fert.lw<-mle2(log.biomass~dnorm(mean=b0+b1*fert.code,sd=s),start=list(b0=4,b1=0,s=1),fixed=list(b1=0),data=c1)
fit.fert.up<-mle2(log.biomass~dnorm(mean=b0+b1*fert.code,sd=s),start=list(b0=4,b1=0,s=1),fixed=list(b1=2),data=c1)
res<-calcslice(fit.fert.lw,fit.fert.up)
critical.value<- -logLik(fit.fert)[[1]]+1.92

plot(res$y~I(2*res$x),type='l',xlab="b1",ylab="NLL")
abline(critical.value,0,lty=3)

summary(fit.fert)


### Exploratory analyses

plot(log.biomass~centered.year,data=c1)

fit.lin<-lm(log.biomass~centered.year,data=c1)
summary(fit.lin)
abline(fit.lin,lwd=2)

lw<-lowess(c1$centered.year,c1$log.biomass,f=0.2)
lines(lw,col='blue',lwd=2)

lw.fert<-lowess(c1$centered.year[c1$fertilized=="fert"],c1$log.biomass[c1$fertilized=="fert"],f=0.2)
lw.nofert<-lowess(c1$centered.year[c1$fertilized=="no.fert"],c1$log.biomass[c1$fertilized=="no.fert"],f=0.2)
lines(lw.fert,col='red',lwd=2)
lines(lw.nofert,col='red',lty=2,lwd=2)

# might want to compress the temporal resoultion into 'year' to avoid 
# small-scale variation in harvest date within season

summary(c1$log.biomass)

fit.lin2<-lm(log.biomass~centered.year*fertilized,data=c1)
summary(fit.lin2)

