#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Lab 2: Fire Part 1
# ELME MLE course
# Colin T. Kremer

# written 2012, updated 5/2014, 5/2015, 7/2017

# Reading:
#	Grace and Keeley 2006, (Grace2006.pdf on wiki)

# Data:
#	Keeley_rawdata_ELME.csv (on wiki)

# Topics covered:
#	- Introduction to Grace and Keeley 2006, species richness/fire data set
#	- saw some different ways of visualizing data (exploratory data analysis)
#	- revisited differences between likelihood, log likelihood, and negative 
#     log likelihood.
#	- Continued introduction to the bbmle package/mle2() function
#	- Fitting various probability distributions to the observed distribution 
#     of species richness in Grace and Keeley 2006.
#	- Investigating results of model fits
#	- Introduction to model comparison

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Getting started ####

# Tell R what folder you're working in (aka 'set your working directory')
# this makes it easy to find data files and other items located in the same folder
# and is the default location where any output you save will be kept.
setwd("/Users/colin/Teaching/ELME/ELME 2017/labs/Lab2 - stochastic and deterministic model components/")

# If we're going to work with data, we some how have to get it from the file and format that it's stored in, and suck it over into R where we can manipulate it and analyze it.  My preference is to work with .csv (comma separated value) files - easily created in excel - but R can read in any number of other kinds of file formats, from text files using various delimiters (the name given to characters seperating entries in a spreadsheet), to more specialized files like netCDF files (sometimes requiring additional packages).
# see also ?read.csv()
data<-read.csv("Keeley_rawdata_ELME.csv")

# This data file can be found on the class site.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Data description ####

# See Grace2006.pdf on the wiki > Lab activities

# from the abstract: "This study investigates patterns of plant diversity following wildfires in fire- prone shrublands of California, seeks to understand those patterns in terms of both local and landscape factors, and considers the implications for fire management. Ninety study sites were established following extensive wildfires in 1993, and 1000-m2 plots were used to sample a variety of parameters."


#### Graphical exploration ####

# get a sense of the structure and layout of this data
head(data)
str(data)


# In general, it can be a really good idea to explore your data by plotting it in a variety of ways.

# For now, we'll focus on just one variable: species richness, the number of species present in a location.
# This is the primary response variable investigated by Grace and Keeley 2006 (the source of this data).

# First, we'll look at histograms, showing the frequency with which different levels of species abundance occured in the sampled plots. Histograms are a great way to get a sense of the variation that exists in a particular variable. This is typically the kind of variation that we try to explain with statistical modeling.
hist(data$rich,col='blue')

# Histograms group (or bin) together similar values. It's important to look at the data using a variety of different bin sizes, otherwise the visuals can be deceptive.
hist(data$rich,30,col='blue')	# change the number of bins
hist(data$rich,30,col='blue',xlab='Species richness')	# add a label to the x axis

# Another way of looking at this data is to plot the density (essentially a continuous smoothing of the histogram)
plot(density(data$rich),xlab='Species richness',lwd=2,ylim=c(0,0.06),main="")

# Important to note that species richness is a discrete variable, so density plots can be deceptive as they produce a smoothed output that appears to be continuous.
# It's also important to play around with the 'bandwidth' setting 'bw' - analogous to changing the bin size in a histogram, for example:
plot(density(data$rich,bw=2),xlab='Species richness',lwd=2,ylim=c(0,0.06),main="")


# A third way of visualizing variation is a box and whisker plot.
boxplot(data$rich,ylab='Species richness')

# boxplots can be created to show variation in a focal variable, as grouped by a covariate.
boxplot(data$rich~data$age,xlab='Age',ylab='Species richness')

# - not a terribly useful example for this data set, I just wanted to show you how to generate multiple boxplots


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### What kind of distribution best describes variation in species richness? ####

# So, we've observed that there's a lot of variation in species richness between sites.  We'd like to quantify this variation, and understand what mechanisms drive it.

library(bbmle)	# load up our mle tools.
source("/Users/colin/Teaching/ELME/ELME 2017/labs/mle.tools_070717.R")

# Steps in using maximum likelihood estimation to fit a probability distribution to data:

#	1) Pick a distribution
#	2) Write out the negative log likelihood function based on this distribution
#	3) Select starting guesses for parameter values
#	4) Use numerical optimization to calculate parameter values with the maximum likelihood given the data (or solve analytically, where possible).


###	1) Pick a distribution

# These data are discrete (number of species present). 
# - What kinds of probability distributions come to mind?






# Since the data are discrete, and all values are >=0, let's start 
# by trying the poisson distribution might be good to use


###	2) Write out the negative log likelihood function based on this distribution

# calculate the likelihood of a particular estimate of the poisson parameter, 
# the mean or lambda, given the observed data.

# poisson likelihood function
poissonL<-function(lambda){
	prod(dpois(data$rich,lambda=lambda))
}

# poisson log likelihood function
poissonLL<-function(lambda){
	sum(dpois(data$rich,lambda=lambda,log=T))	
}

# poisson negative log likelihood function
poissonNLL<-function(lambda){
	-sum(dpois(data$rich,lambda=lambda,log=T))	
}

# Now let's see what these functions look like.

lambda.guesses<-seq(40,60,1)	# lambda values to try
lik<-sapply(lambda.guesses,poissonL)	# determine the corresponding likelihood
log.lik<-sapply(lambda.guesses,poissonLL)	# determine the corresponding log likelihood
neg.log.lik<-sapply(lambda.guesses,poissonNLL)	# determine the corresponding negative log likelihood

#~~~ Digression ~~~#

# What is this crazy sapply function?

# sapply(vector, function.to.apply.to.vector)

# Short answer - see R crash course. Basically, it functions like a simple for-loop,
# applying the specified function to each value in the given vector and returning 
# the results. 

#~~~ End Digression ~~~#

# Visualize the difference between likelihood, log likelihood and negative log likelihood
par(mfrow=c(2,2))		# set up graphing window
plot(density(data$rich,bw=4),xlab='Species richness')
plot(lik~lambda.guesses,main='likelihood',xlab='lambda',type='b')
plot(log.lik~lambda.guesses,main='log likelihood',xlab='lambda',type='b')
plot(neg.log.lik~lambda.guesses,main='negative log likelihood',xlab='lambda',type='b')



###	3) Select starting guesses for parameter values

# How do we come up with a good way of estimating initial values for the 
# parameters of our distribution? There are several different strategies:
#	- plot the data and eyeball it
#	- use the distribution_demos.cdf toolbox
#	- method of moments estimator (will discuss in lecture)
#	- try a whole bunch of values and see what happens.

# For this data, try out the distribution_demos.cdf tool, 
# but we'll stick with the 'method of moments' estimate, in this case just the
# mean of the data:
mean(data$rich)


###	4) Use numerical optimization to calculate parameter values with the maximum 
# likelihood given the data (or solve analytically, where possible).

# Now, to fit the poisson distribution to the data, we'd like to employ a numerical 
# optimization algorithm that will help us find the maximum (minimum) of this log 
# likelihood (negative log likelihood) surface. This should be familiar from our 
# earlier lab activity.

# The primary optimizer that we'll use during this course is called mle2() and is 
# found in the bbmle package.
# - technically it actually calls another, more general function called optim() 
#   (sound familiar?) which does the real optimization, but mle2() has a lot of 
#   nice features that make it convenient to use and easier than talking directly 
#   to optim() ourselves.

# There are a couple of different ways that you can send models to mle2(), which 
# we'll explore now.

#### Method 1: Provide the name of your negative log likelihood calculator

fit.pois.1<-mle2(poissonNLL,start=list(lambda=mean(data$rich)),data=data)

summary(fit.pois.1)

# Maximum likelihood estimation
#
# Call:
# mle2(minuslogl = poissonNLL, start = list(lambda = mean(data$rich)), 
#     data = data)
#
# Coefficients:
#        Estimate Std. Error z value     Pr(z)    
# lambda 49.23333    0.73962  66.566 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#
# -2 log L: 939.4122 

# Checking out the summary of our mle fit, we can see several pieces of information:
#
# - Call: this reminds us of the information that we sent to mle2
#
#	- Coefficients: this is one of the most important sections; it tells us what 
#   the maximum likelihood estimate of each of our coefficients/parameters is. 
#   It also provides an estimate of the standard error associated with the mle 
#   coefficients and their 'significance' (more on this later).
#
#	- Finally, at the end, it shows us what the value of -2 log Likelihood is (also
#   sometimes called the deviance). Why -2*log(Likelihood) you may well ask? This 
#   turns out to be a really useful quantity that we'll run into later on when we 
#   talk about model comparison.



#### Method 2: Provide mle2() with a statistical formula of the form: y ~ distribution(x,p)

# Here we tell mle2 that we want to model some quantity y (our response variable) 
# according to some distribution, having parameter(s) p and potentially incorporating
# other data, x.

# In this case, we want to model species richness, so 'rich' from our data set ends
# up on the left hand side of the '~'. 

# We're investigating how well the Poisson distribution does at describing the 
# distribution of species richness, so our distribution on the right hand side 
# of the '~' is dpois().

# Inside of dpois() we have to provide the mean of the poisson distribution, lambda.
# In this case, we say that lambda = l, where l is coefficient that we want to use 
# maximum likelihood analysis to estimate.

# Given this information, mle2 is smart enough to convert this formula into a form 
# that it can use to calculate negative log likelihoods, so we don't have to write 
# a separate negative log likelihood calculator. 
 
fit.pois.2<-mle2(rich~dpois(lambda=l),start=list(l=mean(data$rich)),data=data)
summary(fit.pois.2)

# We can compare the results of this fit to the results from Method 1 above, and 
# see that the results are identical - we haven't changed the fundamental statistical
# model that we're using, just how we talk to mle2.

# We can also look directly at the parameter estimates
coef(fit.pois.2)
l.mle<-coef(fit.pois.2)[[1]]	# and save them for later use.


# If we can get away with not writing our own negative log likelihood calculators, 
# why would we ever bother to do so?
#	- writing out the explicit negative log likelihood calculator is useful 
#   pedagogically when you're first getting used to mle analyses, because it 
#   makes you think about what you're doing
#	- the negative log likelihood calculator is handy when you want to try to plot
#   your likelihood function/surfaces.
#	- sometimes, if you're working with customized distributions or complicated 
#   likelihoods, it can be difficult or impossible to squeeze everything in to 
#   the shorter statistical formulas, and it's helpful to have your likelihood 
#   calculator be a separate function.


### Remember our first big question:

# Q1: What kind of distribution best describes variation in species richness?

# So, we've determined by way of maximum likelihood analysis that the poisson 
# distribution that best fits our data on variation in species richness has a 
# mean (lambda) of 49.233.  This represents in effect the best job that the 
# poisson distribution can do at describing the data observed. (or, put another 
# way, the most likely value of lambda, given the data, is 49.233)

# It's comforting that this corresponds exactly with the mean of our original data:
mean(data$rich)

# But we might like to take this one step farther to ask the question 'This may be 
# the best that the poisson distribution can do, but is its best actually good enough?' 

# Does the poisson distribution describe the data well?

# Visualize
par(mfrow=c(1,1))
plot(density(data$rich,bw=4),xlab='Species richness',lwd=2,ylim=c(0,0.06),main='')
curve(dpois(x,lambda=l.mle),0,100,col='red',lwd=2,add=T)
legend(5,0.06,legend=c('data distribution','poisson distribution'),lty=1,
       col=c('black','red'),lwd=2,cex=0.5)

# What do you think?  Is this a reasonable/believeable fit?

# Thinking back to some of our lecture topics, and the properties of the poisson
# distribution, can you think of any reasons why the poisson distribution might 
# be having difficulty describing the distribution of species richness in our data?



# One of the strengths of the maximum likelihood approach, is the ease with which 
# we can generate multiple models (corresponding to multiple hypotheses), fit them
# to our data, and then compare them to each other to determine which model(s) 
# is preferred.

# For example, we could imagine perhaps that the species richness data might follow
# a normal distribution.  Indeed, for most traditional statistical tests/approaches, 
# we would be constrained to make this assumption, instead of investigating the 
# poisson distribution as we did above.

### Method 1: Explicit negative log likelihood calculator

# normal negative log likelihood calculator:
normalNLL<-function(mu,s){
	-sum(dnorm(data$rich,mean=mu,sd=s,log=T))	
}

fit.norm.NLL<-mle2(normalNLL,start=list(mu=mean(data$rich),s=sd(data$rich)),data=data)
summary(fit.norm.NLL)


### Method 2: provide statistical formula

fit.norm<-mle2(rich~dnorm(mean=mu,sd=s),start=list(mu=mean(data$rich),s=sd(data$rich)),data=data)
summary(fit.norm)

logLik(fit.norm)

coef(fit.norm)
mu.mle<-coef(fit.norm)[[1]]
sd.mle<-coef(fit.norm)[[2]]


# Visualize the poisson and normal fits to our data
plot(density(data$rich,bw=4),xlab='Species richness',lwd=2,ylim=c(0,0.06))
curve(dnorm(x,mean=mu.mle,sd=sd.mle),0,100,col='blue',lwd=2,add=T)
curve(dpois(x,lambda=l.mle),0,100,col='red',lwd=2,add=T)
legend(5,0.055,legend=c('data distribution','poisson distribution',
        'normal distribution'),lty=1,col=c('black','red','blue'),lwd=2,cex=0.5)


# How do these two models compare to each other?  Which do you think is better?

# What accounts for their differences?


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Formal model comparison ####

# So far we've been using graphs to look at our model fits, and to try to 
# visually assess how good of a job they do. Ideally though, we'd like to 
# have a much more quantitative, objective means of comparing models.

# We'll cover this more in lecture, but for now here's a brief introduction 
# to the essence of this approach:

# We can use the quantitative value of the maximum likelihoods associated 
# with the poisson and the normal models to compare them against each other 
# (in a sense, to find out how much more likely one model is than the other, 
# GIVEN THE SAME Y DATA). One way of doing this is to used Akaike Information 
# Criteria (or AIC) values. 

# In addition to using the information contained in a model's maximum likelihood,
# the AIC value also incorporates a penalty based on the number of parameters 
# that a model uses. In general, models with more parameters (or 'degrees of 
# freedom') are expected to do a better job of fitting data than models with 
# fewer parameters.  Thus, we might expect that models with more parameters should
# have a larger maximum likelihood than models with fewer parameters.  What we 
# really want to know is if the addition of these extra parameters leads to a 
# substantially better model fit, more than we might expect from the extra degrees 
# of freedom alone. This is a justification for adding a penalty based on the 
# number of parameters in a model

### Calculation/derivation of AIC values

# calculate the AIC given the negative log likelihood and the degrees of 
# freedom (= number of parameters)
# -2*NLL+2*number.of.parameters
-2*logLik(fit.pois.1)[[1]]+2*attr(fit.pois.1,"df")

# calculate the AIC directly
AIC(fit.pois.1)

# AIC for normal distribution model
AIC(fit.norm)

# compare these AIC values
AICtab(fit.norm, fit.pois.1)


# What does this mean?  How would you interpret the results?


# More on this in lecture, but we can also calculate 'Akaike weights'
# for each of these models. Sometimes we end up having multiple models 
# that do a pretty comparable job of explaining our data. One approach 
# for handling this situation, if we're trying to make predictions, is 
# to average together the predictions of each of the models in our set 
# of acceptable models. The Akaike weights help us average together models, 
# weighting them based on just how good they ar relative to the rest of 
# the models considered. See EMD page 215-216 for more.
AICtab(fit.norm, fit.pois.1, weights=T)



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#___ Exercise: ___#

# Can you investigate the negative binomial distribution as 
# an alternative to the poisson distribution and the normal distribution?

?dnbinom()
# check out the negative binomial distribution on wikipedia, and in the 
# distribution_demos.cdf tool kit.


# Perform an mle2 fit:

# 1) by way of writing a negative log likelihood calculator (Method 1)

# 2) directly, by writing out your statistical model inside of your call
#   to mle2() (Method 2)


# Compare this model to the normal and poisson models.
#	1) using AIC
#	2) visually

# What do you conclude?

# Which model is best, and why?

# Come together as a class and/or in smaller groups to discuss your conclusions.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Comments on kurtosis ####

# We mentioned earlier that comparing the maximum likelihoods of models 
# (by way of AIC) enables us to compare models against each other. This 
# allowed us to evaluate the performance of the normal, poisson, and 
# negative binomial models against EACH OTHER.  However, these results 
# are only relative - of a specific set of models, the AIC approach can 
# only tell us how good any given model is in comparison to the other 
# models in the set. We never know if our set of models actually contains
# the 'true' best model (if you even think such a model exists). In other 
# words, AIC comparison doesn't let us get a sense of the absolute fit of 
# the model to the data.

# This is a pretty important philosophical (and practical) issue that we'll
# return to.


# Looking at the plots showing model fits, it's pretty clear that the normal 
# and negative binomial distributions do the best job of of describing the 
# distribution of species richness (compared to the poisson distribution).  
# However, the real distribution of species richness is what we would call
# 'platykurtic' 

# < Digression: Find platypus picture... >
# course website > MLE_files > Lab 2 > Student 1927 platykurtic platypus.pdf


# In other words, it has 'shoulders' - the top of the distribution is 
# pretty flat, without a strong mode  (value that is much more commonly 
# observed than surrounding values). As described by Student 1927, 
# platykurtic distributions also have much shorter 'tails', ie, less 
# probability density associated with extreme values than a normal distribution.

# This often suggests that there are other forces at work driving the 
# observed patterns of variation. In this case, we might hypothesize 
# that factors other than just stochasticity are affecting the species 
# richness that we observe (probably this is no surprise to all the 
# ecologists in the group)

# We'll investigate this more in the next lab, but if you're curious, 
# think about relevant mechanisms, and/or check out Grace2006.pdf


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Bonus exercise ####

# Another important variable measured in this data set is the severity 
# of wildfires. Wildfires can be highly variable, and their severity, 
# heat intensity, timing and spatial extent can all affect the richness 
# and composition of plant communities.

# Can you examine wildfire severity (data$firesev) in this data set? 
# What distributions might describe variation in fire severity? You may 
# need to reference Grace and Keeley 2006 (Grace2006.pdf on wiki) to 
# understand how fire severity was quantified. Can you fit various 
# distributions to the empirical distribution of fire severity? Use 
# distribution_demos.R to select potential distributions to try, or 
# reference EMD book page 120.



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Recap ####

# In this lab script, we've covered the following:

#	- Introduction to Grace and Keeley 2006, species richness/fire data set
#	- saw some different ways of visualizing data (exploratory data analysis)
#	- revisted differences between likelihood, log likelihood, and negative 
#     log likelihood.
#	- Continued introduction to the bbmle package/mle2() function
#	- Fitting various probability distributions to the observed distribution 
#     of species richness in Grace and Keeley 2006.
#	- Investigating results of model fits
#	- Introduction to model comparison

#----------------------------------------------------------------#
#----------------------------------------------------------------#






##### "Answers"


# For negative binomial

### Method 1: Explicit negative log likelihood calculator

# negative log likelihood calculator:
nbmNLL<-function(mu,s){
  -sum(dnbinom(data$rich,size=s,mu=mu,log=T))	
}

fit.nbm.NLL<-mle2(nbmNLL,start=list(mu=mean(data$rich),s=sd(data$rich)),data=data)
summary(fit.nbm.NLL)


### Method 2: provide statistical formula

fit.nbm<-mle2(rich~dnbinom(mu=mu,size=s),start=list(mu=mean(data$rich),s=sd(data$rich)),data=data)
summary(fit.nbm)

AICtab(fit.nbm.NLL,fit.nbm)

mu2.mle<-coef(fit.nbm)[[1]]
size.mle<-coef(fit.nbm)[[2]]


# visualize the poisson and normal fits to our data
plot(density(data$rich,bw=4),xlab='Species richness',lwd=2,ylim=c(0,0.06))
curve(dnorm(x,mean=mu.mle,sd=sd.mle),0,100,col='blue',lwd=2,add=T)
curve(dpois(x,lambda=l.mle),0,100,col='red',lwd=2,add=T)
curve(dnbinom(x,mu=mu2.mle,size=size.mle),0,100,col='purple',lwd=2,add=T)
legend(0,0.062,legend=c('data distribution','poisson distribution','normal distribution','negative binomial distribution'),lty=1,col=c('black','red','blue','purple'),lwd=2,cex=0.5)

AICtab(fit.nbm,fit.norm,fit.pois.2)


### Bonus activities


bfit.norm<-mle2(firesev~dnorm(mean=mu,sd=exp(s)),start=list(mu=4.5,s=1),data=data)
bfit.lnorm<-mle2(firesev~dlnorm(meanlog=mu,sdlog=exp(s)),start=list(mu=4.5,s=1),data=data)
bfit.gamma<-mle2(firesev~dgamma(shape=mu,scale=exp(s)),start=list(mu=4.5,s=1),data=data)

AICtab(bfit.norm,bfit.lnorm,bfit.gamma)








