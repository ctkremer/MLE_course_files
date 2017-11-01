#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Lab 5: Logistic regression using MLE
# ELME MLE course
# Colin T. Kremer

# Logistic regression is a specific kind of 'generalized linear model',
# see Bolker's EMD, chapter 9, pages 308-311 for more.

# Another good reference for generalized linear models is 
# 'Data analysis using regression and multilevel/hierarchical models'
# by Andrew Gelman and Jennifer Hill, especially chapters 5 & 6

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Load data ####

# Libraries
library(bbmle)

# load mle tools
source("/Users/colin/Teaching/ELME/ELME 2017/labs/mle.tools_070717.R")

# load data
data<-read.csv("/Users/colin/Teaching/ELME/ELME 2017/labs/Lab5 - generalized linear models/cyanobacteria_data_053114.csv")
data<-na.omit(data)

head(data)

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Exploratory graphical analyses ####

# Take a look at how presence/absence of nitrogen fixing cyanobacteria 
# depends on temperature.
plot(pa~Temperature,data=data)
# What can you conclude from this plot? Anything?


# It helps to look at these data a little differently. I've written a tool
# for this purpose, called regplot. You can read more about it in mle.tools_052815.R
# Basically, it splits the x axis (covariate) up into n little intervals or bins, 
# and calculates how many 1's occur within that interval relative to the number of
# total observations in that interval:  (# of 1's)/(# of 1's and 0's) these values
# provide an estimate of the probability of observing a 1 over a small range of the
# x axis. We're interested in determining whether this probability changes over the
# whole range of the x axis. By repeatedly calculating these local estimates, then 
# plotting them, we can look for patterns:
regplot(data$pa,data$Temperature,30)

# The local 'probability of occurrence' estimates are shown as red triangles.

# What can you conclude from this new plot? Do you think water temperature has an
# effect on the probability that nitrogen fixing cyanobacteria are present in a lake?

# Now let's get quantitative about this.

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Logistic regression ####


# What stochastic distribution/process do you think might be appropriate for 
# modeling this kind of data?


# We've got discrete data, with values that are either 0 or 1, and we're interested
# in understanding how the probability of observing a 1 versus a 0 might change 
# with environmental factors.

# Yup, you probably guessed it - we're looking at data that could come from a
# binomal process/distribution. 
# Remember the parameters we need?
#   n = number of trials (1 in this case), also called 'size' in R
#   p = probability of observing a success (1) on each trial, also called 'prob' in R
?dbinom()


# Let's make a null model as a starting point:
m0<-mle2(pa~dbinom(size=1,prob=p0),start=list(p0=0.5),data=data)
summary(m0)
# What does it mean that p0 = 0.44899?


####__Exercise__####

# Let's try to add in a covariate now; we'll focus on temperature first, and 
# similar to our work in previous labs, let's assume that temperature has a linear
# effect on probability of presence.
mT.1<-mle2(pa~dbinom(size=1,prob=___?____),start=list(p0=0.5, p1= ),data=data)
summary(mT.1)


### Uh-oh. Why are we getting so many error messages?
# ('In dbinom(x = c(1L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 0L,  ... : NaNs produced')


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Link functions ####

# A linear model of the form:  prob = p0 + p1 * Temperature
# maps Temperature (a variable that theoretically can range from -Infinity to Infinity)
# to 'prob', which also ranges from -Infinity to Infinity in this case. BUT, what 
# does it mean to have a probability less than zero, or greater than one? Such 
# values are non-sensical. Fundamentally, that's what causes the mle2() model above
# to freak out, unless you had just the right initial guesses.

# We can use a different functional form for the deterministic part of this model, 
# other than the equation of a line.

# A traditional choice for such a function looks like this:

#   prob = inverse.logit(p0 + p1*Temperature)
#   prob = exp(p0 + p1 * Temperature)/(1 + exp(p0 + p1 * Temperature))

# This approach relies on using something called the 'logit' link function. 
# There's a easy function defined in mle.tools_052914.R for this, called expit(), 
# which is equivalent to the inverse of the logit function:
# expit < -function(x){exp(x)/(1+exp(x))}
expit(1)

# It looks like this:
curve(expit(0.2 + 2.1*x),-2,4)

# One of the useful properties of this function is that no matter what value we use 
# for x, expit(x) will always return a value between 0 and 1. Try for yourself:
expit(100)

# If we use this function instead of the equation for a line we tried in our 
# previous model, we'll be safe - we'll never get predictions for prob that are
# negative or greater than one. 


#~~~ Digression ~~~#

# We've talked about different kinds of deterministic functions that can be used
# in statistical models. The one I just introduced, expit() might look like a 
# nonlinear function (when we plot it, it's not drawing a line, right?)

# This is one of those areas where terminology gets murky. From the perspective 
# of Statisticians, these models are still considered to be 'linear' even though
# we've used the expit() function, because the equation we imbed within a given 
# link function is still linear in terms of our covariates (X's)

#~~~ End digression ~~~#


# Let's try fitting our presence/absence data using this new deterministic function
mT.2<-mle2(pa~dbinom(size=1,prob=expit(p0+p1*Temperature)),start=list(p0=0.5, p1=0),data=data)
summary(mT.2)

# Now we're in business! No warnings.

# Visualize result:
cfs<-coef(mT.2)
regplot(data$pa,data$Temperature,20)
curve(expit(cfs[[1]]+cfs[[2]]*x),2,35,add=T,col='red')


####__Exercise__####

# Assess this model:
AICtab(___?____)
anova(___?___)

# What do you think? Does this seem like a good fit? What can you determine 
# about the biological processes that might be occuring in these lakes? (ie, how 
# would you interpret the results?)


####__Exercise__####

# Try looking into some of the other covariates present in this data set 
# (doy, pH, Nitrogen (lwTN) and Phosphorus (lwTP))

# Can you assemble models with multiple covariates? What about interactions?

# Use model comparison approaches to find the best model for explaining the
# presence of cyanobacteria.


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Other kinds of link functions ####

# There are lots of different functions that have the property of mapping
# values from (-Inf, Inf) to (0,1)

# Here are a few others:
# initialize various functions mapping Reals to (0,1) 
loglog<-function(x){exp(-exp(-x))}
cloglog<-function(x){1-exp(-exp(-x))}
cauchy<-function(x){0.5+(1/pi)*atan(x)}


####__Exercise__####

# Explore how these functions differ from each other (make plots?)

# Run some mle2() models using different link functions; do any of them 
# perform better than the expit() function we used above?



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Fitting Generalized Linear Models in R ####

# We've been using mle2() to fit our logistic regressions, but I want to point
# out a perhaps easier (if less clear) way of fitting the same models using 
# default functions in R.

# Recall from above:

# Probability of presence as a function of temperature
mT.2<-mle2(pa~dbinom(size=1,prob=expit(p0+p1*Temperature)),start=list(p0=0.5, p1=0),data=data)
summary(mT.2)

# Fit the same model using glm() and specifying a binomial distribution & link function
glmT<-glm(pa~Temperature,family=binomial(link=logit),data=data)
summary(glmT)

# How do the parameter estimates compare to each other?
AICtab(glmT,mT.2)

# Everything should be identical; in fact, I glm() is invoking a maximum likelihood
# approach behind the scenes.

# glm() can be used for a variety of different distributions and link functions, see
# the help file for more. It only works on models that are statistically 'linear', 
# and can't deal with some of the crazier distributions, so it's still important to 
# know how to build your own models using something like mle2().



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Beyond binomal data ####

# We've been looking at explaining a discrete y-variable (or dependent variable)
# with only 2 possible values, 0 and 1

# It's possible to model data with more than two outcomes; for example, if you study
# behavioral ecology you might watch an animal's behavior at different points in time,
# recording whether it was sleeping, foraging, or courting a mate.

# These represent 3 different discrete possibilities.

# We can model this kind of data using an extension of the binomial distribution
# called the multinomial distribution.

# Where the binomial distribution is analogous to flipping a coin repeatedly 
# (outcome is heads or tails), the multinomial distribution is analogous to rolling
# a K-sided dice N number of times, given the probability of the dice landing on 
# each side is prob = c(p1, p2, ..., pK)

# When K = 2 and N = 1, the multinomial distribution reduces to the binomial 
# distribution (there are K = 2 possible outcomes, 0 or 1, we flip a coin once as N = 1).

?dmultinom()

# generate random data from multinomial distribution
rmultinom(10, size = 12, prob = c(0.1,0.2,0.8))

# We can write models for mle2() for multinomial data, although the link function
# becomes a little more complicated.

# I'll develop this more before class if I have time, but it's a fairly specialized
# application. Could be good for a group project, depending on the kind of data 
# people bring.

# See http://en.wikipedia.org/wiki/Multinomial_logistic_regression

# Also the Gelman and Hill book


#----------------------------------------------------------------#
#----------------------------------------------------------------#
