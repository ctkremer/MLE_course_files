#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Lab 1: Distributions, random data, likelihood calculations & numerical optimization
# ELME MLE course
# Colin T. Kremer

# written 2012, updated 5/2014, 5/2015, 7/2017

# Reading:
#	EMD
#	- Chpt. 4: Probability distributions
#	- Chpt. 6: Likelihood calculations
#	- Chpt. 7: Optimization

# Topics covered:
#	- the basics of working with probability distributions in R
#	- how to generate random data (from a poisson distribution)
#	- how to use maximum likelihood methods to estimate the 
#	  parameters generating that random data.
#		- the relationship between likelihood, log likelihood, 
#		  and negative log likelihood functions
#	- an introduction to the bbmle package/mle2() function
#	- (dis-)advantages of numerical optimization
#	- and an introduction to various numerical optimization algorithms

#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Probability distributions in R ####

# Probability functions:

# R has a ton of different probability distributions, many of which will be useful to us in 
# this class. Each distribution has a number of R functions associated with it.

# For example, let's think about the poisson distribution. There's a whole family of useful
# functions related to the poisson distribution, and in R the names of these functions all
# end with _pois() - we get different values depending on the character that we put in the 
# blank space _ in front of _pois():

# A random number from a poisson distribution:
# n=1 tells the function we only want 1 random number.
# lambda=5 tells the function what we want the mean of the poisson distribution to be.
rpois(n=1,lambda=5)

# The probability of observing a particular value from a poisson distribution:
# (in the context of a continuous distribution, this gives us the probability density instead)
# x=6 supplies the value of the particular observation (6) that we're asking for the probability of
# lambda=5 tells the function what we want the mean of the poisson distribution to be.
dpois(x=6,lambda=5)

# What is the probability of observing a value less than or equal to (OR greater than or equal to)
# a particular value, given a poisson distribution?

# Here q=6 provides the critical value; This function calculates the probability of observing 
# a value less than or equal to q, if we specify lower.tail=TRUE, OR the probability of 
# observing a value greater than or equal to q given lower.tail=FALSE

# As before, lambda=5 provides the mean of this particular poisson distribution.
ppois(q=3,lambda=5,lower.tail=TRUE)

# When our distribution is discrete, this is the same as adding up the probabilities of 
# each possible value up to the value of q:
dpois(x=0,lambda=5)+dpois(x=1,lambda=5)+dpois(x=2,lambda=5)+dpois(x=3,lambda=5)


#~~~ Digression ~~~#

# (Throughout these lab scripts, I'll include short digressions touching on additional 
# topics; read them if you're interested, but they aren't essential to the labs)

# When applied to a continuous distribution, this is the same as integrating the distribution 
# from -infinity to q, OR from q to +infinity (lower.tail=FALSE)
pnorm(q=-1,mean=0,sd=1,lower.tail=TRUE)
integrate(dnorm,lower=-100,upper=-1)

#~~~ End digression ~~~#


# This one can be a little confusing, but is essentially just the opposite of ppois() - 
# instead of supplying a critical value q and obtaining the probability of observing an 
# equally extreme value from the poisson distribution, this function allows us to supply 
# the probability and get an estimate of q, the critical value.
qpois(0.2650259,lambda=5,lower.tail=TRUE)


#### Probability distribution library: ####

# As mentioned before, this same suite of functions exists for pretty much every distribution 
# that R knows about. Off the top of my head, here's a list of the distributions most commonly
# used in R. See the EMD book for more - Chpt 4, specifically Table 4.1 and Figure 4.17

_unif			  # uniform distribution
_norm()			# normal distribution
_lnorm()		# log normal
_pois()			# poisson
_binom()		# binomial
_nbinom()		# negative binomial
_exp()			# exponential
_beta()			# beta
_betabinom()# beta binomial
_gamma()		# gamma
_chisq()		# chi-squared
_f()			  # F distribution
_t()			  # t distribution

# ... and a ton more (sometimes from specific R packages)

# For example, from the emdbook package:
#	dbetabinom
#	dmvnorm
#	dzinbinom

# You can explore the properties of some of these distributions interactively. Check out 
# 'distribution_demos.R' on the course site (must be run in R Studio).


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Generate random data from a poisson distribution ####

set.seed(49)	# this sets the seed of your computer's random number generator.

#~~~ Digression ~~~#

# (shhh! big secret: computers don't really generate 'new' random numbers when we ask for them. 
# instead they've got big huge lists of 'randomly' generated numbers produced by other processes
# that are saved and kept behind the scenes. When we ask for a random number, our computers go to 
# some spot on that list and start reading out numbers to us - sort of like randomly pointing 
# at a number in a phone book. If you tell your computer where to start reading each time, 
# you'll get the same results each time you run your code, and *hopefully* the same results 
# that I got when I wrote this code - potentially depending on your computer/operating system.)

#~~~ End digression ~~~#

ys<-rpois(5,lambda=7)
ys

# When I did this, I got:
# 6 7 5 8 3
# (did you get something different?)


####__Exercise__####

# Try generating 10 random numbers from a normal distribution


# What's the probability of observing a value as extreme (or more extreme)
# than 6.7 given a chisquared distribution with 3 degrees of freedom?



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Likelihood calculations ####

ys<-c(6,7,5,8,3)

# We just generated some random data from a poisson distribution, given a lambda of 7.
# Now let's go the other way; pretend that we don't know that lambda was 7. 
# Can we use these data to try to guess what value of lambda was used?


# Let's focus on the first value of our new, randomly generated data - in my case, 6

# This is our data, which we take to be true, and precisely known.  
# Now we want to ask the question "How likely is it that this observation 
# (or piece of data) came from a poisson distribution with a mean of lambda=3?"

# To calculate this, we can do the following:
lambda<-3
dpois(x=6,lambda=lambda)

# We find that the likelihood of lambda=3 being the mean of a poisson distribution that
# generated our observation (x=6) is 0.0504091

# What about for each of our other observations?
dpois(x=7,lambda=lambda)
dpois(x=5,lambda=lambda)
dpois(x=8,lambda=lambda)
dpois(x=3,lambda=lambda)

# These values tell us how likely it is that lambda equals 3 given each independent observation.

# How would we figure out how likely it is that lambda=3 given ALL of the data at once?
# - we can essentially treat likelihoods as probabilities.
#	- remember that the joint probability of independent events is equal to 
#   the product of the probabilities of each single event

# so the likelihood that lambda=3 given the first two observations would be:
dpois(x=6,lambda=lambda)*dpois(x=7,lambda=lambda)

# a shorter way of writing this would be to write:
prod(dpois(x=c(6,7),lambda=lambda))

# It follows that the likelihood that lambda=3 given all of our observations would be
prod(dpois(x=c(6,7,5,8,3),lambda=lambda))


# But why lambda=3? It was really just an arbitrary choice; or you could think of it 
# as a hypothesis. The key to this whole maximum likelihood approach is that we can 
# calculate what the likelihood is that lambda equals any given value (not just 3!), 
# based on the data that we observed. By calculating the likelihoods for a bunch of 
# different values of lambda, we can start to get a sense of which values of lambda 
# are most likely, given our data. In fact, we can calculate the value of lambda having 
# the largest (or maximum) likelihood give the data. This would be our 'maximum likelihood 
# estimate' of the parameter lambda.  Let's look at some examples.

# First let's look at the likelihood of different lambdas given just our first data point.
curve(dpois(6,lambda=x),0,20,xlab='lambda value',ylab='Likelihood')

# What value of lambda do you think is most likely?

# We could zoom in a bit on the range (5,7)
curve(dpois(6,lambda=x),5,7,xlab='lambda value',ylab='Likelihood')

# and again
curve(dpois(6,lambda=x),5.9,6.1,xlab='lambda value',ylab='Likelihood')

# Each time we could refine our visual guess as to the value of lambda corresponding to 
# the largest (or maximum) likelihood. But this is a pretty painful way of estimating the 
# maximum likelihood value of lambda; it's slow, and imprecise.

# Fortunately, there are other ways of going about this. In some cases we can derive a 
# mathematical representation or function that:
# 1) explicitly describes the likelihood of any given lambda based on our data, and 
# 2) can be readily differentiated and used to analytically locate its maximum value.

# <Recall the analytical result from lecture regarding poisson distribution and apply here>


# When our likelihood functions are messy (more complicated and/or not easily differentiable),
# it can be more practical (easier and faster) to use numerical techniques and computers to 
# come up with a reasonable estimate of the maximum of the likelihood function. 


# The first step is to write a likelihood calculator - this is a function that helps us 
# calculate the likelihood of different values of lambda, for specified data/observations.

# See the prod(dpois(...,lambda)) bit that we used earlier? we can generalize it within a 
# function, like this.

# poisson likelihood calculator:
poisLik<-function(data,lambda){
	sapply(lambda,function(x){ prod(dpois(data,lambda=x))})
}

# Now it will tell us the likelihood of an number of data points given an estimate of lambda 
# and a poisson distribution
poisLik(data=c(6),lambda=c(5))
poisLik(data=c(6,7),lambda=c(5))

# We can use a numerical 'optimization' algorithm to find the value of lambda that maximizes 
# the likelihood of our data. 

# check out
?optim

# We have to give optim() the function we want to optimize, a starting guess for the value of the 
# parameter we want to maximize with respect to, a method, and instructions on whether we 
# want it to find a maximum or a minimum.

# For example, we'll start with just the first data point:
optim(fn=poisLik,par=c(lambda=4),data=ys[1],method="BFGS",control=list(fnscale=-1))
# this yields an estimate of lambda=6, which makes a lot of sense, as our only data point was also 6

# We can try running this function using more of our data set, say 2 points:
optim(fn=poisLik,par=c(lambda=4),data=ys[1:2],method="BFGS",control=list(fnscale=-1))
# the new estimate for lambda is 6.5. Notice that it took a lot more steps (see $counts)
# to obtain this estimate of lambda than it did when we had just a single data point.

# Now lets try using all of the data:
optim(fn=poisLik,par=c(lambda=4),data=ys,method="BFGS",control=list(fnscale=-1))
# see the $convergence section of the output? when it's value is 1, it means optim failed to 
# find a maximum, which means we don't have a good estimate for lambda. Why did this happen?

# If we give optim a slightly better starting guess for lambda, and allow it to take a lot more
# steps before giving up on finding a maximum, we can still get an answer
optim(fn=poisLik,par=c(lambda=4.6),data=ys,method="BFGS",control=list(fnscale=-1,trace=T,maxit=10000))
# this is a real bummer, however, because we want to be able to analyze big data sets with lots of points.


### Never fear!

# When working with numerical optimizations there are a number of tricks that we can employ
# to make our task easier.

# First let's get a handle on why optim is having such a hard time. 
# An inherent property of likelihoods is that the more data you have to base your likelihood 
# estimates on, the smaller the joint likelihood (across all the data) becomes.  

# We can look at this graphically:

# Let's plot likelihood values across an whole range of possible lambda values, using 
# increasingly more data. As we do so, pay attention to the magnitude of the likelihood values...
curve(poisLik(data=ys[1],lambda=x),0,20,xlab='lambda value',ylab='Likelihood')
curve(poisLik(data=ys[1:2],lambda=x),0,20,xlab='lambda value',ylab='Likelihood',add=T)
curve(poisLik(data=ys[1:3],lambda=x),0,20,xlab='lambda value',ylab='Likelihood',add=T)
curve(poisLik(data=ys[1:4],lambda=x),0,20,xlab='lambda value',ylab='Likelihood',add=T)
curve(poisLik(data=ys[1:5],lambda=x),0,20,xlab='lambda value',ylab='Likelihood',add=T)

# By the time we're using 5 data points, even at its highest, the maximum of our likelihood 
# curve is no more than 4x10^-5 or 0.00004
par(mfrow=c(1,2))   # (splits our plotting window into two panels for comparing plots)
curve(poisLik(data=ys[1],lambda=x),0,20,xlab='lambda value',ylab='Likelihood',main="n = 1")
curve(poisLik(data=ys[1:5],lambda=x),0,20,xlab='lambda value',ylab='Likelihood',main='n = 5')

# These small values for the joint likelihood make it a lot harder for optim to find a maximum,
# because the difference in likelihood associated with lambda = 5 versus lambda = 15 is less
# than 1 x 10^-4 !  optim has to be very careful to detect such a small maximum.


#~~~ Digression: ~~~#

# It's also worth paying attention to the width of the likelihood curve - compare how steep 
# the function is when we use just a single data point, versus when we use all five. In general, 
# the more data we have, the steeper this likelihood curve or surface will get, and the 
# stronger/tighter our inference/conclusions will get - if there's any pattern to be
# found in our data. More on the strength of the likelihood surface in later labs.

#~~~ End digression ~~~#


### Fortunately, there are some tricks that we can use to get around this problem. 

# If, instead of finding the maximum of our likelihood function, we try to find the 
# minimum of the negative log likelihood ... (equal to the maximum of the log likelihood)

# Here's the negative log likelihood calculator for a poisson distribution:
poisNLL<-function(data,lambda){
	-sapply(lambda,function(x) sum(dpois(data,lambda=x,log=T)))
}
# notice we can tell dpois() to take the log of likelihoods with log=T, and here we are 
# finding the sum of the log likelihoods as supposed to the product of the likelihoods

# We can run optim on this poisNLL instead of the poisLik, and it works well for the whole data set
optim(fn=poisNLL,par=c(lambda=4.6),data=ys,method="BFGS",control=list(trace=T,maxit=5000))
# converges rapidly on essentially the same value for lambda (5 steps rather than >5,000 steps).


# Why is this ok to do?

#	(1) taking the log of the likelihoods changes their absolute magnitudes, but not the location 
# of their peaks and valleys.  In otherwords, the maximum likelihood is still the maximum with
# respect to the other likelihood values after you log transform them all.

#	(2) the log of a bunch of things multiplied together is the same as the sum of the log of each
# of the individual terms.

# For example
a<-log(dpois(x=6,lambda=lambda)*dpois(x=7,lambda=lambda))
b<-log(dpois(x=6,lambda=lambda))+log(dpois(x=7,lambda=lambda))
c(a,b)
a==b

#	(3) Together, these properties turn out to be a big advantage computationally, as numbers that
# are really close to zero are difficult for computers to keep track of. Taking the log makes 
# numbers 'larger', as does being able to work with sums rather than the product of decimal 
# numbers typically less than 1 (always a smaller decimal number).

#	(4) Finally, there is an (unfortunate) historical artifact: early numerical algorithms were 
# designed for finding the minimum of a function, and so the default behavior of most 
# numerical optimization functions is to find a minimum. We want to find the maximum of our 
# likelihood functions though.  To accomplish this, we can multiply our log likelihood 
# function by -1, such that all of our maxima now represent minima (found at exactly the same
# locations). We can then ask our numerical algorithms to find the minimum of our negative 
# log likelihood function. This corresponds to finding the maximum of the log likelihood 
# function, which, in turn, is the same as finding the maximum of the likelihood function.

# Let's get a sense of this graphically:

# poisson log likelihood function:
poisLL<-function(data,lambda){
	sapply(lambda,function(x) sum(dpois(data,lambda=x,log=T)))
}

par(mfrow=c(3,1))
curve(poisLik(ys[1:5],x),2,12,xlab='lambda',ylab='Likelihood')
segments(5.8,-200,5.8,200,col='red')
curve(poisLL(ys[1:5],x),2,12,xlab='lambda',ylab='Log Likelihood')
segments(5.8,-200,5.8,200,col='red')
curve(poisNLL(ys[1:5],x),2,12,xlab='lambda',ylab='Negative Log Likelihood')
segments(5.8,-200,5.8,200,col='red')




# Let's stop and think for a moment - where are we going with all of this? 
#	(1) Why are we doing this? 
#	(2) What do the results of this optimization tell us?  
#	(3) What is significant about lambda = 5.8?
#		  - is this a surprising answer?



### Further exploration:

# What happens if we increase sample size?
set.seed(33)

# The following code creates mock data sets of increasingly large sizes, then estimates lambda using likelihood. You can look at the details of the code, or just run it to see the graph at the end.

par(mfrow=c(3,3))						          # set up graphing window
ns<-c(10,20,50,100,200,500,1000,2000)	# vector of sample sizes to try
ls<-rep(NA,length(ns))					      # empty vector for storing lambda estimates
for(i in 1:length(ns)){					      # loop through each sample size
	ys2<-rpois(ns[i],lambda=7)			    # random data from poisson distrib. with sample size ns[i]
	res<-optim(fn=poisNLL,par=c(lambda=20),data=ys2,    # fit a poisson distribution to the data
	           method="BFGS",control=list(maxit=10000))			
	ls[i]<-res$par						          # save the resulting lambda estimate
	hist(ys2,30,col='blue',main=paste('n =',ns[i],", lambda =",round(res$par,2)),
	     xlab='data',xlim=c(0,19))	
}
# plot estimates as a function of sample size
plot(ls~ns,type='l',xlab='Sample size',ylab='lambda estimate')	
points(ls~ns,pch=19,cex=0.7)	# annotate
segments(0,7,2000,7,lty=3)


# What do these graphs mean? 
# How are our estimates of lambda changing, and why? 
# What are the downsides of large sample sizes?




#### A few comments on methods: the bbmle package ####

# We can take a peak at what's happening inside of the optim() function...
# example of tracing out the steps taken en route to convergence on an estimate of
# the minimum negative log likelihood:
lambda.guess<-20
nits<-optim(fn=poisNLL,par=c(lambda=lambda.guess),data=ys,method="BFGS",control=list(trace=T,maxit=100,REPORT=1))$counts
res<-c(lambda.guess,sapply(seq(1,nits[[1]],1),function(x) optim(fn=poisNLL,par=c(lambda=lambda.guess),data=ys,method="BFGS",control=list(maxit=x))$par))
par(mfrow=c(1,1))
plot(res,type='l',xlab='Steps',ylab='lambda estimate',main=paste("lambda guess = ",lambda.guess))


# Up until now, we've been using the function optim() to perform our numerical optimization 
# (minimization/maximization). There is another collection of functions that we'll start 
# using extenstively that customize the optim() function to make it particularly easy to 
# use for the purposes of maximum likelihood analysis. These are kept in an R package 
# called 'bbmle' (Ben Bolker's mle package).

library(bbmle)	# this loads the bbmle package

# If you don't already have it installed, run this:
# install.packages("bbmle")

# We can use the central function of the bbmle package, mle2, to do the same calculations 
# we just investigated above. For example (using a negative log likelihood calculator 
# slightly modified to handle the data differently):

poisNLL2<-function(lambda){
	-sum(dpois(ys,lambda=lambda,log=T))
}

fit1<-mle2(poisNLL2,start=list(lambda=lambda.guess))
summary(fit1)


#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Numerical Optimization: the good, the bad, and the ugly... ####

# We've just investigated how numerical optimization enables us to locate the maximum 
# likelihood estimates of parameters, given data. The examples that we've looked at 
# consist of small amounts of quite tidy data, and the model that we've been fitting to 
# the data is as simple as it gets - one parameter. We can even find the answers to this 
# optimization problem analytically, without resorting to numerics. 

# However, as problems get more complicated - either through the addition of more 
# (and noisier) data and covariates, or the use of more complicated statistical models, 
# the number of parameters, and the difficulty of estimating their values, increases. This 
# quickly takes us into regions where we can no longer analytically arrive at maximum 
# likelihood estimates. Numerical methods are awesome because they let us arrive at 
# solutions to more difficult problems. 

# But they come with their own limitations:
#
#	- the complexity of likelihood surfaces grows quickly with increasing number of parameters 
#   (each parameter adds another dimension to the problem)
#
#	- the time and computational resources required for locating a minima/maxima grows with the
#   addition of more parameters. This limitation is what originally prevented statisticians 
#   from using maximum likelihood techniques. As computing power grows and becomes less 
#   expensive all the time, this is much much less of a limit than it used to be, but not 
#   inconsequential for some problems.
#
#	- numerical methods can be fooled into locating local minima/maxima, rather than the global
#   minima/maxima that we'd like to locate. 

#	To explore how some common numerical optimization algorithms work (with accompanying 
# graphics), check out the optimization_demos.cdf example document. It also provides an 
# example of how numerical optimization algorithms can be fooled by local minima when 
# likelihood curves/surfaces are complex (aka 'woogedy').

# Numerous algorithms exist, and often when you have a complicated problem, it can be 
# advantageous to try out different algorithms, to see if they give consistent answers, 
# or if one of them runs much more quickly or stably than the others.  Here's a short list 
# of some available options that we can use for MLE problems (for more detail, see EMD, 
# chapter 7):
?optim


# "Nelder-Mead"

# Classic, work horse of minimization algorithms. From 1965; robust, but sometimes slow. 
# This is the default algorithm for optim(), and is the algorithm demonstrated in the interactive examples. 

# Wikipedia suggests that this algorithm is also called the 'amoeba method' - how can you go 
# wrong there?
# It's not derivative based, which improves its performance on bumpy likelihood surfaces. 
# Ben Bolker warns that it's unreliable for optimization in 1 dimension (ie, single parameters), 
# and in contrast to optim(), the default method for mle2 is "BFGS"
# http://en.wikipedia.org/wiki/Nelder–Mead_method

# "BFGS"

# a 'quasi-Newton' method, published by Broyden, Fletcher, Goldfarb and Shanno (hence BFGS - 
# intuitive name to remember, huh? I remember it b/c it reminds me of a Roald Dahl book, 
# 'The BFG', standing for Big Friendly Giant. Hey, whatever it takes...)
# This is the default method for mle2, and relies on being able to estimate the derivatives 
# of the negative log likelihood function, making it susceptible to problems when the 
# likelihood surface is bumpy.


# "L-BFGS-B"

# A modification of the "BFGS" method that allows you to specify bounds or constraints on 
# the value of parameters. This can be a really useful feature, because often parameters of 
# distributions are necessarily constrained (for example, in the poisson distribution, 
# lambda must be > 0). This method suffers from the same issues as the BFGS, and, 
# additionally, starting guesses for parameters must satisfy the imposed constraints.

# "SANN" = Simulated Annealing

# based on concepts drawn from physics/thermodynamics - heating up metal and letting it cool 
# very slowly... Belongs to a class of algorithms know as stochastic global optimization 
# methods. Can be quite slow, but uses only values taken directly from the function to be 
# minimized (so there's no requirement that the function be differentiable, or need to 
# estimate the derivative of a function).


# "CG"

# or 'conjugate gradients' method. Haven't used this one much; help file suggests 
# "Conjugate gradient methods will generally be more fragile than the BFGS method, but as 
# they do not store a matrix they may be successful in much larger optimization problems."


# "Brent"

# For one dimensional optimization only. Derivative free method that attempts to bracket a 
# minima. Fast and reliable for simple problems.


# Compare various methods, based on poisson example above.

f1<-mle2(poisNLL2,start=list(lambda=lambda.guess),method="BFGS")

f2<-mle2(poisNLL2,start=list(lambda=lambda.guess),method="L-BFGS-B",lower=c(lambda=0.0001))

f3<-mle2(poisNLL2,start=list(lambda=lambda.guess),method="Nelder-Mead")

f4<-mle2(poisNLL2,start=list(lambda=lambda.guess),method="SANN")

f5<-mle2(poisNLL2,start=list(lambda=lambda.guess),method="CG")

f6<-mle2(poisNLL2,start=list(lambda=lambda.guess),method="Brent",lower=c(1),upper=c(10))

coef(f1)
coef(f2)
coef(f3)
coef(f4)
coef(f5)
coef(f6)

# We can also check out slow each method is, using the function:
# ?system.time()

methods<-c("BFGS","Nelder-Mead","SANN","CG")

# Apply this function to calculate the mean response times for different algorithms:
sapply(methods, function(x) mean(sapply(seq(1,10,1), function(y) system.time(mle2(poisNLL2,start=list(lambda=lambda.guess),method=x))[[1]])))



#----------------------------------------------------------------#
#----------------------------------------------------------------#

#### Recap: ####

# In this lab script, we've covered:

#	- the basics of working with probability distributions in R
#	- how to generate random data (from a poisson distribution)
#	- how to use maximum likelihood methods to estimate the 
#	  parameters generating that random data.
#		- the relationship between likelihood, log likelihood, 
#		  and negative log likelihood functions
#	- an introduction to the bbmle package/mle2() function
#	- (dis-)advantages of numerical optimization
#	- and an introduction to various numerical optimization algorithms

#----------------------------------------------------------------#
#----------------------------------------------------------------#
	

#### 'Answers' ####

pchisq(6.7,df=3, lower.tail = F)

	