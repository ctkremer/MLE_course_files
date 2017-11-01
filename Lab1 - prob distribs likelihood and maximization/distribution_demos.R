
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Interactive probability distributions
#
#	  - Must be run in R Studio
# 	- Developed by Colin T. Kremer for ELME 2014
#   - Checked 2015, 2017

#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Load library
library(manipulate)
library(dplyr)
library(ggplot2)

# load required tools
source("/Users/colin/Teaching/ELME/ELME 2017/labs/mle.tools_070717.R")

#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Examples:

# Poisson distribution
manipulate(poisplot(lambda),lambda=slider(0.01,25))

# Binomial distribution
manipulate(binomplot(size,prob),size=slider(0,50,initial=10),prob=slider(0,1,initial=0.2))

# Negative binomial distribution
manipulate(nbinomplot(mu,size),size=slider(0,50,initial=50),mu=slider(0,30,initial=20))

# Normal distribution
manipulate(normplot(mean,sd),mean=slider(-30,30,initial=0),sd=slider(0.01,5,initial=2))

# Compare how well the normal distribution approximates the poisson
manipulate(pois.normplot(lambda),lambda=slider(0.01,80))

# Beta distribution
manipulate(betaplot(mu,phi),mu=slider(0,1,initial=0.5,step=0.05),phi=slider(0,5,initial=2,step=0.25))

# Gamma distribution
manipulate(gammaplot(shape,scale),shape=slider(0,8,initial=2,step=0.25),scale=slider(0,10,initial=2,step=0.25))

# Exponential distribution
manipulate(expplot(rate),rate=slider(0.05,2,initial=2,step=0.05))




