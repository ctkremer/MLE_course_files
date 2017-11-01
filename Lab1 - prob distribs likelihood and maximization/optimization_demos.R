
#----------------------------------------------------------------#
#----------------------------------------------------------------#

# Interactive optimization demos
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

#### Demonstrations of optimization algorithms ####

# Optimization for the Poisson distribution
#   - use 'guess' to manipulate the starting guess for lambda
#   - use 'nsteps' to control how many steps forward in time the optimizer halts
manipulate(nmplot1(guess,nsteps),guess=slider(0,12,initial=2),nsteps=slider(1,30,initial=2))

# Non-local optima demonstration
#   - change the initial guess to see if you can trap the optimization 
#     algorithm in different valleys
manipulate(nmplot2(guess,nsteps),guess=slider(-4,10,initial=2),nsteps=slider(1,50,initial=2))

# optimization in 2 dimensions via Nelder-Mead
#   - try changing the starting guesses, and manipulating nsteps
manipulate(nmplot3(mu.guess,s.guess,nsteps),mu.guess=slider(0,15,initial=2),s.guess=slider(0.1,8,initial=6),nsteps=slider(3,50,initial=3))
