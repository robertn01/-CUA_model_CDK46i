# Replicating cost-utility analysis model for Palbociclib based on Excel file
# Used the Markov model frame from the UCL training
# 3-state deterministic Markov-model - also enables probabilistic implementation
# Created by Robert Nagy (PhD student)
# 04/08/2019
# v1.2
# Notes: sumbitted model did not applied half-cycle correction
###############################################################

# cleaning the working environment

rm(list = ls())
graphics.off() # remove graphics, plots from memory

# set working directory - please use the preferred file path

setwd("/Users/robertn/Desktop/R_projects/PhD_Markov_model/Markov_model")  # on Linux\ Unix, Mac use '\' rather than '//'
getwd() # check

# Load necessary libraries
# If not installed use the following line first
# install.packages("VGAM")

library(VGAM)

# Set an arbitrary seed for random draf from distribution in the global environment (global variable)

set.seed(55)

# Define the number and names of treatments
# These are Standard of Care with website
# and Standard of Care without website

n.treatments<-2
treatment.names<-c("Palbociclib + Letrozole","Placebo + Letrozole")

# Define the number and names of states of the model
# This is two and they are "Smoking" and "Not smoking"
n.states<-3
state.names<-c("Pre-progressed","Progressed","Dead")

# Define the number of cycles
# This is 480: time horizon is 40 years and cycle length is 1 month
# The code will work for any even n.cycles (need to change the discounting code if
# an odd number of cycles is desired)

n.cycles<-480

# Define simulation parameters
# This is the number of PSA samples to use
n.samples<-1000

# WTP threshold (in £/QALY)
lambda = 30000

####### Define the input parameters of the model #########################################

# The transition matrix is a 3x3 matrix
# Rows sum to 1

# There is one transition matrix for each treatment option and each PSA sample
# Store them in an array with (before filling in below) NA entries: object constructor


transition.matrices<-array(0,dim=c(n.treatments,n.samples,n.states,n.states),
                           dimnames=list(treatment.names,NULL,state.names,state.names))

# First the transition matrix for standard care: Placebo + Letrozole (PFS exp rate= 0.0475; OS exp rate = 0.0165)

# Transitions from Pre-progressed state (PSA - Exponential)
transition.matrices["Placebo + Letrozole",,"Pre-progressed","Progressed"]<- 1-exp(-log(2)/14.6)  #rexp(n.samples,rate = 1/0.0475)#1-exp(-0.0475)
transition.matrices["Placebo + Letrozole",,"Pre-progressed","Dead"]<-1-exp(-log(2)/42)   #rexp(n.samples, rate = 1/0.0165)
transition.matrices["Placebo + Letrozole",,"Pre-progressed","Pre-progressed"]<- 1-(transition.matrices[2,,1,2] + transition.matrices[2,,1,3]) #1-sum(transition.matrices[2,,1,2:3])


#colMeans(transition.matrices[2,,1,])

# Transitions from Progressed state (PSA - Exponential)
transition.matrices["Placebo + Letrozole",,"Progressed","Dead"]<-1-exp(-log(2)/42)   # rexp(n.samples, rate = 1/0.0165)
transition.matrices["Placebo + Letrozole",,"Progressed","Progressed"]<-1-transition.matrices[2,,2,3]

# Transitions from Dead state: same for both therapeutic arm, therefore:
transition.matrices["Placebo + Letrozole",,"Dead","Dead"]<-transition.matrices["Palbociclib + Letrozole",,"Dead","Dead"]<-1

# The transition matrix for new treatment: Palbociclib + Letrozole (PFS exp rate= 0.0274; OS exp rate = 0.0141)

# Transitions from Pre-progressed state (PSA - Exponential)
transition.matrices["Palbociclib + Letrozole",,"Pre-progressed","Progressed"]<-1-exp(-log(2)/25.3) # rexp(n.samples, rate = 1/0.0274)
transition.matrices["Palbociclib + Letrozole",,"Pre-progressed","Dead"]<-1-exp(-log(2)/49.1)   # rexp(n.samples, rate = 1/0.0141) 
transition.matrices["Palbociclib + Letrozole",,"Pre-progressed","Pre-progressed"]<-1-(transition.matrices[1,,1,2] + transition.matrices[1,,1,3])

# Transitions from Progressed state (PSA - Exponential)
transition.matrices["Palbociclib + Letrozole",,"Progressed","Dead"]<-1-exp(-log(2)/49.1)   # rexp(n.samples, rate = 1/0.0141)
transition.matrices["Palbociclib + Letrozole",,"Progressed","Progressed"]<-1-transition.matrices[1,,2,3]


# Now define the QALYS associated with the states per cycle

# There is one for each PSA sample and each state
# Store in an NA array and then fill in below

state.qalys<-array(dim=c(n.samples, n.states),dimnames=list(NULL,state.names))

#  as no information about the sd - both costs and QALYs remain fixed for all the PSA samples

state.qalys[,"Pre-progressed"]<- 0.721/12# rnorm(n.samples,mean=0.06,sd=?)/12
state.qalys[,"Progressed"]<- 0.5052/12
state.qalys[,"Dead"]<- 0


# And finally define the state costs per cycle

state.costs<-array(0,dim=c(n.samples, n.states),dimnames=list(NULL,state.names))

state.costs[,"Pre-progressed"]<- 5 #pre-progressed refers to line 1 treatment: non-drug costs only! I've used arbitrary £5
state.costs[,"Progressed"]<-1200  # here refers to line 2, 3, 4 and BSC (disease-related costs & subsequent treatment-related costs combined average)
state.costs[,"Dead"]<- 0


# Define the treatment costs per cycle

# One for each PSA sample and each treatment
# Treatment costs are actually fixed but this allows flexibility if we
# want to include uncertainty/randomness in the cost

treatment.costs<-array(dim=c(n.treatments,n.samples),dimnames=list(treatment.names,NULL))

# VAT incl/excl?
treatment.costs["Palbociclib + Letrozole",]<-2954
treatment.costs["Placebo + Letrozole",]<-4


###### Simulation ###########################################################

# Build an array to store the cohort vector at each cycle

# Each cohort vector has 3 (=n.states) elements: probabilities (!) being in the given state
# There is one cohort vector for each treatment, for each PSA sample, for each cycle.

cohort.vectors<-array(dim=c(n.treatments,n.samples,n.cycles,n.states),
                      dimnames=list(treatment.names,NULL,NULL,state.names))

# Assume that everyone starts in the progression-free state
cohort.vectors[,,1,"Pre-progressed"]<-1
cohort.vectors[,,1,"Progressed"]<-0
cohort.vectors[,,1,"Dead"]<-0

# Build an array to store the costs and QALYs accrued per cycle

# One for each treatment, for each PSA sample, for each cycle
# These will be filled in below in the main model code
# Then discounted and summed to contribute to total costs and total QALYs

cycle.costs<-array(dim=c(n.treatments,n.samples,n.cycles),
                   dimnames=list(treatment.names,NULL,NULL))

cycle.qalys<-array(dim=c(n.treatments,n.samples,n.cycles),
                   dimnames=list(treatment.names,NULL,NULL))

# Build arrays to store the total costs and total QALYs

# There is one for each treatment and each PSA sample
# These are filled in below using cycle.costs, 
# treatment.costs, and cycle.qalys

total.costs<-array(dim=c(n.treatments,n.samples),
                   dimnames=list(treatment.names,NULL))

total.qalys<-array(dim=c(n.treatments,n.samples),
                   dimnames=list(treatment.names,NULL))


# The remainder of the cohort.vectors will be filled in by Markov updating below

######### Main model code #################################################################

# Loop over the treatment options
for(i.treatment in 1:n.treatments)
{
  # Loop over the PSA samples
  for(i.sample in 1:n.samples)
  {
    # Loop over the cycles
    # Cycle 1 is already defined so only need to update cycles 2:n.cycles
    for(i.cycle in 2:n.cycles)
    {
      # Markov update (Markov trace?)
      # Multiply previous cycle's cohort vector by transition matrix
      # i.e. pi_j = pi_(j-1)*P
      cohort.vectors[i.treatment,i.sample,i.cycle,]<-
        cohort.vectors[i.treatment,i.sample,i.cycle-1,] %*%
        transition.matrices[i.treatment,i.sample,,]
    }
    
    # Now use the cohort vectors to calculate the 
    # total costs for each cycle
    
    cycle.costs[i.treatment,i.sample,]<-
      cohort.vectors[i.treatment,i.sample,,] %*% state.costs[i.sample,] + cohort.vectors[i.treatment,i.sample,,1] * treatment.costs[i.treatment,i.sample]
    
    # And total QALYs for each cycle
    cycle.qalys[i.treatment,i.sample,]<-
      cohort.vectors[i.treatment,i.sample,,] %*% state.qalys[i.sample,]
    
    # Combine the cycle.costs and treatment.costs to get total costs
    
    # Apply the discount factor - now assume 3.5% annual discount
    
    # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
    # Each year acounts for twelve cycles so need to repeat the discount values
    
    disc <- (1/1.035)^rep(c(0:(n.cycles/12-1)),each=12)
    
    total.costs[i.treatment,i.sample]<-cycle.costs[i.treatment,i.sample,]%*%disc
    
      # (1/1.035)^rep(c(0:(n.cycles/12-1)),each=12)
    
    # Combine the cycle.qalys to get total qalys
    
    # Apply the discount factor 
    
    # (1 in first year, 1.035 in second, 1.035^2 in third, and so on)
    # Each year acounts for 12 cycles so need to repeat the discount values
    total.qalys[i.treatment,i.sample]<-cycle.qalys[i.treatment,i.sample,]%*%disc
      # (1/1.035)^rep(c(0:(n.cycles/12-1)),each=12)
  }
}

##### Analysis of results ######################################################


# Average costs

show(average.costs<-rowMeans(total.costs)) # not a big difference in costs

# Average effects (in QALY units)

show(average.effects<-rowMeans(total.qalys))

# Incremental costs and effects relative to 'standard of care' (here: letrozole alone)

incremental.costs<-total.costs["Palbociclib + Letrozole",]-total.costs["Placebo + Letrozole",]

incremental.effects<-total.qalys["Palbociclib + Letrozole",]-total.qalys["Placebo + Letrozole",]

# The ICER 

ICER<-mean(incremental.costs)/mean(incremental.effects)

# Incremental net benefit at the willingness-to-pay threshold, lambda
# Sometimes positive and sometimes negative 
# Need to look at averages and consider probabilities of cost-effectiveness
incremental.net.benefit<-lambda*incremental.effects-incremental.costs

# Average incremental net benefit
# This is positive indicating cost-effectiveness at the lambda threshold
average.inb<-mean(incremental.net.benefit)

# Probability cost-effective
# This is the proportion of samples for which the incremental net benefit is positive

probability.cost.effective<-sum(incremental.net.benefit>0)/n.samples

length(which(incremental.net.benefit>0))/n.samples # same thing

# Now use the BCEA package to analyse the results...


##################### Analysing the results/ BCEA - PRACTICAL #############################

################## 1.Understanding the Markov update loop ##############################

# a)

dim(transition.matrices) #actually this is a 4 dim space: dimensionality VS product of nr of values allowed per dimension..

transition.matrices["Placebo + Letrozole",1,,]

transition.matrices["Palbociclib + Letrozole",1,,]

# b)

transition.matrices["Palbociclib + Letrozole",1,"Pre-progressed",]

colMeans(transition.matrices["Palbociclib + Letrozole",,"Pre-progressed",])

colMeans(transition.matrices["Placebo + Letrozole",,"Pre-progressed",])


colMeans(transition.matrices["Palbociclib + Letrozole",,"Progressed",])

colMeans(transition.matrices["Placebo + Letrozole",,"Progressed",])

# c) Markov update

dim(cohort.vectors) # indicated 'pi' in lecture notes

cohort.vectors["Palbociclib + Letrozole",1,1,] #initially everyone starts at the 'Pre-progressed' state
transition.matrices["Palbociclib + Letrozole",1,,]

# Markov update formula: Pi(t) = Pi(t-1)*P, where P denotes the appropriate transition matrix

# Now, we can perform matrix multiplication to get the cohort vector for the second cycle on 
# "Palbociclib + Letrozole" in the first PSA sample.

cohort.vectors["Palbociclib + Letrozole",1,1,]%*%transition.matrices["Palbociclib + Letrozole",1,,]
cohort.vectors["Palbociclib + Letrozole",1,2,] #2nd cycle for 1st PSA

# Repeating this for each treatment (i.treatment), each PSA sample (i.sample), and for 
# each cycle (i.cycle) we get the line of code within the three nested for loops

cohort.vectors[i.treatment,i.sample,i.cycle,]<-
  cohort.vectors[i.treatment,i.sample,i.cycle-1,]%*%
  transition.matrices[i.treatment,i.sample,,]


############## 2. Using BCEA to analyse the results ##################################################

library(BCEA)

# a)
mm<-bcea(e = t(total.qalys), c = t(total.costs), ref = 2, interventions = treatment.names, Kmax = 100000) #probably ref should be 1, Kmax also modified

# b)

# The "EIB" is expected incremental
# benefit at the lambda wtp, the "CEAC" (cost-effectiveness acceptability curve) is the probability 
# that the reference of "no treatment" has highest net benefit (most cost-effective) 
# at the specified willingness-to-pay, and the ICER is the incremental cost-effectiveness ratio. 
# The last of these can be compared with the standard willingness-to-pay threshold of £X.

summary(mm) #Cost-effectiveness analysis summary

# c)
ceplane.plot(he = mm, wtp = lambda)
rowMeans(total.costs)
mean(total.costs["Palbociclib + Letrozole",]) - mean(total.costs["Placebo + Letrozole",])

rowMeans(total.qalys)
mean(total.qalys["Palbociclib + Letrozole",]) - mean(total.qalys["Placebo + Letrozole",])

# d)

ceac.plot(he = mm)
# We see that this converges to around 75% as the willingness-to-pay increases

######### 3. Using the info.rank() function from BCEA ###############################

# The info.rank() function will be used to estimate the proportion of the total 
# decision uncertainty, quantified by expected value of perfect information (EVPI) , 
# to which each of the uncertain parameters contribute.

# [TBC'd...]
