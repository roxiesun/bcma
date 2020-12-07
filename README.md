# BCMA: Bayesian Causal Mediation Analysis with Latent Mediators and Survival Outcome

This page illustartes a simulation on Bayesian causal mediation analysis with latent mediators and a survival outcome. The methodology and real data application to [ADNI](http://adni.loni.usc.edu/) data can be found in the original paper:
[#](https://).


The proposed BCMA procedure

## Usage

R version 2.8.1 or after. Package 'mvtnorm' and 'matrix sampling' are used.

#### List of R Function


main.sim.R		--- Main program for simulation.

def_con.R		  --- Define constant.

def_rec.R		  --- Allocation of computation variables.

true.par.R		--- Set true values of parameters

iden.R			--- Define value concerning model identification.

init1.R			--- Set initial value for parameters.

initial.R		--- Initialise the latent variable and its indicators.

prior.R			--- Set hyperparameters of prior distribution.

sim_cesore.R		--- Simulate data.

readdata.R		--- Read data.

MCMC.R			--- Run MCMC algorithm.

function.R		--- Define some functions used in MCMC.R.


#### Implementation

step 1: In def_con.R, change constant according to model.

step 2: In true.par.R, 	change true parameter values.

step 3: In iden.R, 	change identification variable according to model specification.

step 4: In init1.R,	change initial value of parameters.

step 5: In prior.R, 	change hyperparameters according to prior input.

step 6: Run the main.sim.R.


