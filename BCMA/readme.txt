R version 2.8.1 or after. Package 'mvtnorm' and 'matrix sampling' are used.
*************************
LIST OF R FUNCTION
******************

main.sim.R		---- Main program for simulation.
def_con.R		---- Define constant.
def_rec.R		---- Allocation of computation variables.
true.par.R		---- Set true values of parameters
iden.R			---- Define value concerning model identification.
init1.R			---- Set initial value for parameters.
initial.R		---- Initialise the latent variable and its indicators .
prior.R			---- Set hyperparameters of prior distribution.
sim_cesore.R		---- Simulate data.
readdata.R		---- Read data.
MCMC.R			---- Run MCMC algorithm.
function.R		---- Define some functions used in MCMC.R.

*****
USAGE
*****

step 1: In def_con.R, 	change constant according to model.
step 2: In true.par.R, 	change true parameter values.
step 3: In iden.R, 	change identification variable according to model specification.
step 4: In init1.R,	change initial value of parameters.
step 5: In prior.R, 	change hyperparameters according to prior input.
step 6: Run the main.sim.R.


*******
Results
*******

The acceptance rates value will be printed on R screen.

Estimation of parameters are saved in lists and written in csv files with corresponding names in the working directory.
