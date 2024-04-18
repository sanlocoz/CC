
##################
# GETTING STARTED 
##################

# this statement removes everything from R's memory
rm(list=ls())

# change working directory
setwd("D:/Wageningen/Period 5/Catchment/Module 2/Practical/recession_based_models")

# load functions (for later)
source("function_select_part.R")
source("function_plot_output.R")

# Load package for calibration.
# You may still need to install the package.
library(minpack.lm)


###
### READ DATA
###

d = read.table("PEQ_Hupsel_hour.dat", header=TRUE)





###########################
# 6.2 DIFFERENTIAL EQUATION
###########################

# Write a function with name model_equation and arguments xi, Pi, ETi, a and b.
# The function should return the right-hand side of the differential equation.

diff_eq = function(xi, Pi, ETi, a, b)
{
  #cat(xi, Pi, ETi, a, b )
  #print("\n")
  return(a * (exp(xi)^(b-2) * (Pi - ETi)) - a* exp(xi)^(b-1))
}






####################################
# 6.3 RUNGE-KUTTA INTEGRATION SCHEME
####################################

# Write a function with name runge_kutta and arguments xi, Pi, ETi, a and b.
# The function should return a new value of x.

runge_kutta = function(xi, Pi, ETi, a, b)
{
  k1 = diff_eq(xi     , Pi, ETi, a, b)
  k2 = diff_eq(xi+k1/2, Pi, ETi, a, b)
  k3 = diff_eq(xi+k2/2, Pi, ETi, a, b)
  k4 = diff_eq(xi+k3  , Pi, ETi, a, b)
  # cat(k1,k2,k3,k4)
  # print(k2)
  return(xi+(k1+2*k2+2*k3+k4)/6)
}




##############
# 6.4 FOR-LOOP
##############

# Write a function with name run_model and arguments forc, a, b and Q0.
# The function should:
# - Make an empty vector for x and fill the first element using the observed discharge at t=0. 
# - Run a for-loop over all time steps (starting at 2).
# - Return a vector with computed discharges.

run_model = function(forc, a, b, Q0)
{
  # Make an empty vector for x and fill the first element using the observed discharge at t=0. 
  x    = c()
  x[1] = log(Q0)
  
  # Run a for-loop over all time steps (starting at 2).
  for(i in 2:length(forc$Q))
  {
    x[i] = max(-5, runge_kutta(x[i-1], forc$P[i], forc$ET[i], a=a, b=b))
    # print(x[i])
  }
  
  # Return a vector with computed discharges.
  return(exp(x))
}







#####
# RUN
#####

# Call function select_part to cut out part of the data frame.
forc = select_part(dataframe=d, start=2012000000, end=2015000000)

# Run the model.
Qmod = run_model(forc, a=0.05, b=2, Q0=forc$Q[1])

# Give the run a logical name.
run_name = "test"

# Plot results (read and run the function below once first).
#plot_output(Qmod, run_name)

# Write the output to file.
write.table(data.frame(forc,Qmod), paste("output/", run_name, ".dat", sep=""), row.names=FALSE)


####################
# HYDROGRAPH FITTING
####################

# Define initial and boundary values for the parameters you want to calibrate.
# You can also calibrate the initial values (if you want).
# a  = runif(500, min=0.01, max=1)
# b  = runif(500, min=0.5, max=4)
# eff  = runif(500, min=0  , max=1)
# 
# KGE = c()
# 
# library(hydroGOF)
# 
# for(i in 1:5000)
# {
#   print(i)
#   forc_copy = forc
#   forc_copy$ET = forc_copy$ET*eff[i]
#   #print(forc_copy$ET)
#   Qmod = run_model(forc_copy, a=a[i], b=b[i], Q0=forc_copy$Q[1])
#   KGE[i]= KGE(sim=forc_copy$Q, obs= Qmod, na.rm=T)
# }
# 
# pars = cbind(a,b,eff,KGE)

#write.table(pars, "output/pars_Hupsel_recess.dat", row.names=FALSE)

parameter = read.table("output/pars_Hupsel_recess.dat", header=TRUE)

best_pars = parameter[order(-parameter$KGE)[1:10],]

Qmod = matrix(ncol=10, nrow=nrow(forc))
storage = matrix(ncol=10, nrow=nrow(forc))
for(i in 1:10)
{
  forc_copy = forc
  forc_copy$ET = forc_copy$ET*best_pars$eff[i]
  #print(forc_copy$ET)
  mod = run_model(forc_copy, a=best_pars$a[i], b=best_pars$b[i], Q0=forc_copy$Q[1])
  Qmod[,i] = mod
  storage[,i] = 1/best_pars$a[i] * 1/(2-best_pars$b[i]) * mod^(2-best_pars$b[i])
  #plot_output(mod, "hde")
}

write.table(Qmod, "output/dischargeR.dat")
write.table(storage, "output/storageR.dat")
write.table(forc_copy$date, "output/date.dat")

