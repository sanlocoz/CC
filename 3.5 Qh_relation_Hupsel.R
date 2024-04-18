
#################
# Getting started 
#################

# Remove everything from R's memory.
rm(list=ls())

# Load the WALRUS package.
library(WALRUS)

# Change working directory to the folder where data-, figures- and output-subfolders 
# are located.
setwd("D:/Wageningen/Period 5/Catchment/Module 3/Practical/WALRUS")


######
# Data
######

# Read daily or hourly precipitation, potential evapotranspiration and discharge data.
data = read.table("data/PEQ_Hupsel_hour.dat", header=TRUE)
data_ET = read.table("data/calibrationet.csv", header=TRUE, sep = ";")
data_ET$date=paste(substr(data_ET$date,7,10),substr(data_ET$date,4,5),substr(data_ET$date,1,2), sep="")
data_ET$date = as.numeric(strptime(as.character(data_ET$date), 
                           format = "%Y%m%d", tz = "UTC"))

# Specify which period of the total data set you want to use as forcing data.
# Use the same date format as in data file (for example yyyymmddhh).
forc = WALRUS_selectdates("data", 2011060000, 2011090000)
forc = WALRUS_selectdates("data", 2012060000, 2012090000)
forc = WALRUS_selectdates("data", 2012000000, 2014043100)



# Preprocessing forcing and specifying for which moments output should be stored. 
# The argument dt is time step size used in the output data file (in hours).
WALRUS_preprocessing(f=forc, dt=1)
WALRUS_preprocessing_calibration()

start_date = forcing_date[1]
end_date = forcing_date[length(forcing_date)-1]

ET_observed = get('data_ET')[get('data_ET')$date >= start_date & get('data_ET')$date <= 
                                 end_date, ]
ET_pot_observed = forc$ETpot

dischargeR = read.table("output/dischargeR.dat")
storageR = read.table("output/storageR.dat")
dateR = read.table("output/dateR.dat")

dv = read.table("output/dv.dat")
ETact = read.table("output/ETact.dat")
discharge = read.table("output/discharge.dat")

storageR = cbind(dateR, storageR)
dischargeR = cbind(dateR, dischargeR)

storageR$x = as.numeric(strptime(as.character(storageR$x), 
                                   format = "%Y%m%d%H", tz = "UTC"))

storageR = get('storageR')[get('storageR')$x >= start_date & get('storageR')$x <= 
                            end_date, ]


dischargeR$x = as.numeric(strptime(as.character(dischargeR$x), 
                                 format = "%Y%m%d%H", tz = "UTC"))

dischargeR = get('dischargeR')[get('dischargeR')$x >= start_date & get('dischargeR')$x <= 
                             end_date, ]

write.table(dischargeR,"output/dischargeRcut.dat")
write.table(storageR, "output/storageRcut.dat")

func_data_perc = function(data){
  d25 = apply(data, 1, quantile, 0.25)
  d50 = apply(data, 1, quantile, 0.5)
  d75 = apply(data, 1, quantile, 0.75)
  
  return(data.frame(d25, d50, d75))
}

###processing
QR = func_data_perc(dischargeR)
SR = func_data_perc(storageR)

QWQ = func_data_perc(discharge[,1:10])
QWE = func_data_perc(discharge[,10:20])
QWC = func_data_perc(discharge[,20:30])

library(hydroGOF)
print(KGE(sim=QWQ$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))
print(KGE(sim=QWE$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))
print(KGE(sim=QWC$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))

print(pbias(sim=QWQ$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))
print(pbias(sim=QWE$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))
print(pbias(sim=QWC$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))

print((rPearson(sim=QWQ$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))^2)
print((rPearson(sim=QWE$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))^2)
print((rPearson(sim=QWC$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))^2)

print(rmse(sim=QWQ$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))
print(rmse(sim=QWE$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))
print(rmse(sim=QWC$d50[8785:20424], obs=forc$Q[8785:20424], na.rm=T))

dvWQ = func_data_perc(dv[,1:10])
dvWE = func_data_perc(dv[,10:20])
dvWC = func_data_perc(dv[,20:30])


EQ = func_data_perc(ETact[2:nrow(ETact),1:10])
EE = func_data_perc(ETact[2:nrow(ETact),10:20])
EC = func_data_perc(ETact[2:nrow(ETact),20:30])

func_to_daily = function(data){
    group_indices <- rep(1:ceiling(length(data)/24), each = 24, length.out = length(data))

    return(as.numeric(tapply(data, group_indices, sum)))
}

EQd = func_to_daily(EQ$d50)
EEd = func_to_daily(EE$d50)
ECd = func_to_daily(EC$d50)


print(KGE(sim=EQd[367:851], obs=data_ET$ET[367:851], na.rm=T))
print(KGE(sim=EEd[367:851], obs=data_ET$ET[367:851], na.rm=T))
print(KGE(sim=ECd[367:851], obs=data_ET$ET[367:851], na.rm=T))

print(pbias(sim=EQd[367:851], obs=data_ET$ET[367:851], na.rm=T))
print(pbias(sim=EEd[367:851], obs=data_ET$ET[367:851], na.rm=T))
print(pbias(sim=ECd[367:851], obs=data_ET$ET[367:851], na.rm=T))

print((rPearson(sim=EQd[367:851], obs=data_ET$ET[367:851], na.rm=T))^2)
print((rPearson(sim=EEd[367:851], obs=data_ET$ET[367:851], na.rm=T))^2)
print((rPearson(sim=ECd[367:851], obs=data_ET$ET[367:851], na.rm=T))^2)

print(rmse(sim=EQd[367:851], obs=data_ET$ET[367:851], na.rm=T))
print(rmse(sim=EEd[367:851], obs=data_ET$ET[367:851], na.rm=T))
print(rmse(sim=ECd[367:851], obs=data_ET$ET[367:851], na.rm=T))

QWQ = func_data_perc(discharge[,1:10])
QWE = func_data_perc(discharge[,10:20])
QWC = func_data_perc(discharge[,20:30])
write.table(QWQ, "output/QWQ.dat")
write.table(QWE, "output/QWE.dat")
write.table(QWC, "output/QWC.dat")
write.table(forc, "output/forc.dat")
write.table(dvWQ, "output/dvWQ.dat")
write.table(dvWE, "output/dvWE.dat")
write.table(dvWC, "output/dvWC.dat")
write.table(SR, "output/SR.dat")
write.table(QR, "output/QR.dat")
#####################
# Change Q-h-relation
#####################

# Build a function in which the discharge is computed as a function of stage height.
# The arguments pars and hSmin are included even though they are not used in the function.
func_Q_hS_Hupsel = function(x, pars, hSmin)
{
  h = x/1000
  if(h <= 0)
  {
    0
  }else if(h < 0.2)
  {
    (10 ^(0.1645 + 2.0144 * log10(h) + 0.1342 * log10(h)^2)) *0.5538
  }else if(h <= 1.5)
  {
    (10 ^ (0.2741 + 2.3236 * log10(h) + 0.3511 * log10(h)^2)) *0.5538
  }else{
    2.769 + (10 ^(0.2741 + 2.3236 * log10(h-1.5) + 0.3511 * log10(h-1.5)^2)) *0.5538
  }
}

# Then set this function as the current stage-discharge relation.
set_func_Q_hS(func_Q_hS_Hupsel)

# And to check if you did it right:
show_func_Q_hS()

# To check the value for a specific stage height manually:
func_Q_hS_Hupsel(x=1000)


# ##########################
# # Monte Carlo calibration
# ##########################
# 
# # Load package to compute Kling-Gupta efficiency.
# library(hydroGOF)
# 
# # Make many (in this example 1000) random parameter sets within certain limits.
# cW  = runif(2000, min=1, max=2000)
# cG  = runif(2000, min=1e5, max=2e7)
# cQ  = runif(2000, min=1  , max=100)
# 
# # Make an empty vector to store mean sums of squares during for-loop.
# SS_P  = c()
# KGE_P = c()
# KGE_ET = c()
# KGE_COMP = c()
# 
# # Run a for-loop over all parameter sets and run WALRUS in every iteration.
# for(i in 1:2000)
# {
#   print(i)
#   # look up the parameter set for this iteration
#   parameters = data.frame(cW=cW[i], cG=cG[i], cQ=cQ[i], cV=4, cS=4,  
#                           dG0=1250, cD=1500, aS=0.01, st="loamy_sand")
# 
#   # run WALRUS
#   modeled    = WALRUS_loop(pars=parameters)
#   
#   ET_cal  = modeled$ETact
#   ET_cal  = ET_cal[2:length(ET_cal)]
#   group_indices <- rep(1:ceiling(length(ET_cal)/24), each = 24, length.out = length(ET_cal))
#   
#   ET_daily = as.numeric(tapply(ET_cal, group_indices, sum))
# 
#   
#   # Compute and store the mean sum of squares.
#   SS_P [i]     = mean((Qobs_forNS-modeled$Q)^2)
#   KGE_P [i]     = KGE(sim=modeled$Q, obs=Qobs_forNS, na.rm=T)
#   KGE_ET [i]     = KGE(sim=ET_daily, obs=ET_observed$ET, na.rm=T)
#   KGE_COMP [i]     = 0.5*KGE_P[i] + 0.5*KGE_ET[i]
# }
# 
# # Write the results to file: parameter values and belonging sum of squares.
# write.table(data.frame(cbind(cW, cG, cQ, SS_P, KGE_P, KGE_ET, KGE_COMP)), "output/pars_Hupsel_ETComp.dat")
# 


# ###############################
# # Parameters and initial values
# ###############################
# parameter = read.table("data/pars_Hupsel_ETComp.dat", header=TRUE)
# 
# #running for Q driven
# best_pars = parameter[order(-parameter$KGE_P)[1:10],]
# best_pars = rbind(best_pars,parameter[order(-parameter$KGE_ET)[1:10],])
# best_pars = rbind(best_pars,parameter[order(-parameter$KGE_COMP)[1:10],])
# 
# # make empty vector for modelled discharge
# Qmod = matrix(ncol=30, nrow=nrow(forc)+1)
# ETact = matrix(ncol=30, nrow=nrow(forc)+1)
# dV = matrix(ncol=30, nrow=nrow(forc)+1)
# 
# for(i in 1:30)
# {
#   print(i)
#   
#   # Define the parameters (cW, cV, cG, cQ, cS), initial conditions (dG0) and 
#   # catchment characteristics (cD, aS, soil type).
#   parameters = data.frame(cW=best_pars$cW[i], cG=best_pars$cG[i], cQ=best_pars$cQ[i], cV=4, cS=4,  
#                           dG0=1250, cD=2000, aS=0.1, st="loamy_sand")
#   
#   #####
#   # Run
#   #####
#   
#   # Run the model. 
#   mod = WALRUS_loop(pars=parameters)
#   
#   # save modelled Q in the matrix
#   Qmod[,i] = mod$Q
#   ETact[,i] = mod$ETact
#   dV[,i] = mod$dV
# }
# 
# write.table(Qmod, "output/discharge.dat")
# write.table(ETact, "output/ETact.dat")
# write.table(dV, "output/dv.dat")
# ######
# # Plot
# ######
# 
# # plot observed discharge
# plot(c(0,forc$Q), type="l")
# 
# # plot lines for eack parameter set
# for(i in 1:10)
# {
#   lines(Qmod[,i], col="dodgerblue")
# }
# 
# # draw observed discharge on top
# lines(c(0,forc$Q))
# 
# WALRUS_postprocessing(o=mod, pars=parameters, n=name)
# 
# # Define the parameters (cW, cV, cG, cQ, cS), initial conditions (dG0) and 
# # catchment characteristics (cD, aS, soil type).
# pars = data.frame(cW=200, cV=4, cG=5e6, cQ=10, cS=4,  
#                   dG0=1250, cD=1500, aS=0.01, st="loamy_sand")
# 
# #####
# # Run
# #####
# 
# # Run the model. 
# mod = WALRUS_loop(pars=pars)
# 
# 

# ##########################
# # Output files and figures
# ##########################
# 
# # Give the run a logical name. This will be used for the output data files and figures.
# name = "Qh_relation_Hupsel"
# 
# # Postprocessing: create datafiles and show figures.
# WALRUS_postprocessing(o=mod, pars=pars, n=name)
# 
# 
# #######################
# # Compare Q-h-relations
# #######################
# 
# # Set back to original function
# set_func_Q_hS(NULL)
# func_Q_hS_default = show_func_Q_hS()
# 
# # Rewrite functions for plotting.
# func_Q_hS_Hupsel_for_plot  = Vectorize(function(x){return(func_Q_hS_Hupsel(x))})
# func_Q_hS_default_for_plot = Vectorize(function(x){return(func_Q_hS_default(x,pars=pars,hSmin=0))})
# 
# # Plot two curves and add a legend.
# par(mfrow=c(1,1))
# curve(func_Q_hS_default_for_plot, col="dodgerblue", 0, 1600, ylim=c(0,5), xlab="h [mm]", ylab="Q [mm/h]")
# curve(func_Q_hS_Hupsel_for_plot, col="purple", add=TRUE)
# legend(c("Default","Hupsel"),col=c("dodgerblue","purple"), lty=1, bty="n", x="topleft")
# 

