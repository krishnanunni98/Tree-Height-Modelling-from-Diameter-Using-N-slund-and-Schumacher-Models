##########################################################################
#                                                                        #
                           #
#  Modeling of tree height                           
#  (UEF/For, No. 3513089, 8 ECTS)                                        #
#                                                                        #
                                  #
#                                                                        #
##########################################################################


rm(list=ls()) 			#tk: clean the workspace

setwd("c:/users/krish/downloads") 	#tk: set directory according to your folder location: same way as in the 3rd exercise

#=============================================#
#=# 1) Definition of user defined functions #=#
#=============================================#
#tk: Define a function to calculate Naeslund's height easily; see Lectures: Eq. 15 in section 4.3

Naeslund_h_f = function(parameters, d13.cm){
  value = -999
  value = 1.3 + (d13.cm^2)/((parameters[1] + parameters[2]*d13.cm)^2)
  value
}

#===========================================#
#tk: ...Naeslund with bias correction (Tailor's series)

Naeslund_h_f_bc = function(parameters, d13.cm, sigma){
  value = -999
  value = 1.3 + (d13.cm^2)/((parameters[1] + parameters[2]*d13.cm)^2) + (3*(d13.cm^2)*sigma^2)/((parameters[1] + parameters[2]*d13.cm)^4)
  value
}

#===========================================#
#tk: define a function for linear mixed effect model by Schumacher

Schumacher_lme = function(c0, c1, u0.j, bias.correction, d13.cm){
  value = -999  
  value = exp( c0 + c1*(d13.cm+3)^-1 + u0.j + bias.correction )
  value
}

#=================================#
#=# 2) Compilation of data etc. #=#
#=================================#
htdata <- read.table("tree_list.dat", header=TRUE, sep= "\t")  #tk: read the tree data
htdata <- htdata[order(htdata$plot_id, -htdata$d13.cm), ]      #tk: order by plot_id and then by diameter (largest to smallest)

# Create a new data.frame with sample trees only
moddat=subset(htdata, htdata$h.m>0)

# Calculate transformations
moddat$y = moddat$d13.cm/sqrt(moddat$h.m-1.3)

# Determine hierarchical id-variables as factors
moddat$plot_id = as.factor(moddat$plot_id)
moddat$tree_id = as.factor(moddat$tree_id)

names(moddat)
plot(moddat$d13.cm, moddat$h.m, xlab="d13.cm", ylab="h.m") #tk: plot diameter at breast height against height

#=======================================#
#=# 3) Estimation of model parameters #=#
#=======================================#

# Fitting linear regression models with the OLS using a forest standwise set of data
moddat1=subset(moddat, moddat$plot_id==10031051) #tk: select this plot (subjective decision)
moddat1 #tk: there are 19 trees on this plot

# Naeslund's (1937) stand height model

# Linear model form; model needs to be fit to estimate the parameters 
# that are used in the function and needed to draw curves
Naeslund.lm = lm(moddat1$y ~ moddat1$d13.cm, data=moddat1)     #tk: fit the model using the data from plot 10031051 only!
parameters = as.numeric(coefficients(Naeslund.lm))             #tk: extract the estimated regression coefficients
coef(summary(Naeslund.lm))[, "Std. Error"]                     #tk: inspect standard errors of coefficients 
sigma = summary(Naeslund.lm)$sigma                             #tk: extract residual standard error
obsdat1 = subset(htdata, htdata$plot_id == 10031051)           #tk: select all observations from one plot, i.e. all tallied trees         
dx.cm = seq(min(obsdat1$d13.cm), max(obsdat1$d13.cm), by=.1)   #tk: create a sequence from smallest to largest d13 with 0.1 cm intervals; this is needed to produce graphs

plot(moddat1$d13.cm, moddat1$h.m,                              #tk: plot d13 against height in modelling data
     xlim=c(0, max(moddat1$d13.cm+5)), 
     ylim=c(0, max(moddat1$h.m+5)),
     xlab="d13.cm",
     ylab="h.m")

lines(dx.cm, Naeslund_h_f(parameters, dx.cm),                  #tk: add blue curve for Naeslund_h_f model by calling it within the lines-function
      col="blue", lty="dashed", lwd=2)

lines(dx.cm, Naeslund_h_f_bc(parameters, dx.cm, sigma),        #tk: add green curve for bias corrected Naeslund's model
      col="green", lty="dashed", lwd=2)

# Nonlinear model form
Naeslund.nls = nls(h.m ~ 1.3 + d13.cm^2/(c0 + c1*d13.cm)^2,    #tk: Fit nonlinear model; Eq. 15 in section 4.3
                   start=c(c0=1, c1=.1), data= moddat1)

parameters = as.numeric(coefficients(Naeslund.nls))            #tk: extract the parameters

lines(dx.cm, Naeslund_h_f(parameters, dx.cm),                  #tk: add red curve to the previous plot
      col="red", lty="dashed", lwd=2)

# General linear mixed-effects models
# A variance component model: ln(h.ji) = c0 + c1/(d.ji+3) + u0.j + e0.ji

library(nlme) 						#tk: this package is needed to fit LME models

moddat2 = subset(htdata, htdata$h.m>0) 			#tk: create a subset with trees for which the height has been measured (1749 trees)
moddat2$lnh = log(moddat2$h.m)         			#tk: create natural logarithmic transformation for height
moddat2$invd3 = 1/(moddat2$d13.cm+3)   			#tk: create inverse transformation. 
#tk: When d13 -> 0, 1/d -> infinity. 
#tk: Therefore, small diameters have large influence on the estimated parameters of the model.
#tk: To decrease the effects of small diameters, intercept +3 was applied. 

#tk: Fit the lme model
Schum.lme0 = lme( lnh ~ invd3, data= moddat2, random=~1 | plot_id )               #tk: technically these two are the same
Schum.lme = lme( log(h.m) ~ I(1/(d13.cm+3)), data= moddat2, random=~1 | plot_id ) #tk: ... i.e., it is possible to create the transformations within the model fitting function

#-# Completed by 20/03/2013

obsdat1 = subset(htdata, htdata$plot_id==10031051) 	#tk: same plot as before as an example (all 27 trees)
obsdat1c = subset(obsdat1, obsdat1$h.m>0)          	#tk: a subset with only trees with measured heights (19 trees)

summary(Schum.lme)                               	#tk: model summary
VarCorr(Schum.lme)                               	#tk: Inspect variance and correlation components
var.u0=as.numeric(VarCorr(Schum.lme)[1] )        	#tk: within plot variance
var.e0=as.numeric(VarCorr(Schum.lme)[2] )        	#tk: between plot variance
c0=as.numeric(fixef(Schum.lme))[1]               	#tk: intercept of model
c1=as.numeric(fixef(Schum.lme))[2]               	#tk: coefficient for predictor
n.c=nrow(obsdat1c)                               	#tk: number of observations (=19)

#tk: calculate random effect for plot 10031051 and bias correction
u0.j = var.u0/(var.u0+var.e0/n.c) * (mean(log(obsdat1c$h.m))-mean(c0+c1*(obsdat1c$d13.cm+3)^-1)) #tk: predictor for random effect; see equation 27.1 in section 5.2.2
var.pe.u0.j = var.u0/(var.u0+var.e0/n.c) * var.e0/n.c
bias.correction = 1/2*(var.pe.u0.j+var.e0)

ranef(Schum.lme)                                             	#tk: Check that u0.j was correct; 0.204357170 for plot 10031051 
dx.cm = seq(min(obsdat1$d13.cm), max(obsdat1$d13.cm), by=.1) 	#tk: create a sequence from smallest to largest d13 with 0.1 cm intervals
hx.m = Schumacher_lme(c0, c1, u0.j, bias.correction, dx.cm)  	#tk: predict the localised heights using the function that was defined in the beginning of this R-code
plot(obsdat1c$d13.cm, obsdat1c$h.m, xlim=c(0,max(obsdat1c$d13.cm+5)), ylim=c(0,max(obsdat1c$h.m+5)), xlab="d13.cm", ylab="h.m")	#tk: plot observed d13 against height
lines(dx.cm, hx.m, col="blue", lty=1, lwd=2)                 	#tk: add line for the localised Schumacher's model

parameters = as.numeric(coefficients(Naeslund.nls))          	#tk: extract coefficients of the Naeslunds nonlinear model; these lines 144-146 are technically the same as lines 99-102
hx.m = Naeslund_h_f(parameters, dx.cm)                       	#tk: predict with the Naeslund's function
lines(dx.cm, hx.m, col="red", lty=1, lwd=2)                  	#tk: ...and add a red curve to the plot

#tk: End of the R script