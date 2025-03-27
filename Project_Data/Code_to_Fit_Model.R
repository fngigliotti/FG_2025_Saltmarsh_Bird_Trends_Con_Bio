############################################################
###R Code to Fit JAGS Model for Population Trend Analysis###
############################################################

##Gigliotti et al. 2025: Short-term stabilization of population 
##trends raises questions about the future of saltmarsh birds 
##in eastern North America

############################
###Load Required Packages###
############################

library(jagsUI)
#Could fit with other MCMC samplers as well if desired

###############
###Load Data###
###############

#R object containing all data needed for fitting models
load(file = "Unscaled_Saltmarsh_Bird_Data.RDA")

############################
###Process Data for Model###
############################

#Thin to focal species (e.g., saltmarsh sparrow)
saltmarsh.bird.data$nspecies<-NULL
saltmarsh.bird.data$y<-saltmarsh.bird.data$y[,"SALS",,]
saltmarsh.bird.data$y.time<-saltmarsh.bird.data$y.time[,"SALS",,,]

##No sites are dropped when thinning to saltmarsh sparrow, but thinning to other
##species will require thinning of model covariates as well (to match number of sites
##in survey data)

#Collapse y.time down to just first detection of focal species at any point

saltmarsh.bird.data$y.time<-apply(saltmarsh.bird.data$y.time,c(1:3),function(x) min(x, na.rm = T))
#Fill Inf values with NA (no individuals detected on survey)
saltmarsh.bird.data$y.time[saltmarsh.bird.data$y.time=="Inf"]<-NA
saltmarsh.bird.data$ttd<-saltmarsh.bird.data$y.time #rename for model

#Create detection array
d.array<-1*!is.na(saltmarsh.bird.data$ttd)
saltmarsh.bird.data$d<-d.array

#Initialize sites as occupied by species if species ever observed there
zi.st<-apply(d.array, c(1,3), max, na.rm = TRUE)

#Initialize abundance data
Nst<-apply(saltmarsh.bird.data$y, c(1,3), max, na.rm = TRUE) + 2
Nst[Nst=="-Inf"]<-2 #abundance initialized at 2 at unobserved points 
#                   (to prevent issue with parent values in JAGS processor)

#Fill in 0s at 'unoccupied' sites
Nst[zi.st==0]<-0

#Initialize ttd data
ttdst <-d.array
ttdst[] <- NA
ttdst[d.array == 0] <- saltmarsh.bird.data$nwindows + 1

#Initals (can initialize other parameters if desired to expedite convergence time)
inits <- function(){list(N=Nst, zi = zi.st, ttd=ttdst)}

###############################################
###Center and scale all covariates for model###
###############################################

#UVVR
saltmarsh.bird.data$UVVR$Mean<-as.numeric(scale(saltmarsh.bird.data$UVVR$Mean))
saltmarsh.bird.data$UVVR$Trend<-as.numeric(scale(saltmarsh.bird.data$UVVR$Trend))

##latitudes
saltmarsh.bird.data$latitudes<-as.numeric(scale(saltmarsh.bird.data$latitudes))

##SLR
saltmarsh.bird.data$tides<-as.numeric(scale(saltmarsh.bird.data$tides))

##max yearly tide
saltmarsh.bird.data$max.tide<-scale(saltmarsh.bird.data$max.tide)

##patch area
saltmarsh.bird.data$patch.area<-as.numeric(scale(saltmarsh.bird.data$patch.area))

##road density
saltmarsh.bird.data$patch.road.dens<-as.numeric(scale(saltmarsh.bird.data$patch.road.dens))

##vegetation
saltmarsh.bird.data$veg.props<-scale(saltmarsh.bird.data$veg.props)

#noise
noise<-apply(saltmarsh.bird.data$noise,c(2,3),as.numeric)
noise<-array(scale(noise),dim = c(saltmarsh.bird.data$nsites,2,9))
dimnames(noise)<-dimnames(saltmarsh.bird.data$noise)
saltmarsh.bird.data$noise<-noise

#time since sunrise
sunrise<-array(scale(saltmarsh.bird.data$sunrise),dim = c(saltmarsh.bird.data$nsites,2,9))
dimnames(sunrise)<-dimnames(saltmarsh.bird.data$sunrise)
saltmarsh.bird.data$sunrise<-sunrise

#sky
sky<-array(scale(saltmarsh.bird.data$sky),dim = c(saltmarsh.bird.data$nsites,2,9))
dimnames(sky)<-dimnames(saltmarsh.bird.data$sky)
saltmarsh.bird.data$sky<-sky

#temperature 
temp<-array(scale(saltmarsh.bird.data$temp),dim = c(saltmarsh.bird.data$nsites,2,9))
dimnames(temp)<-dimnames(saltmarsh.bird.data$temp)
saltmarsh.bird.data$temp<-temp

#Julian day
j.day<-array(scale(saltmarsh.bird.data$j.day),dim = c(saltmarsh.bird.data$nsites,2,9))
dimnames(j.day)<-dimnames(saltmarsh.bird.data$j.day)
saltmarsh.bird.data$j.day<-j.day

##################################
###Specify Parameters Monitored###
##################################

params<-c("alpha.lambda", "beta.lambda", "sd.patch",
          "sd.observer", "alpha.p", "beta.p",
          "mean.det", "omega", "N.tot", "N.trend.global", "N.trend.state",
          "state.trend", "mu.state", "sd.state", "log.lik")

###############
###Fit Model###
###############
model.fit<-jagsUI::autojags(data = saltmarsh.bird.data, inits = inits,
                 parameters.to.save = params, model.file = 'Saltmarsh_Bird_Population_Change_Model.txt',
                 n.chains = 4, n.adapt=1000, iter.increment=5000,
                 n.burnin=10000, n.thin=5, Rhat.limit = 1.1,
                 max.iter=100000, parallel = TRUE)

#######################
###End of Model Code###
#######################