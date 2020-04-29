########## Start of header ################
# Title:  Univariate Linear Mixed Model - Muttama 2019 0-5cm - trial 30m grid
#  
# Description: 
# Testing Linear Mixed Model for single response variable - SOC
#
#
# Project: Muttama_data
# Author: Kate Coelli  kate.coelli@sydney.edu.au
# Date: 7/04/2020
# Last Updated: 29/04/2020
#
########## End of header ################

#################################################
################# packages ######################
#################################################

library(geoR)
library(moments)
library(dplyr)
library(epiR)
library(ModelMetrics)
library(lattice)
library(colorspace)
library(raster)

# library(gstat)
# library(sp)
# library(graphics)


# library(rgdal)
# library(raster)

# library(ggplot2)
# library(ggspatial)
# library(rasterVis)

# library(graphics)

#################################################
########### data import and selection ###########
#################################################


#import dataframe of covariates (also includes variable to be predicted)
#this dataframe was created using Muttama_data_cube_splined script
setwd("C:/Users/kcoe7598/Dropbox (Sydney Uni)/R/soil_carbon_modelling_PhD/Muttama_data")
Muttama_point_df<- read.csv("data/Muttama_SOC_and_covariates.csv") #if not already in workspace

#filter data for Univariate mixed model = 1 layer, 1 time point
#select relevant columns for LMM

pred.soil<- Muttama_point_df%>%
  filter(Year == "2013" & Month == 6 & Upper == "0")%>%
  dplyr::select(Easting, Northing,  BD, radiometrics, DEM, rain_annual_avge, EVI_max, EVI_min, NDVI_max, NDVI_min, clay., SOC) # - replace clay and BD??? currently SLGA


#################################################
############### explore the data ################
#################################################

#view "skewness" of data
hist(pred.soil$SOC) #distribution is normal - transform if data not normally distributed
boxplot(pred.soil$SOC) #further demonstrates normality

#identify which columns are the covariates
str(pred.soil)

#if skewness is between -1.5-1.5 then does not need transforming
skewness(pred.soil[,1:11])#check skewness of covariates

# All covariates are between -1.5 and 1.5 so don't need transforming

#look at correlation between covariates and the response variable(SOC)
#cor(pred.soil[predicted], pred.soil[predictors])

cor(pred.soil[,12], pred.soil[,1:11]) #pearson's correlation 

cor(pred.soil[,12], pred.soil[,1:11], method = "spearman") #spearman correlation

str(pred.soil)

#visualise relationships between covariates and variable to be predicted
pairs(cbind(pred.soil[,12], pred.soil[,1:2]))# SOC and Eastings and Northings

plot(pred.soil$rain_annual_avge, pred.soil$SOC)#annual average for the 5 years prior to the sample year

pairs(pred.soil)#all covariates and SOC


#create new dataset to alter
pred0<- pred.soil


#rescaled for variogram- now between 0 and 100 -> scale is in kilometers from lowest easting or northing respectively
pred0$x<-(pred0$Easting-min(pred0$Easting))/1000
pred0$y<-(pred0$Northing-min(pred0$Northing))/1000


#################################################
############### create variogram ################
#################################################

#look at structure of pred0 to find column numbers
str(pred0)

#leave out eastings and northings as x and y are built from these 
OC.geodata <- as.geodata(pred0,coords.col=13:14,data.col=12, covar.col=c(3:11,13:14), covar.names = "obj.names")#coordinates are also covariates

summary(OC.geodata)
 
plot(OC.geodata)

#Experimental semivariogram
#create variogram with all covariates
OC.var <- variog(OC.geodata, trend = ~ BD+radiometrics+DEM+rain_annual_avge+EVI_max+EVI_min+NDVI_max+NDVI_min+clay.+x+y, option = "bin")

#view variogram details
str(OC.var)

#format gridplot
par(mfrow=c(1,1))

#plot variogram
plot(OC.var,pch = 19, col= "blue", main="",xlim=c(0,50),ylim=c(0,0.5),xlab="distance(km)")
#need to check - semivariance should not decreases as you get further away from the measurement


#################################################
################# fit LMM model #################
#################################################

#use variogram plot and values to estimate range and sill values for maximum likelihood test
#best guess
range <-  10 #eyeball distance where ~logarithmic curve would flatten out
nugget <- 0.2 #use ~semivariance of point closest to y-axis
sill <- 0.4 #eyeball semivariance where ~logarithmic curve starts to flatten our
partial_sill<- sill-nugget

#sensitivity test
range_50<- 0.5*range
range_150<- 1.5*range
partial_sill_50<- 0.5*partial_sill
partial_sill_150<- 1.5*partial_sill

#fitting likfit function for combination of different range and partial sill to test sensitivity
ini.<- expand.grid(c(partial_sill, partial_sill_50, partial_sill_150), c(range, range_50, range_150))
par(mfrow=c(3,3))#plot parameters - alter depending on how many iterations in sensitvity test (rows in ini.)

#fit using model = exponential
for(i in 1:nrow(ini.)){
  lmm.exp <- likfit(OC.geodata, trend= ~BD+radiometrics+DEM+rain_annual_avge+EVI_max+EVI_min+NDVI_max+NDVI_min+clay.+x+y, #insert all covariates in trend
                  lambda=1,  ini.cov.pars= ini.[i,], fix.nugget = FALSE, nugget = nugget, 
                  lik.method = "REML", cov.model = "exp")
  print(summary(lmm.exp))
  plot(OC.var,pch = 19, col= "blue", main=paste0("set", i), xlab="distance(m)", xlim=c(0,50),ylim=c(0,0.8))#
  lines(lmm.exp, col="red")
 }


#fit using model = spherical
for(i in 1:nrow(ini.)){
  lmm.sph <- likfit(OC.geodata, trend= ~BD+radiometrics+DEM+rain_annual_avge+EVI_max+EVI_min+NDVI_max+NDVI_min+clay.+x+y, #insert all covariates in trend
                    lambda=1,  ini.cov.pars= ini.[i,], fix.nugget = FALSE, nugget = 0.18, 
                    lik.method = "REML", cov.model = "sph")
  print(summary(lmm.sph))
  plot(OC.var,pch = 19, col= "blue", main=paste0("set", i),xlab="distance(m)",xlim=c(0,50),ylim=c(0,0.8))
  lines(lmm.sph, col="green")
}



##### sensitivity test made very little difference so just use initial estimates when refitting the model
# I just chose exponential since both similar

#refit model

lmm.exp<- likfit(OC.geodata, trend= ~BD+radiometrics+DEM+rain_annual_avge+EVI_max+EVI_min+NDVI_max+NDVI_min+clay.+x+y, #insert all covariates in trend
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")

summary(lmm.exp)
# want model with smaller AIC
# AIC = variation not explained + number of parameters
# penalty for more parameters -> simpler better


#fitted variogram

par(mfrow = c(1,1))
plot(OC.var,pch = 19, col= "blue", main="",xlab="distance(m)", xlim=c(0,50),ylim=c(0,0.6))#make plot limits sensible
lines(lmm.exp, col="red")


###Variable selection
#variables for p-value function
lmm.exp
pred0
OC.geodata

### function to extract regression coefficients from likfit model - Thanks Liana :) ###

pvalsummary <- function(likfitmodel, dataframe, geodata){
  
  ##Code for getting P-values from geoR object
  
  fitmodel <- likfitmodel
  datas <- dataframe
  
  ## get coefficients
  coefficients <- fitmodel$beta  #geoR- dont get P-val for fixed effects
  ## get standard errors
  se_error <- sqrt(diag(fitmodel$beta.var))
  ## get t values
  t_value <- coefficients/se_error
  ## and probabilities
  t_prob <- 2 * pt(-abs(t_value), df = (nrow(datas) -length (coefficients)))
  #Confidence intervals
  tcrit<-abs(qt(0.025,lower.tail=F,df=(nrow(datas) -length (coefficients))))
  U95 <- coefficients + (tcrit*se_error)
  L95 <- coefficients - (tcrit*se_error)
  
  ## make pretty
  coef_mat <- cbind(coefficients, L95, U95,  se_error, t_value, t_prob)
  colnames(coef_mat) <- c("Estimate","Lower 95CI", "Upper 95CI",
                          "Std.Err","t value", "Pr(>|t|)" )
  rownames(coef_mat) <- colnames(model.matrix(fitmodel$trend, geodata))
  printCoefmat(coef_mat)
  
}


#use function

pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)


#drop covariate with highest p-value = y
#then rerun lmm.exp
#stop when all pvalues are below 0.1

lmm.exp<- likfit(OC.geodata, trend= ~BD+radiometrics+DEM+rain_annual_avge+EVI_max+EVI_min+NDVI_max+NDVI_min+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")


summary(lmm.exp)


###Variable selection

pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)

#drop covariate with highest p-value = radiometrics

lmm.exp<- likfit(OC.geodata, trend= ~BD+DEM+rain_annual_avge+EVI_max+EVI_min+NDVI_max+NDVI_min+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")


summary(lmm.sph)


###Variable selection

pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)

#drop covariate with highest p-value = EVI_min

lmm.exp<- likfit(OC.geodata, trend= ~BD+DEM+rain_annual_avge+EVI_max+NDVI_max+NDVI_min+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")

summary(lmm.exp)


###Variable selection
pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)


#drop covariate with highest p-value =EVI_max.

lmm.exp<- likfit(OC.geodata, trend= ~BD+DEM+rain_annual_avge+NDVI_max+NDVI_min+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")
summary(lmm.exp)


###Variable selection
pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)

#drop covariate with highest p-value =rain_anual_avge.

lmm.exp<- likfit(OC.geodata, trend= ~BD+DEM+NDVI_max+NDVI_min+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")
summary(lmm.exp)


###Variable selection

pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)


#drop covariate with highest p-value =NDVI_max

lmm.exp<- likfit(OC.geodata, trend= ~BD+DEM+NDVI_min+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")
summary(lmm.exp)


###Variable selection
pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)




#drop covariate with highest p-value = Bulk Density.

lmm.exp<- likfit(OC.geodata, trend= ~DEM+NDVI_min+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")
summary(lmm.exp)


###Variable selection

pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)

#drop covariate with highest p-value = NDVI_min.

lmm.exp<- likfit(OC.geodata, trend= ~DEM+clay.+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")
summary(lmm.exp)


###Variable selection

pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)




#drop covariate with highest p-value = clay.

lmm.exp<- likfit(OC.geodata, trend= ~DEM+x, 
                 lambda=1,  ini.cov.pars= c(partial_sill, range), fix.nugget = FALSE, nugget = nugget, 
                 lik.method = "REML", cov.model = "exp")
summary(lmm.exp)


###Variable selection

pvalsummary(likfitmodel = lmm.exp, dataframe = pred0, geodata = OC.geodata)



#final coefficient matrix


#               Estimate    Std.Err t value Pr(>|t|)   
# (Intercept) -0.7746122  0.6759904 -1.1459 0.257086   
# DEM          0.0042002  0.0014601  2.8766 0.005817 **
#   x            0.0448558  0.0143289  3.1305 0.002862 **






################################################
############### Cross-Validation ###############
################################################


#######Cross-validation  by kriging in geoR#######

xv.REML <- xvalid(OC.geodata, model=lmm.exp)

#Plot observed vs predicted
plot(xv.REML$data,xv.REML$predicted, main = "Observed vs Predicted SOC - GeoR", xlab = "Observed SOC (%)", ylab = "Predicted SOC (%)")
abline(0,1)

#r - correlation
cor(xv.REML$data,xv.REML$predicted)

#Lins concorodance correlation coefficient
concord<-epi.ccc(xv.REML$data,xv.REML$predicted)
concord$rho.c

#SSPE - want mean to be 1 and median to be 0.445
sspe<-(xv.REML$error)^2/xv.REML$krige.var
summary(sspe)

#Bias
mean(xv.REML$error)


#RMSE
sqrt(mean(xv.REML$error^2))



########## Cross validation by kriging in geostat ########
#adapted from liana's code

#convert data to gstat dataset
# dataset
head(pred0)#identify coordinate columns
pred0.gstat<-SpatialPointsDataFrame(pred0[1:2],pred0)

#Convert the geoR model to gstat
exp.gstat <- as.vgm.variomodel(lmm.exp)

#THINK JUST SIGNIFICANT COVARS BUT NEED TO CHECK
out <- krige.cv(SOC ~ DEM+x, pred0.gstat, exp.gstat, nfold=nrow(pred0)) #leave one out CV - can also set folds

#Plot observed vs predicted
plot(out$observed, out$var1.pred, main = "Observed vs Predicted SOC - Geostat", xlab = "Observed SOC (%)", ylab = "Predicted SOC (%)" )
abline(0,1,col = "red")


#r - correlation
cor(out$observed, out$var1.pred)

#Lins concorodnace correlation coefficient
concord <- epi.ccc(out$observed, out$var1.pred)
concord$rho.c      #concordance correlation coefficient


#SSPE (theta)
sspe <- (out$residual)^2/out$var1.var
summary(sspe)

#RMSE  
dfout <- as.data.frame(cbind(out@coords, out@data))   #remove NAs
dfout <- (dfout[complete.cases(dfout),])
summary(is.na(dfout))

rmse(dfout$observed, dfout$var1.pred)


################################################
################## kriging ####################
################################################

#File for prediction
head(Muttama_grid)
Muttama_grid <- read.table("data/Muttama_covariates_grid_trial.csv", header=T, sep=",") #incase not already in workspace


#filter data for 1 layer, 1 time point
#select relevant columns

Muttama_grid<- Muttama_grid%>%
  #filter(Year == "2013" & Month == 6 & Upper == "0")%>%
  select(Easting, Northing,  BD, radiometrics, DEM, rain_annual_avge, EVI_max, EVI_min, NDVI_max, NDVI_min, clay.) # - replace clay and BD??? currently SLGA


#rescaled for variogram- now between 0 and 100 -> scale is in kilometers from lowest easting or northing respectively
Muttama_grid$x <- (Muttama_grid$Easting-min(Muttama_grid$Easting, na.rm=T))/1000
Muttama_grid$y <- (Muttama_grid$Northing-min(Muttama_grid$Northing, na.rm=T))/1000

# dataset
head(pred0) 

#see colour palettes available 
hcl_palettes(plot=T)

xyplot(pred0$Northing~pred0$Easting, groups=c(pred0$SOC.<1, between(pred0$SOC, 1, 2), pred0$SOC>2),
       pch=16, lwd = 3, col = diverging_hcl(3, palette = "Blue-Red"), xlab="Easting (m)", ylab="Northing (m)", 
       key = list(space="right", points=list(col=diverging_hcl(3, palette = "Blue-Red"), pch=16, lty=c(3,2), lwd=0),
                                text=list(c("<1"," 1<2", "2<"), border=TRUE)), scales=list(tck=-0.5)
       )


# set as gstat objects for kriging
#data
head(pred0.gstat)

#grid for predictions
head(Muttama_grid)#identify coordinate columns

Muttama_grid.gstat<-SpatialPointsDataFrame(Muttama_grid[1:2],Muttama_grid)#gstat object


######### predict onto grid points ##########
eblup<-krige(SOC~DEM+x, pred0.gstat,
             newdata=Muttama_grid.gstat,model = exp.gstat)

eblup

summary(eblup$var1.pred)
spplot(eblup["var1.pred"])

################################################
############## convert to raster ###############
################################################


raster_Muttama_trial_grid<-eblup
gridded(raster_Muttama_trial_grid) <- TRUE
raster_Muttama_trial_grid<-raster(raster_Muttama_trial_grid)
raster_Muttama_trial_grid

crs(raster_Muttama_trial_grid)<- UTM
raster_Muttama_trial_grid<- projectRaster(from = raster_Muttama_trial_grid, crs = UTM) #choose a raster for "to" in the desired projection - NDVI in this case


writeRaster(raster_Muttama_trial_grid,
            "data/raster_Muttama_trial_grid.tif", 
            overwrite=T)

