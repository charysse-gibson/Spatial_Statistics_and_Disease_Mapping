######################################################
## Behavioral Risk Factors Explaining Excess Stroke Mortality Rates in Texas Counties: A Stroke Belt Study
## BST 6200: Spatial Statistics and Disease Mapping
## Charysse Gibson
## 4/29/2020
######################################################

setwd('C:/Users/chary/OneDrive/Documents/SLU/Sping 2020/BST 6200 spatial Statistics and Disease Mapping/Final project/')

library( tigris )
library( sp )
library( ggplot2 )
library( ggthemes )
library( tmap )
library(tmaptools)
library( dplyr )
library( tidyverse )
library( sf )
library( GISTools )
library( spdep )
library( SpatialEpi )
library( rgeos )

Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

##----import data----

TX <- counties("Texas" , cb=TRUE )

TXPopulation <- read.csv("TXPop.csv", stringsAsFactors = FALSE)
TXStroke <- read.csv("TXStroke.csv", stringsAsFactors = FALSE)
TXSmoking <- read.csv("TXSmoking.csv", stringsAsFactors = FALSE)
TXPhysical <- read.csv("TXPhysical.csv", stringsAsFactors = FALSE)
TXDrinking <- read.csv("TXDrinking.csv", stringsAsFactors = FALSE)

TXPopulation$Value <- as.integer(TXPopulation$Value)

# check number of counties
# nrow(TX@data)
# nrow(TXStroke)
# nrow(TXSmoking)
# nrow(TXPhysical)
# nrow(TXDrinking)

TX@data <-
  left_join(x=TX@data,y=TXPopulation,by=c('NAME'='County'))
TX@data <-
  rename(TX@data, population=Value)

TX@data <-
  left_join(x=TX@data,y=TXStroke,by=c('NAME'='County'))
TX@data <-
  rename(TX@data, stroke=Value)
TX@data$stroke_count <- as.integer(ceiling((TX@data$stroke * TX@data$population)/100000))

TX@data <-
  left_join(x=TX@data,y=TXSmoking,by=c('NAME'='County'))
TX@data <-
  rename(TX@data, smoking=Value)
# TX@data$smoking <- TX$smoking*100
TX@data$smoking <- TX$smoking

TX@data <-
  left_join(x=TX@data,y=TXPhysical,by=c('NAME'='County'))
TX@data <-
  rename(TX@data, physical=Value)
# TX@data$physical <- TX$physical*100
TX@data$physical <- TX$physical

TX@data <-
  left_join(x=TX@data,y=TXDrinking,by=c('NAME'='County'))
TX@data <-
  rename(TX@data, drinking=Value)
# TX@data$drinking <- TX$drinking*100
TX@data$drinking <- TX$drinking

# change from st to sf 
TX.sf = st_as_sf( TX )

##----descriptive stats----

library(tableone)
myvars <- c("stroke","smoking","physical","drinking")
tab1 <- CreateTableOne(vars = myvars, data = TX@data)
tab1
summary(tab1) #includes details and standardized mean differences


##----chloropleth maps----

# windows(30,30)
qtm(TX.sf, fill="stroke", text="NAME", text.size=.7,
    format="World_wide", style="classic", 
    main.title='Figure 1. Chloropleth Map of Stroke Mortality Rate (per 100,000)',
    main.title.size=0.9, text.root=5, fill.title="Stroke Mortality Rate (per 100,000)",
    fill.style="fixed",fill.breaks=seq(0,200,40))

# windows(6.5,6.5)
tm_shape( TX.sf ) +
  tm_fill(col="stroke", title = "Stroke Mortality Rate (per 100,000)", style="fixed", breaks=seq(0,200,25)) +
  tm_borders( col="black" )
# windows(6.5,6.5)
tm_shape( TX.sf ) +
  tm_fill(col="stroke", title = "Stroke Mortality Rate (per 100,000)", style="fixed", breaks=seq(0,200,25)) +
  tm_borders( col="black" )
# windows(6.5,6.5)
tm_shape( TX.sf ) +
  tm_fill(col="stroke", title = "Stroke Mortality Rate (per 100,000)", style="fixed", breaks=seq(50,200,25)) +
  tm_borders( col="black" )

# windows(4.5,4.5)
tm_shape( TX.sf ) +
  tm_fill(col="smoking", title = "% Smokers", style="fixed", breaks=seq(10,40,5) ) +
  tm_borders( col="black" )

# windows(4.5,4.5)
tm_shape( TX.sf ) +
  tm_fill( col="physical", title = "% Physically Inactive", style="fixed" , breaks=seq(10,40,5) ) +
  tm_borders( col="black" )

# windows(4.5,4.5)
tm_shape( TX.sf ) +
  tm_fill( col="drinking", title = "% Excessive Drinking", style="fixed" , breaks=seq(10,40,5) ) +
  tm_borders( col="black" )

##----Moran's I----

TX.nb = poly2nb( TX.sf , queen=FALSE )
TX.lw = nb2listw( TX.nb )

#----stroke----

TX.stroke.moran.test = moran.test( TX.sf$stroke, TX.lw )
print( TX.stroke.moran.test )

I = TX.stroke.moran.test$estimate[[1]]
expI = TX.stroke.moran.test$estimate[[2]]
varI = TX.stroke.moran.test$estimate[[3]]
z.stat = ( I - expI ) / sqrt(varI)
z.stat
p.value = pnorm( z.stat , lower.tail=FALSE )

# Monte Carlo Simulation of Moran's I
TX.stroke.moran.sim = moran.mc( TX.sf$stroke, TX.lw , 99999 )
str( TX.stroke.moran.sim )
simulated.I = TX.stroke.moran.sim$res
# histogram of simulation
windows(9,5)
hist( simulated.I , xlim=c(-1,1) , breaks=seq(-1,1,0.01) , col="blue", xlab = 'Simulated Moran\'s I', main='Monte Carlo Simulation for Moran\'s I')
abline( v=TX.stroke.moran.sim$statistic , col='red')

#----smoking----

TX.smoking.moran.test = moran.test( TX.sf$smoking, TX.lw )
print( TX.smoking.moran.test )

I = TX.smoking.moran.test$estimate[[1]]
expI = TX.smoking.moran.test$estimate[[2]]
varI = TX.smoking.moran.test$estimate[[3]]
z.stat = ( I - expI ) / sqrt(varI)
z.stat
p.value = pnorm( z.stat , lower.tail=FALSE )

# Monte Carlo Simulation of Moran's I
TX.smoking.moran.sim = moran.mc( TX.sf$smoking, TX.lw , 99999 )
str( TX.smoking.moran.sim )
simulated.I = TX.smoking.moran.sim$res
# histogram of simulation
windows(9,5)
hist( simulated.I , xlim=c(-0.5,0.5) , breaks=seq(-0.5,0.5,0.01) , col="blue", xlab = 'Simulated Moran\'s I', main='Monte Carlo Simulation for Moran\'s I')
abline( v=TX.smoking.moran.sim$statistic , col='red')


#----physcial inactivity----

TX.physical.moran.test = moran.test( TX.sf$physical, TX.lw )
print( TX.physical.moran.test )

I = TX.physical.moran.test$estimate[[1]]
expI = TX.physical.moran.test$estimate[[2]]
varI = TX.physical.moran.test$estimate[[3]]
z.stat = ( I - expI ) / sqrt(varI)
z.stat
p.value = pnorm( z.stat , lower.tail=FALSE )

# Monte Carlo Simulation of Moran's I
TX.physical.moran.sim = moran.mc( TX.sf$physical, TX.lw , 99999 )
str( TX.physical.moran.sim )
simulated.I = TX.physical.moran.sim$res
# histogram of simulation
windows(9,5)
hist( simulated.I , xlim=c(-0.5,0.5) , breaks=seq(-0.5,0.5,0.01) , col="blue", xlab = 'Simulated Moran\'s I', main='Monte Carlo Simulation for Moran\'s I')
abline( v=TX.physical.moran.sim$statistic , col='red')

#----excessive drinking----

TX.drinking.moran.test = moran.test( TX.sf$drinking, TX.lw )
print( TX.drinking.moran.test )

I = TX.drinking.moran.test$estimate[[1]]
expI = TX.drinking.moran.test$estimate[[2]]
varI = TX.drinking.moran.test$estimate[[3]]
z.stat = ( I - expI ) / sqrt(varI)
z.stat
p.value = pnorm( z.stat , lower.tail=FALSE )

# Monte Carlo Simulation of Moran's I
TX.drinking.moran.sim = moran.mc( TX.sf$drinking, TX.lw , 99999 )
str( TX.drinking.moran.sim )
simulated.I = TX.drinking.moran.sim$res
# histogram of simulation
windows(9,5)
hist( simulated.I , xlim=c(-0.5,0.5) , breaks=seq(-0.5,0.5,0.01) , col="blue", xlab = 'Simulated Moran\'s I', main='Monte Carlo Simulation for Moran\'s I')
abline( v=TX.drinking.moran.sim$statistic , col='red')

##----conditional autoregressive model----

library( nimble )
library( methods )
str(TX@data)

TX.nb = poly2nb( TX.sf )
MO.net = nb2lines( TX.nb , coords=coordinates(TX) )
TX.lw = nb2listw( TX.nb )

k = length(TX@data$NAME)      ## There are 254 counties

num = rep(0,k) 
for (i in 1:k) num[i] = length( TX.lw$neighbours[[i]] )
adj = c() 
for (i in 1:k) adj = c(adj,TX.lw$neighbours[[i]] )

L = length(adj)

TX.Code = nimbleCode({
  alpha ~ dflat()
  beta1 ~ dnorm( 0 , 0.001)
  beta2 ~ dnorm( 0 , 0.001)
  beta3 ~ dnorm( 0 , 0.001)
  tau ~ dgamma( 1 , 0.001 )
  tau1 ~ dgamma( 1 , 0.001 )
  for (i in 1:L)
    weights[i]  <-  1
  s[1:k] ~ dcar_normal(adj[1:L],weights[1:L],num[1:k],tau,zero_mean=1)
  for (i in 1:k) {
    log(theta[i])  <-  alpha + beta1*x1[i] + beta2*x2[i] + beta3*x3[i] + s[i] + v[i] 
    y[i] ~ dpois( n[i]*theta[i] )
    v[i] ~ dnorm( 0 , tau1 )
  }  
})

n  = TX@data$population
y  = TX@data$stroke_count
x1 = TX@data$smoking
x2 = TX@data$physical
x3 = TX@data$drinking

TX.Consts = list( k=k , L=L , adj=adj , num=num , n=n , x1=x1, x2=x2, x3=x3 )
TX.Data = list( y = y )
TX.Inits = list( alpha = 0 , beta1 = 0 , beta2 = 0 , beta3 = 0 , tau = 1 , s=rep(0,k) , tau1 = 1 )

TX.Model = nimbleModel( TX.Code,
                        data = TX.Data,
                        constants = TX.Consts,
                        inits = TX.Inits )

TX.Model$initializeInfo()

compile.TX.Model = compileNimble( TX.Model )

TX.Conf = configureMCMC( TX.Model, print = TRUE )

TX.Conf$addMonitors(c("alpha","beta1","beta2","beta3","tau","tau1","theta"))

TX.MCMC = buildMCMC( TX.Conf )

compile.TX.MCMC = compileNimble( TX.MCMC, project = TX.Model )

niter = 110000
nburn =  100000
set.seed(1)

start.time = proc.time()

samples = runMCMC( compile.TX.MCMC, niter = niter, nburnin = nburn,
                   inits = TX.Inits, nchains = 1, samplesAsCodaMCMC = TRUE )

stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

head( samples )

StrokeEst = apply( samples[,7:260] , 2 , mean )
TX.sf$StrokeEst = StrokeEst*100000

# windows( 15 , 10 )
samples1 = samples
par(mfrow = c(4, 2), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
ts.plot(samples1[ , 'alpha'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(alpha), main = expression(alpha))
ts.plot(samples1[ , 'beta1'], xlab = 'iteration', col="darkgreen" , lwd=1.5 ,
        ylab = expression(beta), main = expression(beta[1]))
ts.plot(samples1[ , 'beta2'], xlab = 'iteration', col="darkgreen" , lwd=1.5 ,
        ylab = expression(beta), main = expression(beta[2]))
ts.plot(samples1[ , 'beta3'], xlab = 'iteration', col="darkgreen" , lwd=1.5 ,
        ylab = expression(beta), main = expression(beta[3]))
ts.plot(samples1[ , 'tau'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(tau), main = expression(tau))
ts.plot(samples1[ , 'tau1'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(tau[1]), main = expression(tau[1]))
ts.plot(samples1[ , 'theta[1]'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[1]), main = expression(theta[1]))
ts.plot(samples1[ , 'theta[2]'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[2]), main = expression(theta[2]))

# windows(6.5,6.5)
tm_shape( TX.sf ) +
  tm_polygons( col="StrokeEst" , border.col="black" , title="Stroke Mortality Rates (per 100,000)", breaks=seq(0,200,40)) 

#windows( 10 , 10 )
tauSample = samples[,"tau"]
tau1Sample = samples[,"tau1"]
par(mfrow = c(2, 1), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
plot( as.numeric(tauSample) , as.numeric(tau1Sample) , pch="." , log="xy" )
plot( as.numeric(tauSample) , as.numeric(tau1Sample) , pch="." )

windows (9,9)
brks = seq(-2,4,0.1)
par(mfrow = c(3, 1))
hist( samples[,"beta1"] , breaks=brks , xlab = expression(beta[1]), main = '') 
abline( v=0 , col='red')
hist( samples[,"beta2"] , breaks=brks , xlab = expression(beta[2]), main = '')
abline( v=0 , col='red')
hist( samples[,"beta3"] , breaks=brks , xlab = expression(beta[3]), main = '')
abline( v=0 , col='red')

setwd("C:/Users/chary/OneDrive/Documents/SLU/Sping 2020/BST 6200 spatial Statistics and Disease Mapping/R Scripts/")
source("DBDA2E-utilities.R")  ## Set of utilities written by John Kruschke ## 

# windows( 9 , 6 )
par(mfrow = c(1, 3))
plotPost( samples1[,"beta1"] , main=expression(beta[1]))
plotPost( samples1[,"beta2"] , main=expression(beta[2]))
plotPost( samples1[,"beta3"] , main=expression(beta[3]))
HDIofMCMC( samples1[,"beta1"] )
HDIofMCMC( samples1[,"beta2"] )
HDIofMCMC( samples1[,"beta3"] )
