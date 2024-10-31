

# Code for the manuscript: Managing the global wetland methane-climate feedback: A review of potential options
# Authors: Emily Ury, Eve-Lyn Hinckley, Daniele Visioni and Brian Buma
# Last update: October 31, 2024



library(terra)
library(sp)
library(rgdal)
library(raster)
library(bioclim)
library(geodata)
library(vtable)
library(ggplot2)

###

# Code to estimate fraction of feedbacks in quantiles
# parameters
se <- 0.2   #SD
fr <- 0.03  #fractional wetland pixels to remove

# iteration sample size
b <- c(1,2,3,5,7,9,11,16,21,26,31,40,50,60,75,99)
b <- 100-b
b <- rep(b, each=2)

################


#load data
# ssp 2.6
ssp <- 2.6
x2010 <- rast("Data/ch4_flux_annual_sum_2010_rcp26.tif")
x2099 <- rast("Data/ch4_flux_annual_sum_2099_rcp26.tif")

# ssp4.5
ssp <- 4.5
x2010 <- rast("Data/ch4_flux_annual_sum_2010_rcp45.tif")
x2099 <- rast("Data/ch4_flux_annual_sum_2099_rcp45.tif")

# Biomes
biomes <- rast("Data/Koppen_CZ.tif")

# load wetland raster, to get ride of infintesimal wetland areas
wetland <- rast("Data/Wetland_area_2020.tif")
wetland <- resample(wetland, x2010)
a.m <- cellSize(wetland, unit="m")			# cell size
a.km <- cellSize(wetland, unit="km")		# cell size
prop.wet <- wetland/a.km				        # proportion wetland
mask.n <- prop.wet > fr
mask.n <- mask.n*1
mask.n[mask.n == 0] <- NA
prop.wet[prop.wet>1] <- 1				#no cells >100% wetland, results from different data layers
par(mfrow=c(2,1)); plot(mask.n); plot(prop.wet)	#visualize

#checks
c10 <- x2010*a.m / 1e12
c99 <- x2099*a.m / 1e12
cellStats(raster(c10), sum, na.rm=T)		# working, produce same numbers as Zhang's publication.  Compare to masked (next section of code)
cellStats(raster(c99), sum, na.rm=T)		# for sensitivity test of removing low % wetland locations.

# these are the masked total emissions.  Compare to prior two lines of code for the unmasked to estimate bias associated with trimming small wetland areas
c10.m <- mask(c10, mask.n); cellStats(raster(c10.m), sum, na.rm=T)	# this is with the reduction of low CH4 areas, compare to test impact of reduction
c99.m <- mask(c99, mask.n); cellStats(raster(c99.m), sum, na.rm=T)	# this is with the reduction of low CH4 areas, compare to test impact of reduction


# setup iteration
index <- 1:length(b)

# setup storage
g <- list()					                     # for the global rasters
ch4 <- matrix(ncol=9, nrow=length(b))		 # for tg global, tg hotspot, area hotspot, percent (i)
nam <- c("Total feedback Tg","Quantile","Total feedback Tg norm",
         "Total in hotspot","Fraction Tg","Total area","Total area in hotspot",
         "Fraction area","SSP")
colnames(ch4) <- nam

bor.ras <- list()
boreal.store <- matrix(ncol=8, nrow=length(b))
nam <- c("Wetland area hotspot","Wetland area total","Fraction area hotspot","Methane hotspot",
         "Methane total","Fraction hotspot","ssp","percentile")
colnames(boreal.store) <- nam

temper.ras <- list()
temper.store <- matrix(ncol=8, nrow=length(b))
nam <- c("Wetland area hotspot","Wetland area total","Fraction area hotspot","Methane hotspot",
         "Methane total","Fraction hotspot","ssp","percentile")
colnames(temper.store) <- nam

trop.ras <- list()
trop.store <- matrix(ncol=8, nrow=length(b))
nam <- c("Wetland area hotspot","Wetland area total","Fraction area hotspot","Methane hotspot",
         "Methane total","Fraction hotspot","ssp","percentile")
colnames(trop.store) <- nam



for (i in index) {
  
  # add uncertainty factor
  semap <- rast(ext(c99), resolution=res(c99)); crs(semap) <- crs(c99)
  
  # reload
  c10 <- c10.m
  c99 <- c99.m
  
  # adjust emission map
  semap[] <- rnorm(ncell(semap), mean=1, sd=se)	# assign multiplyer
  c99 <- c99*semap
  semap[] <- rnorm(ncell(semap), mean=1, sd=se)	# assign multiplyer
  c10 <- c10*semap
  
  # calculate feedbacks
  feedback <- c99 - c10						    # simple difference
  feedback[feedback == 0] <- NA				
  feedback[feedback < 0] <- NA				# areas where the SE map pushes it below zero
  feedback <- mask(feedback,mask.n)		# mask out areas where wetlands are a small fraction of cell. 
  
  # store total feedback
  ch4[i,1] <- cellStats(raster(feedback), sum, na.rm=T)	# around XX extra Note Zhang et al. fig 1 is comparing to 2000)
  
  temp <- feedback
  
  # find hotspots
  brks <- quantile(values(temp, mat=F), probs=b[i]/100, na.rm=T)
  
  # calculate Tg in hotspot
  hot <- (temp > brks[[1]]); hot.i <- hot*1
  
  g[[i]] <- hot.i		# store the raster
  
  hot.i[hot.i == 0] <- NA
  temp.feedback <- mask(feedback, hot.i)
  temp.feedback <- cellStats(raster(temp.feedback), sum, na.rm=T)
  
  # how many Tg overall?
  temp.total <- feedback							# multply cells by m2 of wetlands in cells, bring down to Tg units
  temp.total <- cellStats(raster(temp.total), sum, na.rm=T)
  
  # fraction captured
  temp.frac <- temp.feedback/temp.total * 100
  
  # store methane sums
  ch4[i,3] <- temp.total
  ch4[i,4] <- temp.feedback
  ch4[i,5] <- temp.frac
  
  
  # store area
  t <- mask(wetland, hot.i)	#multiply wetland area raster by mask
  ch4[i,6] <- cellStats(raster(wetland), sum)
  ch4[i,7] <- cellStats(raster(t), sum)
  ch4[i,8] <- ch4[i,7]/ch4[i,6]
  
  ch4[i,2] <- b[i]
  
  ch4[i,9] <- ssp
  
  ###########
  # BIOMES
  biomes <- resample(biomes, feedback ,method="near")
  
  # make individual masks
  boreal.num <- cbind(1,1)
  temperate.num <- cbind(2,1)
  tropical.num <- cbind(3,1)
  
  boreal <- classify(biomes, boreal.num, others=NA)
  temperate <- classify(biomes, temperate.num, others = NA)
  tropical <- classify(biomes, tropical.num, others = NA)
  
  # trim to region
  boreal.norm <- mask(feedback,boreal)
  temperate.norm <- mask(feedback,temperate)
  tropical.norm <- mask(feedback,tropical)
  
  # find hotspots
  brks.boreal <- quantile(values(boreal.norm, mat=F), probs=c(b[i]/100), na.rm=T); brks.boreal
  brks.temperate <- quantile(values(temperate.norm, mat=F), probs=c(b[i]/100), na.rm=T); brks.temperate
  brks.tropical <- quantile(values(tropical.norm, mat=F), probs=c(b[i]/100), na.rm=T); brks.tropical
  
  #### make new objects for the high emission hotspots as 1, else as 0
  # biome specific
  bor <- boreal.norm > brks.boreal[[1]]
  bor <- bor * 1
  
  temper <- temperate.norm > brks.temperate[[1]]
  temper <- temper * 1
  
  trop <- tropical.norm > brks.tropical[[1]]
  trop <- trop * 1
  
  #  Boreal
  bor.feedback.all <- mask(feedback, boreal)
  bor.mask.hot <- bor; bor.mask.hot[bor.mask.hot==0] <- NA
  bor.feedback.hot <- mask(feedback,bor.mask.hot)
  
  t <- mask(wetland, bor.mask.hot)	#multiply wetland area raster by mask
  boreal.store[i,1] <- cellStats(raster(t), sum)
  t <- mask(wetland, boreal)	#multiply wetland area raster by mask
  boreal.store[i,2] <- cellStats(raster(t), sum)
  boreal.store[i,3] <- boreal.store[i,1]/boreal.store[i,2]
  
  #sum up methane for storage
  temp <- bor.feedback.hot 
  boreal.store[i,4] <- cellStats(raster(temp), sum)
  temp <- bor.feedback.all 
  boreal.store[i,5] <- cellStats(raster(temp), sum)
  boreal.store[i,6] <- boreal.store[i,4]/boreal.store[i,5] * 100
  
  boreal.store[i,7] <- ssp
  boreal.store[i,8] <- b[i]
  bor.ras[[i]] <- bor
  
  #  Temperate
  temper.feedback.all <- mask(feedback, temper)
  temper.mask.hot <- temper; temper.mask.hot[temper.mask.hot==0] <- NA
  temper.feedback.hot <- mask(feedback,temper.mask.hot)
  
  t <- mask(wetland, temper.mask.hot)	#multiply wetland area raster by mask
  temper.store[i,1] <- cellStats(raster(t), sum)
  t <- mask(wetland, temper)	#multiply wetland area raster by mask
  temper.store[i,2] <- cellStats(raster(t), sum)
  temper.store[i,3] <- temper.store[i,1]/temper.store[i,2]
  
  #sum up methane for storage
  temp <- temper.feedback.hot 		#multply cells by m2 of wetlands in cells, bring down to Tg units
  temper.store[i,4] <- cellStats(raster(temp), sum)
  temp <- temper.feedback.all 
  temper.store[i,5] <- cellStats(raster(temp), sum)
  temper.store[i,6] <- temper.store[i,4]/temper.store[i,5] * 100
  
  temper.store[i,7] <- ssp
  temper.store[i,8] <- b[i]
  temper.ras[[i]] <- temper
  
  #  Tropical
  trop.feedback.all <- mask(feedback, trop)
  trop.mask.hot <- trop; trop.mask.hot[trop.mask.hot==0] <- NA
  trop.feedback.hot <- mask(feedback,trop.mask.hot)
  
  t <- mask(wetland, trop.mask.hot)			#multiply wetland area raster by mask
  trop.store[i,1] <- cellStats(raster(t), sum)
  t <- mask(wetland, trop)				#multiply wetland area raster by mask
  trop.store[i,2] <- cellStats(raster(t), sum)
  trop.store[i,3] <- trop.store[i,1]/trop.store[i,2]
  
  #sum up methane for storage
  temp <- trop.feedback.hot 		#multply cells by m2 of wetlands in cells, bring down to Tg units
  trop.store[i,4] <- cellStats(raster(temp), sum)
  temp <- trop.feedback.all 
  trop.store[i,5] <- cellStats(raster(temp), sum)
  trop.store[i,6] <- trop.store[i,4]/trop.store[i,5] * 100
  
  trop.store[i,7] <- ssp
  trop.store[i,8] <- b[i]
  trop.ras[[i]] <- trop
  
  
  print(i/max(index))*100
}


# Organize outputs
# global
ch4 <- as.data.frame(ch4)

# summary 
means <- aggregate(ch4, list(Percent=ch4$"Quantile"), mean)
maxes <- aggregate(ch4, list(Percent=ch4$"Quantile"), max)
mins <- aggregate(ch4, list(Percent=ch4$"Quantile"), min)
means <- cbind(means, maxes$"Total in hotspot")
means <- cbind(means, mins$"Total in hotspot")

# temp storage, while running other scenario
tempstore.global <- means



# Boreal
boreal.store <- as.data.frame(boreal.store)
# summary (automated, clunky
means <- aggregate(boreal.store, list(Percent=boreal.store$"percentile"), mean)
maxes <- aggregate(boreal.store, list(Percent=boreal.store$"percentile"), max)
mins <- aggregate(boreal.store, list(Percent=boreal.store$"percentile"), min)

means <- cbind(means, maxes$"Methane hotspot")
means <- cbind(means, mins$"Methane hotspot")

# temp storage, while running other scenario
tempstore.boreal <- means

# Temperate
temper.store <- as.data.frame(temper.store)
means <- aggregate(temper.store, list(Percent=temper.store$"percentile"), mean)
maxes <- aggregate(temper.store, list(Percent=temper.store$"percentile"), max)
mins <- aggregate(temper.store, list(Percent=temper.store$"percentile"), min)

means <- cbind(means, maxes$"Methane hotspot")
means <- cbind(means, mins$"Methane hotspot")

# temp storage, while running other scenario
tempstore.temper <- means

# Tropical
trop.store <- as.data.frame(trop.store)
means <- aggregate(trop.store, list(Percent=trop.store$"percentile"), mean)
maxes <- aggregate(trop.store, list(Percent=trop.store$"percentile"), max)
mins <- aggregate(trop.store, list(Percent=trop.store$"percentile"), min)

means <- cbind(means, maxes$"Methane hotspot")
means <- cbind(means, mins$"Methane hotspot")

# temp storage, while running other scenario
tempstore.trop <- means

##########
##########
###stop###
##########
##########

# only run once. for ssp2.6
both.global <- tempstore.global
both.boreal <- tempstore.boreal
both.temper <- tempstore.temper
both.trop <- tempstore.trop

# run second time. combines ssp2.6 and ssp4.5
both.global <- rbind(both.global, tempstore.global)
both.boreal <- rbind(both.boreal, tempstore.boreal)
both.temper <- rbind(both.temper, tempstore.temper)
both.trop <- rbind(both.trop, tempstore.trop)

# add names
colnames(both.global) <- c("percent","total_methane_nonnorm","quantile","total_methane","hotspot_methane",
                           "fraction","total_area","hotspot_area","fraction_area","ssp","max_total_hotspot","min_total_hotspot")

nam <- c("percent","hotspot_area","total_area","fraction_area","hotspot_methane","total_methane",
         "fraction_methane","ssp","percentile","max_hotspot_methane","min_hotspot_methane")

colnames(both.boreal) <- nam
colnames(both.temper) <- nam
colnames(both.trop) <- nam

nrow(both.global)

ssp <- rep(2.6, 16)
ssp <- c(ssp, rep(4.5, 16))
both.global$ssp <- ssp
both.boreal$ssp <- ssp
both.temper$ssp <- ssp
both.trop$ssp <- ssp

##########
both.global[1:5,]

#### Supporting Info Fig S2

par(mfrow=c(1,2))

plot(both.global$hotspot_methane~both.global$fraction_area,bty="l", col=both.global$ssp*5,pch=16, ylim=c(0,90), 
     xlab="Fraction of global wetland area",ylab="Additional methane emissions (Tg/yr)", xlim=c(0,0.5))
arrows(x0=both.global$fraction_area, y0=both.global$max_total_hotspot, x1=both.global$fraction_area, y1=both.global$min_total_hotspot, 
       length=0.05, angle=90, code=3, col=both.global$ssp*5, lwd=1)
legend("topright",legend=c("Global","RCP 2.6","SSP 4.5"), 
       pch=c(16,NA,NA), 
       lwd=c(NA,2,2), col=c("black",2.6*5, 4.5*5), bg="white", bty="n")

# lines for scenarios in paper table

plot(both.boreal$hotspot_methane~both.boreal$fraction_area, bty="l",col=both.boreal$ssp*5, pch=1, ylim=c(0,90),
     xlab="Fraction of biome wetland area",ylab="", xlim=c(0,0.5))
arrows(x0=both.boreal$fraction_area, y0=both.boreal$max_hotspot_methane, x1=both.boreal$fraction_area, y1=both.boreal$min_hotspot_methane, 
       length=0.05, angle=90, code=3, col=both.boreal$ssp*5, lwd=1)

points(both.temper$hotspot_methane~both.temper$fraction_area, col=both.temper$ssp*5, pch=17)
arrows(x0=both.temper$fraction_area, y0=both.temper$max_hotspot_methane, x1=both.temper$fraction_area, y1=both.temper$min_hotspot_methane, 
       length=0.05, angle=90, code=3, col=both.temper$ssp*5, lwd=1)

points(both.trop$hotspot_methane~both.trop$fraction_area, col=both.trop$ssp*5, pch=16)
arrows(x0=both.trop$fraction_area, y0=both.trop$max_hotspot_methane, x1=both.trop$fraction_area, y1=both.trop$min_hotspot_methane, 
       length=0.05, angle=90, code=3, col=both.trop$ssp*5, lwd=1)

legend("topright",legend=c("Boreal","Temperate","Tropical","SSP 2.6","SSP 4.5"), 
       pch=c(1,17,16,NA,NA), 
       lwd=c(NA,NA,NA,2,2), col=c("black","black","black",2.6*5, 4.5*5), bg="white", bty="n")



























#Plotting
theme_set(theme_bw())

#total hotspot only
ggplot(data=both, aes(x=percent, y=total_hotspot, ymin=min_total_hotspot, ymax=max_total_hotspot,
                      fill=as.factor(ssp))) + 
  geom_line() + 
  geom_ribbon(alpha=0.5) + 
  ylab("Total in hotspot (Tg)") +
  xlab("Percent of global wetland area")

#fraction
ggplot(data=both, aes(x=fraction_area, y=total_hotspot/total_feedback_norm, ymin=max_total_hotspot/total_feedback_norm, ymax=min_total_hotspot/total_feedback_norm,
                      fill=as.factor(ssp))) + 
  geom_line() + 
  geom_ribbon(alpha=0.5) + 
  ylab("Fraction in hotspot (%)") +
  xlab("Percent of global wetland area")




