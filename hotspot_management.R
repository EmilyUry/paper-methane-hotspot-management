
# Code for the manuscript: Managing the global wetland methane-climate feedback: A review of potential options
# Authors: Emily Ury, Eve-Lyn Hinckley, Daniele Visioni and Brian Buma
# Last update: October 29, 2024



library(terra)
library(sp)
library(rgdal)
library(raster)
library(bioclim)
library(geodata)
library(vtable)

####################

# Ingest Koppen-Geiger climate zone maps 1991-2020
# Beck, et al. (2023). High-resolution (1 km) Köppen-Geiger maps for 1901-2099 based on constrained CMIP6 projections. Scientific Data, 10(1), 724. 
# These biomes have been regrouped into 3 climate regions based on Table S2 in the supplement
biomes <- rast("Data/Koppen_CZ.tif")
# 1 = Boreal
# 2 = Temperate
# 3 = Tropical


# Load climate data for current conditions and two future scenarios
# Data available from: https://www.worldclim.org/data/cmip6/cmip6climate.html
current.clim <- geodata::worldclim_global(var = 'bio', res = 10, download = F, path = "/Data")

# SSP 2.6
future.clim2.6 <- rast("Data/wc2.1_10m_bioc_HadGEM3-GC31-LL_ssp126_2081-2100.tif")
# SSP 4.5
future.clim4.5 <- rast("Data/wc2.1_10m_bioc_HadGEM3-GC31-LL_ssp245_2081-2100.tif")


## Load global wetland methane projection from 
# Zhang. et al. (2017). Emerging role of wetland methane emissions in driving 21st century climate change. Proceedings of the National Academy of Sciences of the United States of America, 114(36), 9647–9652. 
# Monthly methane fluxes have been summed to snapshots for the following years
# 1961 (baseline), 2010, 2099

x1961_2.6 <- rast("Data/ch4_flux_annual_sum_1961_rcp26.tif")
x2010_2.6 <- rast("Data/ch4_flux_annual_sum_2010_rcp26.tif")
x2099_2.6 <- rast("Data/ch4_flux_annual_sum_2099_rcp26.tif")

x1961_4.5 <- rast("Data/ch4_flux_annual_sum_1961_rcp45.tif")
x2010_4.5 <- rast("Data/ch4_flux_annual_sum_2010_rcp45.tif")
x2099_4.5 <- rast("Data/ch4_flux_annual_sum_2099_rcp45.tif")


######################################

# parameters
b <- 0.95				# hotspot quantile
fr <- 0.03			# fractional threshold removal of cells with minimal proportion wetland area (cells with <3% total wetland area)
n <- 100				# monte-carlo simulations
se <- 0.20			# uncertainty on emission map: generates a random adjustment with mean of 1 and SD of se


# load wetland raster from:
# Fluet-Chouinard et al.  (2023). Extensive global wetland loss over the past three centuries. Nature, 614(7947), 281-286.
wetland <- rast("Data/Wetland_area_2020.tif")

# Overlay wetland map with wetland methane emissions map
# remove very small wetland areas (causing mismatch)
wetland <- resample(wetland, x1961_4.5)		#to a methane emission raster, as the standard
a.m <- cellSize(wetland, unit="m")			#cell size in m
a.km <- cellSize(wetland, unit="km")		#cell size in km
prop.wet <- wetland/a.km				#proportion wetland
prop.wet[prop.wet > 1] <- 1
mask.n <- prop.wet > fr					#adjust for the fractional removal above
mask.n <- mask.n*1					#convert
mask.n[mask.n == 0] <- NA
par(mfrow=c(2,1)); plot(mask.n); plot(prop.wet)	#visualize 

#checks
c61 <- x1961_4.5*a.m / 1e12
c10 <- x2010_4.5*a.m / 1e12
c99 <- x2099_4.5*a.m / 1e12
cellStats(raster(c61), sum, na.rm=T)		#working, produce same numbers as Zhang's publication.  Compare to masked (next section of code)
cellStats(raster(c10), sum, na.rm=T)
cellStats(raster(c99), sum, na.rm=T)		#for sensitivity test of removing low % wetland locations.

#these are the masked total emissions.  Compare to prior two lines of code for the unmasked to estimate bias associated with trimming small wetland areas
c61.m <- mask(c61, mask.n); cellStats(raster(c61.m), sum, na.rm=T)
c10.m <- mask(c10, mask.n); cellStats(raster(c10.m), sum, na.rm=T)
c99.m <- mask(c99, mask.n); cellStats(raster(c99.m), sum, na.rm=T)

# reported difference in methane emissions using the 3% wetland area cutoff
# emissions in 1961 without cutoff:
round(cellStats(raster(c61), sum, na.rm=T),0)
# emisssions in 1961 with cutoff:
round(cellStats(raster(c61.m), sum, na.rm=T),0)
# percent difference:
round((cellStats(raster(c61), sum, na.rm=T) - cellStats(raster(c61.m), sum, na.rm=T))/cellStats(raster(c61), sum, na.rm=T)*100, 0)


###############################################

# spatial stability of hotspots including uncertainty


#store
g <- list()					#for the global rasters
bor.ras <- list()
temp.ras <- list()
trop.ras <- list()


ch4sum <- matrix(ncol=12, nrow=n)	#for the sum
colnames(ch4sum) <- c("Global Hotspot Tg","Global Total Tg", "Global Fraction in hotspot",
                      "Boreal Hotspot Tg","Boreal Total Tg", "Boreal Fraction in hotspot",
                      "Temperate Hotspot Tg","Temperate Total Tg", "Temperate Fraction in hotspot",
                      "Tropical Hotspot Tg","Tropical Total Tg", "Tropical Fraction in hotspot")

area.store <- matrix(ncol=12, nrow=n)
colnames(area.store) <- c("Global Hotspot area","Global Total area", "Global Fraction in hotspot",
                          "Boreal Hotspot area","Boreal Total area", "Boreal Fraction in hotspot",
                          "Temperate Hotspot area","Temperate Total area", "Temperate Fraction in hotspot",
                          "Tropical Hotspot area","Tropical Total area", "Tropical Fraction in hotspot")

store_sulfate_curr <- matrix(ncol=16, nrow=n)
colnames(store_sulfate_curr) <- c("Global total Tg","Global with S Tg","Global hotspot total Tg","Global hotspot with S Tg","Boreal total Tg","Boreal with S Tg",
                                  "Boreal hotspot total Tg,","Boreal hotspot with S Tg","Temperate total Tg","Temperate with S Tg","Temperate hotspot total Tg","Temperate hotspot with S Tg",
                                  "Tropical total Tg","Tropical with S Tg","Tropical hotspot total Tg","Tropical hotspot with S Tg")


store_sulfate_future <- matrix(ncol=16, nrow=n)
colnames(store_sulfate_future) <- c("Global total Tg","Global with S Tg","Global hotspot total Tg","Global hotspot with S Tg","Boreal total Tg","Boreal with S Tg",
                                    "Boreal hotspot total Tg,","Boreal hotspot with S Tg","Temperate total Tg","Temperate with S Tg","Temperate hotspot total Tg","Temperate hotspot with S Tg",
                                    "Tropical total Tg","Tropical with S Tg","Tropical hotspot total Tg","Tropical hotspot with S Tg")



#### CHANGE HERE FOR EACH SCENARIO
## Run for climate scenario 2.6
future.clim <- future.clim2.6
c61 <- x2010_2.6*a.m / 1e12		
c10 <- x2010_2.6*a.m / 1e12
c99 <- x2099_2.6*a.m / 1e12

## Run for climate scenario 4.5
# future.clim <- future.clim4.5
# c61 <- x2010_4.5*a.m / 1e12		
# c10 <- x2010_4.5*a.m / 1e12
# c99 <- x2099_4.5*a.m / 1e12

# initialize - total run takes 10-15 minutes
index <- 1:n
for (i in index) {
  
  # add uncertainty factor
  semap <- rast(ext(c99), resolution=res(c99)); crs(semap) <- crs(c99)
  
  # adjust emission map
  semap[] <- rnorm(ncell(semap), mean=1, sd=se)	#assign multiplyer
  c99t <- c99*semap
  semap[] <- rnorm(ncell(semap), mean=1, sd=se)	#assign multiplyer
  c61 <- c61*semap
  semap[] <- rnorm(ncell(semap), mean=1, sd=se)	#assign multiplyer
  c10 <- c10*semap
  
  # set cells that were randomly subtracted below zero to zero
  c99[c99<0] <- 0
  c61[c61<0] <- 0
  c10[c10<0] <- 0
  
  # calculate feedbacks
  feedback <- c99 - c61						                  # Simple difference
  feedback[feedback == 0] <- NA				              # remove zeros
  feedback <- mask(feedback,mask.n)		              # Mask out areas where wetlands are a small fraction of cell 
  feedback.norm <- feedback/wetland * 1e12/1e6			# wetland area normalized, g/m2
  
  # find hotspots
  brks <- quantile(values(feedback.norm, mat=F), probs=c(b), na.rm=T); brks
  
  # generate the masked raster for the global hotspots
  hot.95 <- (feedback.norm > brks[[1]]); hot.95 <- hot.95*1
  g[[i]] <- hot.95
  
  # sum emissions from hotspots (Tg)
  hot.95[hot.95 == 0] <- NA
  temp.feedback <- mask(feedback.norm, hot.95)
  temp.feedback <- temp.feedback * (wetland * 10^6) / 10^12          # multiply cells by m2 of wetlands in cells to get Tg
  temp.feedback <- cellStats(raster(temp.feedback), sum, na.rm=T)
  
  # total global emissions (Tg)
  temp.total <- feedback.norm * (wetland * 10^6) / 10^12             # multiply cells by m2 of wetlands in cells to get Tg
  temp.total <- cellStats(raster(temp.total), sum, na.rm=T)
  
  # fraction captured
  temp.frac <- temp.feedback/temp.total * 100
  
  # store methane sums
  ch4sum[i,1] <- temp.feedback
  ch4sum[i,2] <- temp.total
  ch4sum[i,3] <- temp.frac
  
  # store area sums
  t <- mask(wetland, hot.95)	                       # multiply wetland area raster by mask
  area.store[i,1] <- cellStats(raster(t), sum)
  area.store[i,2] <- cellStats(raster(wetland), sum)
  area.store[i,3] <- area.store[1,1]/area.store[1,2]
  
  # BIOMES
  biomes <- resample(biomes, feedback.norm,method="near")
  
  # make individual masks
  boreal.num <- cbind(1,1)
  temperate.num <- cbind(2,1)
  tropical.num <- cbind(3,1)
  
  boreal <- classify(biomes, boreal.num, others=NA)
  temperate <- classify(biomes, temperate.num, others = NA)
  tropical <- classify(biomes, tropical.num, others = NA)
  
  # subset feedback maps
  boreal.feedback <- mask(feedback.norm,boreal)
  temperate.feedback <- mask(feedback.norm,temperate)
  tropical.feedback <- mask(feedback.norm,tropical)
  
  # find hotspots
  brks.boreal <- quantile(values(boreal.feedback, mat=F), probs=c(b), na.rm=T); brks.boreal
  brks.temperate <- quantile(values(temperate.feedback, mat=F), probs=c(b), na.rm=T); brks.temperate
  brks.tropical <- quantile(values(tropical.feedback, mat=F), probs=c(b), na.rm=T); brks.tropical
  
  #### make new objects for the high emission hotspots as 1, else as 0
  # global only
  globe <- feedback.norm > brks[[1]] * 1
  globe <- globe * 1
  
  # biome specific
  bor <- boreal.feedback > brks.boreal[[1]]
  bor <- bor * 1
  
  temper <- temperate.feedback > brks.temperate[[1]]
  temper <- temper * 1
  
  trop <- tropical.feedback > brks.tropical[[1]]
  trop <- trop * 1
  
  #########################
  
  #  Calculate emissions from the hotspot - add back up to Tg/year
  #  Compare to overall emissions to find percentage
  
  #  Boreal
  bor.feedback.all <- mask(feedback.norm, boreal)
  bor.mask.hot <- bor; bor.mask.hot[bor.mask.hot==0] <- NA
  bor.feedback.hot <- mask(feedback.norm,bor.mask.hot)
  
  t <- mask(wetland, bor.mask.hot)	                   # multiply wetland area raster by mask
  area.store[i,4] <- cellStats(raster(t), sum)
  t <- mask(wetland, boreal)	                         # multiply wetland area raster by mask
  area.store[i,5] <- cellStats(raster(t), sum)
  area.store[i,6] <- area.store[1,4]/area.store[1,5]
  
  #sum up methane for storage
  temp <- bor.feedback.hot * (wetland * 10^6) / 10^12  # multply cells by m2 of wetlands in cells, bring down to Tg units
  ch4sum[i,4] <- cellStats(raster(temp), sum)
  temp <- bor.feedback.all * (wetland * 10^6) / 10^12
  ch4sum[i,5] <- cellStats(raster(temp), sum)
  ch4sum[i,6] <- ch4sum[i,4]/ch4sum[i,5] * 100
  
  bor.ras[[i]] <- bor
  
  
  #  Temperate
  temp.feedback.all <- mask(feedback.norm, temperate)
  temp.mask.hot <- temper; temp.mask.hot[temp.mask.hot==0] <- NA
  temp.feedback.hot <- mask(feedback.norm,temp.mask.hot)
  
  t <- mask(wetland, temp.mask.hot)	
  area.store[i,7] <- cellStats(raster(t), sum)
  t <- mask(wetland, temperate)	
  area.store[i,8] <- cellStats(raster(t), sum)
  area.store[i,9] <- area.store[1,7]/area.store[1,8]
  
  temp <- temp.feedback.hot * (wetland * 10^6) / 10^12  
  ch4sum[i,7] <- cellStats(raster(temp), sum)
  temp <- temp.feedback.all * (wetland * 10^6) / 10^12
  ch4sum[i,8] <- cellStats(raster(temp), sum)
  ch4sum[i,9] <- ch4sum[i,7]/ch4sum[i,8]* 100
  
  temp.ras[[i]] <- temper
  
  #  Tropical 
  trop.feedback.all <- mask(feedback.norm, tropical)
  trop.mask.hot <- trop; trop.mask.hot[trop.mask.hot==0] <- NA
  trop.feedback.hot <- mask(feedback.norm,trop.mask.hot)
  
  t <- mask(wetland, trop.mask.hot)	
  area.store[i,10] <- cellStats(raster(t), sum)
  t <- mask(wetland, tropical)	
  area.store[i,11] <- cellStats(raster(t), sum)
  area.store[i,12] <- area.store[1,10]/area.store[1,11]
  
  temp <- trop.feedback.hot * (wetland * 10^6) / 10^12  
  ch4sum[i,10] <- cellStats(raster(temp), sum)
  temp <- trop.feedback.all * (wetland * 10^6) / 10^12
  ch4sum[i,11] <- cellStats(raster(temp), sum)
  ch4sum[i,12] <- ch4sum[i,10]/ch4sum[i,11]* 100
  
  trop.ras[[i]] <- trop
  
  ####
  #  Characterizing hotspots - sulfate applications
  
  # upscale to methane resolution, using the boreal mask (doesn't matter which one)
  current.clim <- crop(current.clim, bor)
  current.clim <- resample(current.clim, bor, method="bilinear")
  future.clim <- crop(future.clim, bor)
  future.clim <- resample(future.clim, bor, method="bilinear")
  
  nam <- c("Annual mean temp","Mean diurnal range","Isothermality (Mean monthly max-min)","Temp seasonality (SD*100)","Max temp warmest month",
           "Min temp coldest month","Temp annual range","Mean temp, wettest quarter","Mean temp, driest quarter",
           "Mean temp, warmest quarter","Mean temp, coldest quarter","Annual precip","Precip, wettest month",
           "Precip, driest month","Precip seaonality","Precip, wettest quarter","Precip, driest quarter",
           "Precip, warmest quarter","Precip, coldest quarter")
  
  # pull out temperature CHANGE map 
  temp.change.map <- (future.clim[[1]] - current.clim[[1]])
  
  # Global sulfate deposition map (mg S m-2 yr-1)
  # Rubin et al. (2023). Global nitrogen and sulfur deposition mapping using a measurement–model fusion approach. Atmospheric Chemistry and Physics, 23(12), 7091–7102. 
  sulf <- rast("Data/Total_S_dep.tif")
  sulf <- resample(sulf, bor, method="bilinear")
  sulf <- sulf/1000000 * 10000

  # uncertainty parameters
  v <- 10.7
  k <- 10.6
  
  #calculations
  sulf.perc <- (rnorm(1, 38.6, v/2)*sulf)/(sulf+rnorm(1,8.71, k/2))		# assumes that the +/- in the paper is 2 SD's, per standard, but it's not explicit
  sulf.perc <- cellStats(raster(sulf.perc), max) - sulf.perc					 
  
  # Temp adjustment 
  # Per van Bodegom and Stams et al, Q10 study, methane increases 5.4x faster above 14C
  # assume a baseline of sulf.perc at 30% at 14C, getting relative decline based on the Q10 from van Bodegem
  
  dec <- function(x,y) {ifelse(x>14, -1.525*x + 51.35, y)}
  sulfate_potential_curr <- calc(raster(current.clim[[1]]), function(x){dec(x,30)})    # calculate percent reduction above 14
  sulfate_potential_curr <- sulfate_potential_curr - 30     				                   # converts to the reduction in effectiveness
  sulfate_potential_curr <- rast(sulfate_potential_curr)					                     # essentially the benefit is gone in hot areas
  
  adj.mask <- current.clim[[1]] > 14						     # assume above 14, no benefit
  adj.mask <- (adj.mask*1)
  adj.mask[adj.mask == 0] <- NA
  sulfate_potential_curr <- mask(sulfate_potential_curr, adj.mask)
  sulfate_potential_curr[is.na(sulfate_potential_curr[]) == T] <- 0
  sulfate_potential_curr <- sulf.perc + sulfate_potential_curr

  #same for future climate
  sulfate_potential_future <- calc(raster(future.clim[[1]]), function(x){dec(x,30)}) 	# calculate percent reduction above 14
  sulfate_potential_future<- sulfate_potential_future- 30     					              # converts to the reduction
  sulfate_potential_future<- rast(sulfate_potential_future)
  
  adj.mask <- future.clim[[1]] > 14
  adj.mask <- (adj.mask*1)
  adj.mask[adj.mask == 0] <- NA
  sulfate_potential_future<- mask(sulfate_potential_future, adj.mask)
  sulfate_potential_future[is.na(sulfate_potential_future[]) == T] <- 0
  
  sulfate_potential_future<- sulf.perc + sulfate_potential_future
  
  # convert to percentage
  sulfate_potential_curr <- sulfate_potential_curr/100
  sulfate_potential_future <- sulfate_potential_future/100
  
  # remove zeros
  sulfate_potential_curr[sulfate_potential_curr[] < 0] <- 0
  sulfate_potential_future[sulfate_potential_future[] < 0] <- 0
  
  # no evidence suppression can go higher than 40%, so cap that.
  sulfate_potential_curr[sulfate_potential_curr > .40] <- .40
  sulfate_potential_future[sulfate_potential_future> .40] <- .40
  
  # change in effectiveness - Figure for change over time, declining potential of this method
  sulfate_potential_change <- (sulfate_potential_curr-sulfate_potential_future)	  # boreal sees larger reduction in sulfate effectiveness

  ## Estimate application options
  ## Note that this is working off total emissions, not just additional feedback emissions
  
  ###    Global application of sulfate, 2010
  globe_sulfate_10 <- mask(c10, mask.n)
  change <- (globe_sulfate_10*sulfate_potential_curr)		                          # change by percentage
  change.sum <- (globe_sulfate_10 - change)			                                  # adjust estimates of original
  store_sulfate_curr[i,1] <- cellStats(raster(globe_sulfate_10), sum, na.rm=T)		# original emissions
  store_sulfate_curr[i,2] <- cellStats(raster(change.sum), sum, na.rm=T)			    # adjusted estimate
  
  # hotspot focus - global
  globe_sulfate_10.hot <- mask(globe_sulfate_10, hot.95)
  change <- (globe_sulfate_10.hot*sulfate_potential_curr)	
  change.sum <- (globe_sulfate_10.hot - change)	
  store_sulfate_curr[i,3] <- cellStats(raster(globe_sulfate_10.hot), sum, na.rm=T)		
  store_sulfate_curr[i,4] <- cellStats(raster(change.sum), sum, na.rm=T)		
  
  ###   boreal application of sulfate, 2010
  boreal_sulfate_10<- mask(c10, mask.n)
  boreal_sulfate_10[boreal_sulfate_10==0] <- NA
  boreal_sulfate_10<- mask(boreal_sulfate_10,  boreal)
  change <- (boreal_sulfate_10*sulfate_potential_curr)
  change.sum <- (boreal_sulfate_10- change)			
  store_sulfate_curr[i,5] <- cellStats(raster(boreal_sulfate_10), sum, na.rm=T)		# for sensitivity test of removing low % wetland locations.
  store_sulfate_curr[i,6] <- cellStats(raster(change.sum), sum, na.rm=T)			    # for sensitivity test of removing low % wetland locations.
  
  # hotspot focus - boreal
  boreal_sulfate_10.hot <- mask(globe_sulfate_10, bor.mask.hot)
  change <- (boreal_sulfate_10.hot*sulfate_potential_curr)	
  change.sum <- (boreal_sulfate_10.hot - change)		
  store_sulfate_curr[i,7] <- cellStats(raster(boreal_sulfate_10.hot), sum, na.rm=T)		
  store_sulfate_curr[i,8] <- cellStats(raster(change.sum), sum, na.rm=T)		
  
  ###   temperate application of sulfate, 2010
  temp_sulfate_10 <- mask(c10, mask.n)
  temp_sulfate_10[temp_sulfate_10==0] <- NA
  temp_sulfate_10<- mask(temp_sulfate_10,  temper)	                            
  change <- (temp_sulfate_10*sulfate_potential_curr)	                          
  change.sum <- (temp_sulfate_10- change)			                                  
  store_sulfate_curr[i,9] <- cellStats(raster(temp_sulfate_10), sum, na.rm=T)		
  store_sulfate_curr[i,10] <- cellStats(raster(change.sum), sum, na.rm=T)			  
  
  # hotspot focus - temperate 
  temp_sulfate_10.hot <- mask(temp_sulfate_10, temp.mask.hot)
  change <- (temp_sulfate_10.hot*sulfate_potential_curr)	
  change.sum <- (temp_sulfate_10.hot - change)		
  store_sulfate_curr[i,11] <- cellStats(raster(temp_sulfate_10.hot), sum, na.rm=T)		
  store_sulfate_curr[i,12] <- cellStats(raster(change.sum), sum, na.rm=T)		
  
  ###   tropical application of sulfate, 2061
  trop_sulfate_10<- mask(c10, mask.n)
  trop_sulfate_10[trop_sulfate_10==0] <- NA
  trop_sulfate_10<- mask(trop_sulfate_10,  trop)			
  change <- (trop_sulfate_10*sulfate_potential_curr)		
  change.sum <- (trop_sulfate_10- change)				
  store_sulfate_curr[i,13] <- cellStats(raster(trop_sulfate_10), sum, na.rm=T)		
  store_sulfate_curr[i,14] <- cellStats(raster(change.sum), sum, na.rm=T)			
  
  # hotspot focus - tropical
  trop_sulfate_10.hot <- mask(trop_sulfate_10, trop.mask.hot)
  change <- (trop_sulfate_10.hot*sulfate_potential_curr)		
  change.sum <- (trop_sulfate_10.hot - change)			
  store_sulfate_curr[i,15] <- cellStats(raster(trop_sulfate_10.hot), sum, na.rm=T)		
  store_sulfate_curr[i,16] <- cellStats(raster(change.sum), sum, na.rm=T)			
  

  #################
  ###    Global application of sulfate, 2099
  globe_sulfate_99 <- mask(c99, mask.n)
  globe_sulfate_99[globe_sulfate_99 ==0] <- NA
  change <- (globe_sulfate_99*sulfate_potential_future)	
  change.sum <- (globe_sulfate_99 - change)			
  store_sulfate_future[i,1] <- cellStats(raster(globe_sulfate_99), sum, na.rm=T)		
  store_sulfate_future[i,2] <- cellStats(raster(change.sum), sum, na.rm=T)	
  
  # hotspot focus - global
  globe_sulfate_99.hot <- mask(globe_sulfate_99, hot.95)
  change <- (globe_sulfate_99.hot*sulfate_potential_future)
  change.sum <- (globe_sulfate_99.hot - change)		
  store_sulfate_future[i,3] <- cellStats(raster(globe_sulfate_99.hot), sum, na.rm=T)		
  store_sulfate_future[i,4] <- cellStats(raster(change.sum), sum, na.rm=T)		
  
  ###   boreal application of sulfate, 2099
  boreal_sulfate_99 <- mask(c99, mask.n)						                            # new copy of original 2099 map
  boreal_sulfate_99[boreal_sulfate_99 ==0] <- NA
  boreal_sulfate_99 <- mask(boreal_sulfate_99,  boreal)	                        # get to just boreal focus
  change <- (boreal_sulfate_99*sulfate_potential_future)	                      # change by percentage
  change.sum <- (boreal_sulfate_99 - change)			                              # adjust estimates of original
  store_sulfate_future[i,5] <- cellStats(raster(boreal_sulfate_99), sum, na.rm=T)		# for sensitivity test of removing low % wetland locations.
  store_sulfate_future[i,6] <- cellStats(raster(change.sum), sum, na.rm=T)			# for sensitivity test of removing low % wetland locations.
  
  # hotspot focus - boreal
  boreal_sulfate_99.hot <- mask(globe_sulfate_99, bor.mask.hot)
  change <- (boreal_sulfate_99.hot*sulfate_potential_future)
  change.sum <- (boreal_sulfate_99.hot - change)		
  store_sulfate_future[i,7] <- cellStats(raster(boreal_sulfate_99.hot), sum, na.rm=T)		
  store_sulfate_future[i,8] <- cellStats(raster(change.sum), sum, na.rm=T)		
  
  ### temperate zone application of sulfate, 2099
  temp_sulfate_99 <- mask(c99, mask.n)
  temp_sulfate_99[temp_sulfate_99 ==0] <- NA
  temp_sulfate_99 <- mask(temp_sulfate_99,  temper)	                            # get to just temperate zone focus
  change <- (temp_sulfate_99*sulfate_potential_future)	
  change.sum <- (temp_sulfate_99 - change)		
  store_sulfate_future[i,9] <- cellStats(raster(temp_sulfate_99), sum, na.rm=T)	
  store_sulfate_future[i,10] <- cellStats(raster(change.sum), sum, na.rm=T)			
  
  # hotspot focus - temperate 
  temp_sulfate_99.hot <- mask(temp_sulfate_99, temp.mask.hot)
  change <- (temp_sulfate_99.hot*sulfate_potential_future)	
  change.sum <- (temp_sulfate_99.hot - change)			
  store_sulfate_future[i,11] <- cellStats(raster(temp_sulfate_99.hot), sum, na.rm=T)		
  store_sulfate_future[i,12] <- cellStats(raster(change.sum), sum, na.rm=T)	
  
  ###   tropical application of sulfate, 2099
  trop_sulfate_99 <- mask(c99, mask.n)
  trop_sulfate_99[trop_sulfate_99 ==0] <- NA
  trop_sulfate_99 <- mask(trop_sulfate_99,  trop)	
  change <- (trop_sulfate_99*sulfate_potential_future)
  change.sum <- (trop_sulfate_99 - change)			
  store_sulfate_future[i,13] <- cellStats(raster(trop_sulfate_99), sum, na.rm=T)		
  store_sulfate_future[i,14] <- cellStats(raster(change.sum), sum, na.rm=T)			
  
  # hotspot focus - tropical
  trop_sulfate_99.hot <- mask(trop_sulfate_99, trop.mask.hot)
  change <- (trop_sulfate_99.hot*sulfate_potential_future)	
  change.sum <- (trop_sulfate_99.hot - change)			
  store_sulfate_future[i,15] <- cellStats(raster(trop_sulfate_99.hot), sum, na.rm=T)		
  store_sulfate_future[i,16] <- cellStats(raster(change.sum), sum, na.rm=T)		
  
  
  print(i/max(index)*100)
}


# Data outputs
area.store			          	# size of hotspots and fraction of total area
ch4sum				              # size of emissions (Tg/year) from hotspots and total area
store_sulfate_curr		      # sulfate effect size
store_sulfate_future		    # future sulfate effect size


### Data reported in table 2
st(as.data.frame(area.store))
st(as.data.frame(ch4sum))


# Sulfate post-process
temp.sulf <- store_sulfate_curr

store <- matrix(nrow=nrow(temp.sulf), ncol=8)
store[,1] <- temp.sulf[,1] - temp.sulf[,2]
store[,2] <- temp.sulf[,3] - temp.sulf[,4]
store[,3] <- temp.sulf[,5] - temp.sulf[,6]
store[,4] <- temp.sulf[,7] - temp.sulf[,8]
store[,5] <- temp.sulf[,9] - temp.sulf[,10]
store[,6] <- temp.sulf[,11] - temp.sulf[,12]
store[,7] <- temp.sulf[,13] - temp.sulf[,14]
store[,8] <- temp.sulf[,15] - temp.sulf[,16]

colnames(store) <- c("Global sulfate reduction", "Global hotspot sulfate reduction", "Boreal sulfate reduction",
                     "Boreal hotspot sulfate reduction","Temperate sulfate reduction","Temperate hotspot sulfate reduction",
                     "Tropical sulfate reduction","Tropical hotspot sulfate reduction")

## Data reported in table 3 - 2010
st(as.data.frame(store))
st(as.data.frame(temp.sulf))

################################################

temp.sulf <- store_sulfate_future

store <- matrix(nrow=nrow(temp.sulf), ncol=8)
store[,1] <- temp.sulf[,1] - temp.sulf[,2]
store[,2] <- temp.sulf[,3] - temp.sulf[,4]
store[,3] <- temp.sulf[,5] - temp.sulf[,6]
store[,4] <- temp.sulf[,7] - temp.sulf[,8]
store[,5] <- temp.sulf[,9] - temp.sulf[,10]
store[,6] <- temp.sulf[,11] - temp.sulf[,12]
store[,7] <- temp.sulf[,13] - temp.sulf[,14]
store[,8] <- temp.sulf[,15] - temp.sulf[,16]

colnames(store) <- c("Global sulfate reduction", "Global hotspot sulfate reduction", "Boreal sulfate reduction",
                     "Boreal hotspot sulfate reduction","Temperate sulfate reduction","Temperate hotspot sulfate reduction",
                     "Tropical sulfate reduction","Tropical hotspot sulfate reduction")

## Data reported in Table 3 - 2099
st(as.data.frame(store))
st(as.data.frame(temp.sulf))



###############################


# Figures for the paper:
par(mfrow=c(2,1))

# Plot global feedback methane emissions
temp <- x2099_4.5 - x1961_4.5; temp[temp == 0] <- NA

library(RColorBrewer)

## Figure 1
tiff(filename = "Figure1.tiff", height = 7.5, width = 6.75, units = "in", res = 800, compression = "lzw")

{
par(mfrow = c(2,1)) 
colors <- c(colorRampPalette(c("gray95","gold2", "darkorange4"))(50))
  
plot(temp, main=" ",
     col=colors)
title("(a) Feedback intensity (2099-1961 rates; g CH4/m2 per yr), RCP 4.5", adj = 0, cex.main = 1)

#subsets
e1 <- extent(c(-82.55, -59.72, -19.48, 3.63))  #NW S America
e1.crop <- crop(temp, e1)
e2 <- extent(c(75.67, 130.21, 7.66, 32.91))    #SE Asia
e2.crop <- crop(temp, e2)
e3 <- extent(c(-94.1, -72.3, 22.4, 58.4))
e3.crop <- crop(temp, e3)			 #E N America

e <- ext(c(-160,-100,-65,0))
inset(e1.crop, e=e, col=colors)
lines(e1)
e <- ext(c(10,110,-85,-35))
inset(e2.crop, e=e, col=colors)
lines(e2)
e <- ext(c(-55,-20,10,85))
inset(e3.crop, e=e, col=colors)
lines(e3)

gstack <- rast(c(g))
gstack.mean <- app(gstack, "mean")

#subsets
e1 <- extent(c(-82.55, -59.72, -19.48, 3.63))  #NW S America
e1.crop <- crop(gstack.mean, e1)
e2 <- extent(c(75.67, 130.21, 7.66, 32.91))    #SE Asia
e2.crop <- crop(gstack.mean, e2)
e3 <- extent(c(-94.1, -72.3, 22.4, 58.4))
e3.crop <- crop(gstack.mean, e3)			 #E N America

colors <- c(colorRampPalette(c("gray95","lightgoldenrod", "midnightblue"))(50))

plot(gstack.mean, main=" ", adj = 0, col=colors)
title("(b) Hotspots (>95th percentile), RCP 4.5", adj = 0, cex.main = 1)

e <- ext(c(-160,-100,-65,0))
inset(e1.crop, e=e, col=colors)
lines(e1)
e <- ext(c(20,110,-85,-40))
inset(e2.crop, e=e, col=colors)
lines(e2)
e <- ext(c(-55,-20,10,85))
inset(e3.crop, e=e, col=colors)
lines(e3)

}

dev.off()

##########################################

## graphical abstract
library(sf)
library(ggplot2)
coasts <- st_read("Data/Coastline/ne_110m_coastline.shp", quiet = TRUE)
map <- results <- as.data.frame(gstack.mean, xy = TRUE)

tiff(filename = "GraphicalAbstract.tiff", height = 2, width = 4, units = "in", res = 800, compression = "lzw")

ggplot() +
  geom_tile(data = map, aes(x=x, y = y, fill = mean)) +
  scale_fill_gradient(low = "gray90", high = "midnightblue", na.value = "#00000000") +
  geom_sf(data = coasts, color = "gray30", fill = NA,linewidth = 0.1) +
  theme_void(base_size = 10) +
  xlab(" ") +
  ylab(" ") +
  ylim(-50,100) +
  xlim(-150, 170) +
  theme(legend.position = "none", 
        plot.title = element_text(size = 10, hjust = 0.5, face = "plain"), 
        axis.title = element_text(size = 10, face = 'plain'))
  
dev.off()

























