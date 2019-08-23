load(file="data/low density points.RData")
load("data/kernel_fits_trunc.RData")
load("data/Posteriors.RData")
plb<-read.csv("data/Plb_water_points6.csv")
dst<-read.csv("data/Distance2.csv")
rain.all<-read.csv("data/nrain_20year.csv")


################################################################################
#Add number of rainy days data to waterpoint dataset - Note, I should eventually merge this information into the original dataset so that this section isn't needed
rain <- rain.all[plb$r_raster3,] #Identify which rainfall cells each waterpoint falls on 
plb <- cbind(plb,rain[,4:23],dst) #Add number of rainy days and distance to nearest centre to end of waterpoint dataset
plb$Property <- as.character(plb$Property)
# fix a few other issues
plb$Purpose<-factor(plb$Purpose, levels = c(levels(plb$Purpose), "Unmanaged"))
plb$Property<-as.character(plb$Property)
plb$Property[is.na(plb$Property)]<-"Unclassified"
################################################################################