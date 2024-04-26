# Pulls Darren's dataset together.  
  # Note: 
    # Spatial location data in Plb_water_points6.csv largely does not match that in art_nat dataset
  
source("modelData/19_pprocess_10Dec.R")
plb<-read.csv("modelData/Plb_water_points6.csv")
dst<-read.csv("modelData/Distance2.csv")
rain.all<-read.csv("modelData/nrain_20year.csv", header = TRUE)

#Add number of rainy days data to waterpoint dataset
rain <- rain.all[plb$r_raster3,] #Identify which rainfall cells each waterpoint falls on 
plb <- cbind(plb,rain[,4:23],dst) #Add number of rainy days and distance to nearest centre to end of waterpoint dataset
plb$Property <- as.character(plb$Property)

# Append costs to dataframe
# modified cost function -- returns conversion, and annual maintenance costs
#New cost model with two types of irrigation
cost.waterpoint2 <- function(plb) { 
  cost.WP <- matrix(NA,nrow(plb),2) # matrix for conversion and annual maintenance cost
  purpose <- plb$Purpose
  type <- plb$Type
  Property <- plb$Property
  Area <- plb$Area
  perimeter <- plb$Perimeter/100 #Perimeter of fences in 100s of metres
  disc<-0 # calculate discounting at time=0 only (i.e. maintenance, per annum)
  for (cc in 1:nrow(plb)){
    if (purpose[cc] == "Artificial" & type[cc] == 1) #Tanks
    {tank.installation <- plb[cc,37]*1.50 + 32*100 + 8500 #Installation cost = Distance (nearest town * $1.50) + Labour (2 men, 2 days @ $100/hr) + Materials ($8500)
    tank.maintenance <- (plb[cc,37]*1.50 + 8*100 + 500)*(1+0.025)^disc #Maintenance = Distance (nearest town * $1.50) + Labour ( 1 man, 1 day @ $100/hr) + Materials ($500) with a discount of 5% p.a 
    tank.repair <- tank.installation/50*(1+0.025)^disc #Repair = Installation cost / 50 subject to discount rate
    cost.WP[cc, 1] <- tank.installation
    cost.WP[cc, 2] <- tank.maintenance + tank.repair} #Installation cost (Distance + Labour + Materials)
    
    if (purpose[cc] == "Artificial" & type[cc] == 2) #Dams
    {tank.installation <- plb[cc,37]*1.50 + 32*100 + 8500 #Installation cost = Distance (nearest town * $1.50) + Labour (2 men, 2 days @ $100/hr) + Materials ($8500)
    tank.maintenance <- (plb[cc,37]*1.50 + 8*100 + 500)*(1+0.025)^disc #Maintenance = Distance (nearest town * $1.50) + Labour ( 1 man, 1 day @ $100/hr) + Materials ($500) with a discount of 5% p.a 
    tank.repair <- tank.installation/50*(1+0.025)^disc #Repair = Installation cost / 50 subject to discount rate
    cost.WP[cc, 1] <- tank.installation
    cost.WP[cc, 2] <- tank.maintenance + tank.repair} #Installation cost (Distance + Labour + Materials)
    
    if (purpose[cc] == "Irrigation" & !Property[cc] == "Shalemar") #Irrigation - Pay for hay production to shut down for a period
    {cost.WP[cc, 1] <- 0
     cost.WP[cc, 2] <- Area[cc]*100*30*65} #Convert km2 to ha and multiply by profit gained for hay per ha 
    
    if (purpose[cc] == "Irrigation" & Property[cc] == "Shalemar") #Irrigation - Fence year round irrigation
    {fence.installation <- plb[cc,37]*1.50 + 32*100*perimeter[cc] + 3000*perimeter[cc] #Installation = Distance (nearest town * $1.50) + Labour (Perimeter * time per unit @ $100/hr) + Materials (Cost per unit length * Perimeter)
    fence.maintenance <- (8*100*26+100)*(1+0.025)^disc #Maintenance = Labour (1 day per week @ $100 per/h for 6 months) + cost ($100)
    fence.repair <- fence.installation/10*(1+0.025)^disc #Repair = Installation cost / 50 subject to discount rate
    cost.WP[cc, 1] <- fence.installation
    cost.WP[cc, 2] <- fence.maintenance + fence.repair}
    if (purpose[cc] == "Dwelling") {cost.WP[cc,] <- 0} #Dwelling
    if (purpose[cc] == "Natural") {cost.WP[cc,] <- 0} #Don't manage natural waterpoints
  }
 
  return(cost.WP)	
}

cMat<-cost.waterpoint2(plb=plb)
colnames(cMat)<-c("Conversion", "Maintenance")
plb<-cbind(plb, cMat)


write.csv(plb, file="mapData/waterbodies.csv") # To be converted from Oz Albers to WGS offline in QGIS..
