############################################################################################################################
#Cost-effective management of cane toads - last edited 18/11/14
###########################################################################################################################

#Water point data is different to Tingley et al. (2013) -> updated from field visit
#This code divides the corridor into chunks and places the centre of a single barrier on each interval. 
#At each location, the code cycles through a range of fixed budgets. 
#Natural, artificial or irrigated areas can be managed. Future irrigation can be included or excluded
#There is an option to truncate the dispersal kernel, although it does not redistribute the cut off area to the shoulder of the kernel
#The code cycles through a 20 year time series of rainfall data, rather than using averages
#Irrigated areas plus La Grange have a unique detection radius and emit a unique number of propogules
#The cost model has been further developed from previous versions to include upfront (including labor, distance from major centres),
#maintenance and replacement costs with a discounting rate of 5% per annum for future management

#HPC script
args=(commandArgs(TRUE))

for(i in 1:length(args)) {
   eval(parse(text=args[[i]]))
}

#Specify if running on JCU computer or laptop
JCU <- TRUE
#Read in required files from either JCU or laptop
if (JCU==TRUE) {
setwd("~/Results_10Dec")
source("~/code/19_pprocess_10Dec.R")
load(file="~/Data/low density points.RData")
load("~/Data/kernel_fits_trunc.RData")
load("~/Data/Posteriors.RData")
plb<-read.csv("~/Data/Plb_water_points6.csv")
dst<-read.csv("~/Data/Distance2.csv")
rain.all<-read.csv("~/Data/nrain_20year.csv")
} else {
setwd("C:/Darren/1_PhD/Project proposal/11_Cane toads/Rcode/evo-dispersal-read-only/evo-dispersal-read-only/temp")
source("C:/Darren/1_PhD/Project proposal/11_Cane toads/13_FINAL_versions_data_drafts/Running_the_model/19_pprocess_10Dec.R")
load(file="C:/Darren/1_PhD/Project proposal/11_Cane toads/Rcode/evo-dispersal-read-only/evo-dispersal-read-only/Data/low density points.RData")
load("C:/Darren/1_PhD/Project proposal/11_Cane toads/Rcode/evo-dispersal-read-only/evo-dispersal-read-only/Data/kernel_fits_trunc.RData")
load("C:/Darren/1_PhD/Project proposal/11_Cane toads/Rcode/evo-dispersal-read-only/evo-dispersal-read-only/Data/Posteriors.RData")
plb<-read.csv("C:/Darren/1_PhD/Project proposal/11_Cane toads/Rcode/evo-dispersal-read-only/evo-dispersal-read-only/Data/Plb_water_points6.csv")
dst<-read.csv("C:/Darren/1_PhD/Project proposal/11_Cane toads/Rcode/evo-dispersal-read-only/evo-dispersal-read-only/Data/Distance2.csv")
rain.all<-read.csv("C:/Darren/1_PhD/Project proposal/11_Cane toads/Rcode/evo-dispersal-read-only/evo-dispersal-read-only/Data/nrain_20year.csv", header = TRUE)
}

################################################################################
#Add number of rainy days data to waterpoint dataset - Note, I should eventually merge this information into the original dataset so that this section isn't needed
rain <- rain.all[plb$r_raster3,] #Identify which rainfall cells each waterpoint falls on 
plb <- cbind(plb,rain[,4:23],dst) #Add number of rainy days and distance to nearest centre to end of waterpoint dataset

################################################################################
#Define how to run the model

future.irrigation <- FALSE #Set to TRUE to include future irrigation, FALSE otherwise
current.irrigation <- TRUE #Set to TRUE to include current irrigation, FALSE otherwise
IPA <- TRUE #Set to TRUE to run with IPA active, FALSE otherwise
Properties <- FALSE #Set to TRUE to remove waterpoints from a specific property, FALSE otherwise. Need to specify property name below
Name = as.character("Thangoo") #Doesn't work yet
Dwelling <- TRUE #Set to TRUE to include dwellings in simulation, FALSE to remove
trunc.kernel <- TRUE #Set to TRUE to run the model with a truncated kernel, FALSE to use the original model from Tingley et al. (2013)
trunc.dist <- 55000 #Set the truncation distance in metres. Choose between 30000, 55000 and 78000
scenario <- 1 #Set to 1 to include rainfall variability, 2 to fix ndays to 180 and 3 to set ndays to 50 
failure.rate <- 0.95 #Set failure rate of managed waterbodies (i.e. the probability that a managed water body becomes available for colonisation
plot.spread <- TRUE #Set to TRUE for plots to be drawn during simulations, FALSE otherwise
reps <- 10 #Set number of simulations for each combination of budget and location

#Set management decisions
gens <- 50 #Set length of management horizon
barriers <- 21 #Set the number of bins in the landscape. Note there one less barrier (i.e. if barrier is set to 3, there will be 2 barrier locations
rmnat <- FALSE #Set to TRUE if natural waterbodies are managed, FALSE otherwise
rmirr <- FALSE #Set to TRUE if irrigated points are managed, FALSE otherwise
rmdwel <- FALSE #Always set to FALSE because we assume dwellings can't be managed
budget <- seq(0,6000000,200000) #Set the budget to spent on at each barrier location

###############################################################################
#Modify waterbody dataset based on above management decisions
plb <- modify.plb(plb, future.irrigation, current.irrigation, IPA, Properties, Name, Dwelling)

pairs_pdist<-pdist.fast(X=plb$POINT_X,Y=plb$POINT_Y,maximum=500000,space.size= 500000) #Calculates distance between all combinations of water bodies

###################################################################################
#Extract relevant data from csv file
ID <- as.numeric(rownames(plb))
X <- plb$POINT_X #X coordinates of waterpoints
Y <- plb$POINT_Y #Y coordinates of waterpoints
Pres <- plb$ARRIVE_MCP #Presence/absence of toads at time=0, 0 = absent, 1 = present, 2 = targets
age <- rep(0,length(X)) #Track number of years each waterpoint is colonised SHOULD THIS BE EQUAL TO PRESENT?
Man <- rep(0,length(X)) #Identify which waterbodies are managed
nats <- which(plb$art_nat==0) #which water bodies (IDs) are natural 
arts <- which(plb$art_nat==1) #which water bodies (IDs) are artificial
irr <- which(plb$Name=="Irrigation") #which water bodies (IDs) are irrigation
dwell <- which(plb$Old_name =="Dwelling") #which water bodies (IDs) are dwellings

###################################################################################
#Create 100 years of rainfall data for each site by looping across 20 yeras 5 times
rain <- plb[,17:36]
rain <- cbind(rain,rain,rain,rain,rain)

#Identify where to place the centre of a barrier in the landscape
corridor <- corridor.smoother(X=X,Y=Y,plb=plb,interval=barriers,plot.barriers=TRUE)

#Calculate how much it costs to manage each waterpoint for the duration of the management program
cost.WP <- cost.waterpoint2(gens,plb=plb,dst=dst,plot.cost=TRUE)
plb <- cbind(plb,cost.WP) #Add cost to waterpoint dataset

# calculate n.pairs using pdist
n.pairs<-do.call("c",lapply(pairs_pdist,nrow)) #Finds how many pairs there are for each point

#Set up matrix to store u and v for each waterpoint
u <- matrix(NA,nrow(plb),3)
colnames(u) <- c("u","v","trc")

#Create table to store results
spread.table<-as.matrix(cbind(ID,X,Y,Pres,n.pairs,u,age,Man),nrow=length(age),ncol=8)

###################################################################################
#Set detection radii for current irrigation plus the La Grange
r.samp <- detection.radius(plb)

##################################################################################
#Run simulations
output <- vector("list", length=length(budget)*nrow(corridor)*reps) # length = no. of points Manoved x no. of locations where points are Manoved x no. of simulations

ptm <- proc.time() #record how long it takes to run simulations

for (kk in 1:length(budget)){ #loop through range of budgets
	for(ii in 2:(nrow(corridor)-3)){ #loop through each corridor location excluding first and last points
		mod <- knock.out.nn.budget(X=corridor[ii,"X"], Y=corridor[ii,"Y"], spread.table=spread.table, n=budget[kk], natural=nats, irrigation=irr, Pres, dwelling=dwell, cost.WP)
		target <- which(mod$spread.table[,"Pres"]==2) #Identify which waterpoints terminate simulations
		mod$spread.table[mod$spread.table[,"Pres"]==2,"Pres"] <- 0 #Make sure they're vacant at start of simulations
		for (jj in 1:reps){ 
			lambda.samp <- 10^rnorm(1, mean=sample.lambda, sd=sample.lambda.sd) #mean number of propagules dispersing from a point -> becomes delta
			managed <- mod$spread.table[,"Man"]
			temp <- spread.pilb.failure(pop=mod$spread.table, gens=gens, pairs=mod$pairs.mod, target=target, delta=lambda.samp, r=r.samp, trunc.dist=trunc.dist) #spread function
			occ <- sum(temp$popmatrix[,"age"] > 0) #record number of waterpoints colonised
			nn <- sum(temp$popmatrix[,"Man"] > 0) #record number of waterpoints managed 
			temp <- c(temp, list(pars=cbind(lambda=lambda.samp, r=r.samp)), list(Removed=nn), list(location=ii), list(colonised=occ), list(cost=budget[kk]))
			output[[nrow(corridor)*reps*(kk-1)+reps*(ii-1)+jj]] <- temp
		print(jj)
		}
	}
}

proc.time() - ptm #Stop timer

save(output, file=paste("IPAT_T55_S1_IF_DF",File.ID,".RData", sep=""))


