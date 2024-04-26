############ Removes rows from the waterpoints dataset based on management decisions
modify.plb <- function(plb, future.irrigation, IPA, Dwelling, Properties){ 
    #Future irrigation
    if (future.irrigation == FALSE) {
      plb <- subset(plb,!(plb$Purpose == "Irrigation" & plb$Type == 2))
    }
    #IPA areas
    if (IPA == TRUE) {
      plb <- subset(plb,!(plb$Property == "Frazier Downs" & plb$IPA==0))
    }
    #Properties
    if (!is.na(Properties)) {
      plb <- subset(plb,!plb$Property %in% Properties)
    }
    #Dwellings
    if (Dwelling == FALSE) {
      plb <- subset(plb,!plb$Purpose == "Dwelling")
    }
    plb$FID <- c(1:nrow(plb))
    rownames(plb) <- plb$FID 
    return(plb)
  }

#Sets unique detection radii for irrigation and La Grange
detection.radius <- function(plb) {
	r.samp<-rep(10^2,nrow(plb)) #Detection radius -> becomes r
	area<-sqrt((plb$Area*1000000)/pi) #Convert area from square kilometers to square meters then calculate radius of a circle of the same area
	r.samp[which(plb$Status == "Current" & plb$Property == "Shamrock Gardens")] <- area[which(plb$Status == "Current" & plb$Property == "Shamrock Gardens")]
	r.samp[which(plb$Status == "Current" & plb$Property == "Shalemar")] <- area[which(plb$Status == "Current" & plb$Property == "Shalemar")]
	r.samp[which(plb$Status == "Current" & plb$Property == "Pardoo")] <- area[which(plb$Status == "Current" & plb$Property == "Pardoo")]
	r.samp[which(plb$Status == "Current" & plb$Property == "Wallal Downs")] <- area[which(plb$Status == "Current" & plb$Property == "Wallal Downs")]
	r.samp[which(plb$Name == "La Grange")] <- area[which(plb$Name == "La Grange")]
	return(r.samp)
}

#Function called by spread function to set unique lambda values for irrigation points. The number of dispersers coming out of 
#an irrigated area is modified to be equal to the average number of disperses per unit area in other areas. 
unique.lambda <- function(occp,delta){
	irr1 <- occp[,1] == which(plb$Status == "Current" & plb$Property == "Shamrock Gardens")
	occp[,10][irr1] <- delta*0.03674617*3.1233 #Mean number of disperses per WB x density of waterbodies x area of irrigation
	irr2 <- occp[,1] == which(plb$Status == "Current" & plb$Property == "Shalemar")
	occp[,10][irr2] <- delta*0.03674617*13.688
	irr3 <- occp[,1] == which(plb$Status == "Current" & plb$Property == "Pardoo")
	occp[,10][irr3] <- delta*0.03674617*0.809
	irr4 <- occp[,1] == which(plb$Status == "Current" & plb$Property == "Wallal Downs")
	occp[,10][irr4] <- delta*0.03674617*8.9721
	town <- occp[,1] == which(plb$Name == "La Grange")
	occp[,10][town] <- delta*0.03674617*0.68
	return(occp)
}

#takes x, y, indices of a matrix and places them into a vector
flatmat<-function(x, y, size){ #size is size of one side of matrix
	(y-1)*size+x
}

neighbours.init<-function(space.size, cell.size){
	lth<-space.size/cell.size
	neigh<-array(0, dim=c(9,2,lth,lth))	
	for (x in 1:lth){
		for (y in 1:lth){
			out<-matrix(nrow=9, ncol=2)
			ifelse((x-1)==0, X<-NA, X<-x-1) #correct for indices going to zero
			ifelse((y-1)==0, Y<-NA, Y<-y-1)
			ifelse((y+1)>lth, Yp<-NA, Yp<-y+1) # correct for indices going to > space
			ifelse((x+1)>lth, Xp<-NA, Xp<-x+1)
			# Assign neighbours
			out[1,]<-c(x, y) #record target cell
			out[2,]<-c(Xp, y)
			out[3,]<-c(X, y)
			out[4,]<-c(x, Y)
			out[5,]<-c(Xp, Y)	
			out[6,]<-c(X, Y)
			out[7,]<-c(x, Yp)
			out[8,]<-c(Xp, Yp)
			out[9,]<-c(X, Yp)
			neigh[,,x,y]<-out
		}	
	}
	neigh	#return the array
}


pdist.fast<-function(X, Y, maximum, space.size){
	X<-X-min(X) # start coordinates at zero
	Y<-Y-min(Y)
	lth<-space.size/maximum
	if (lth%%1!=0) {
		print("Error: maximum must divide into space.size perfectly")
		return(NULL)	
	}
	out<-vector("list", length=length(X)) #list to take neightbour.ID and distance
	points<-vector("list", length=lth^2)#matrix of lists of point IDs
	gX<-X%/%maximum+1 #collapse to grid refs
	gY<-Y%/%maximum+1
	neigh<-neighbours.init(space.size, maximum)
	for (i in 1:length(X)){ #throw point IDs into grid cell list
		temp<-flatmat(gX[i], gY[i], lth)
		points[[temp]]<-c(points[[temp]], i)
	}
	for (i in 1:length(X)){
		nb<-neigh[,,gX[i], gY[i]] #find neighbouring cells
		nb<-subset(nb, is.na(apply(nb,1,sum))==F)
		temp<-flatmat(nb[,1], nb[,2], lth) #find relevant points
		snk.ID<-unlist(points[temp])
		dists<-sqrt((X[i]-X[snk.ID])^2+(Y[i]-Y[snk.ID])^2)
		temp<-cbind(snk.ID, dists)
		temp<-subset(temp, temp[,"dists"]<=maximum)
		out[[i]]<-temp
	}	
	out
}

#The kernel
dcncross.trunc<-function(x, u, v) {  #stuart's cauchy-normal distribution in 2D
	disp.old <-(x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v))) #1D dispersal kernel
	disp <- disp.new <- rep(0,length(x))
	maxd <- 65000
	for (i in 1:length(u)) {
		x1 <- x[i]
		u1 <- u[i]
		v1 <- v[i]
		if (x1 < maxd) {
		integrand <- function(x1) {(x1*u1^v1*v1*sqrt(v1^v1*(u1^2*v1+x1^2)^(-2-v1)))}
		area <- integrate(integrand, lower = 0, upper = maxd)
		disp.new[i] <- disp.old[i]/area$value
		disp[i] <- disp.new[i]/(2*pi*x1)} #2d dispersal kernel
	}
return(disp)
} 

dcncross.trunc2<-function(x, u, v) {  #stuart's cauchy-normal distribution in 2D
	
	disp <- rep(0,length(x))
	for (i in 1:length(disp)) {
		if (disp[i] < 100){
		disp[i] <- (x[i]*u[i]^v[i]*v[i]*sqrt(v[i]^v[i]*(u[i]^2*v[i]+x[i]^2)^(-2-v[i])))/(2*pi*x[i])}
	}
return(disp)
} 	

#The kernel
dcncross<-function(x, u, v) {  #stuart's cauchy-normal distribution in 2D
	(x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v)))/(2*pi*x)
} 

#The kernel
dcncross.T<-function(x, u, v, tr) {  #stuart's cauchy-normal distribution in 2D
	(x*u^v*v*sqrt(v^v*(u^2*v+x^2)^(-2-v)))/tr/(2*pi*x)
} 

# outputs data.  With or without a plot as well.
output<-function(pop, gen, id, plot=T){
	dname<-paste(id, "pop", gen, ".txt", sep="")
	write.table(pop, dname, sep="\t", row.names=F)
	if (plot==T){
		pname<-paste(id, "plot", gen, ".png", sep="")
		plotter(pop, pname)
	}	
}

# spreads the population over gen generations, and compares predictions to observed spread
spread<-function(pop, gens, pairs, delta, r, obs){ #pairs is a list from pdist.fast
  preds<-NULL	
for (i in 1:gens){
		#if (i%%5==0) output(pop, i, K)
		occp<-subset(pop, pop[,"Pres"]==1) #collect occupied sites
		occp<-cbind(occp, lambda=rpois(nrow(occp), delta))
		gma<-sum(occp[,"lambda"])
		if (length(occp[,1])==length(pop[,1])) {
			print(paste("No more vacant opportunities at generation", i))
			break()	
		}
		potl<-pairs[occp[,"ID"]] #collect relevant parts of pair list
		potl<-do.call("rbind", potl)
		src.ID<-rep(occp[,"ID"], times=occp[,"n.pairs"])
		potl<-cbind(src.ID, potl)
		lambda<-rep(occp[,"lambda"], times=occp[,"n.pairs"])
		U<-rep(occp[,"u"], times=occp[,"n.pairs"]) #expand source specific kernel parameters
		V<-rep(occp[,"v"], times=occp[,"n.pairs"])
		recruits<-lambda*dcncross(potl[,"dists"]+0.05, U, V) #calculate densities attributable to each pair
		recruits<-(pi*r^2*neigh.corr(pairs, r))/gma*tapply(recruits, potl[,"snk.ID"], sum) #sum densities from colonised waterbodies over all waterbodies and convert to proportion
		failures<-1-sum(recruits)
		if(failures<0) failures<-0 # catches the approximately statement (primarily happens at large r)
		recruits<-c(recruits, failures) #add on the failures
		recruit.ID<-as.integer(names(recruits)[-length(recruits)])
		recruits<-rmultinom(1, gma, recruits)[-length(recruits)]
		recruits<-cbind(recruit.ID, recruits)
		recruits<-subset(recruits, recruits[,"recruits"]>2)
		pop[match(recruits[,"recruit.ID"], pop[,"ID"]), "Pres"]<-1
		pop[which(pop[,"Pres"]==1), "age"]<-1+pop[which(pop[,"Pres"]==1), "age"]
		preds<-rbind(preds,output_val(pop, i))
	#plotter(pop, file.name=paste("gen",i,".png", sep=""))
	#print(obs_pred_cf(preds,obs))
	}
		obs_pred_cf(preds,obs)
}

# calculates the product of all elements in a vector
vec.prod<-function(vec){
  last<-length(vec)
  cumprod(vec)[last]
}

#plots opportunities and colonised populations
plotter<-function(popmatrix, file.name="temp.png", gen){
	png(filename=file.name, width=7, height=7, units="cm", res=150, pointsize=6)
	plot(popmatrix[,2], popmatrix[,3], xlab="False easting (kms)", ylab="False northing (kms)", pch=19)
	occp<-subset(popmatrix, popmatrix[,"Pres"]==1)
	points(occp[,2], occp[,3], pch=19, col="red")
	legend('topleft', legend=paste("Time =", gen), bty="n", pch=NA, cex=1.5)
	dev.off()
}

output_val<-function(pop, gens){
		pop<-pop[,1:4]
		pop<-cbind(pop,rep(gens, nrow(pop)))
		colnames(pop)[5]<-"Generation"	
		pop
}

obs_pred_cf<-function(preds, obs){
	preds[,5]<-preds[,5]+2007
	preds<-as.data.frame(preds)
	temp<-merge(preds, obs, by.x=c(2, 3, 5), by.y=c(9, 10, 15))
	#browser()
	temp<-temp[,"Pres"]==temp[,"OCCUPIED"]
	temp<-na.exclude(temp)
	sum(temp)
}

# spreads the population over gens generations or until target sites are reached. Includes a failure rate of managed waterpoints
spread.pilb.failure<-function(pop, gens, pairs, target, delta, r, trunc.dist, failure.rate){ #pairs is a list from pdist.fast  
	for (i in 1:gens){ #loop through time
		pop[,"Man"] <- managed	
		#Wet season starts and toads can spread
		u<-(rain[,i]-1)/364
		u<-3*(u-u^2) + u^3
		u<-rain[,i]+3*rain[,i]*(1-u)
		u<-floor(u) 
			if (scenario==1){u<=u}
			if (scenario==2){u[]<-180} #Set to 180 for pessermistic scenario
			if (scenario==3){u[]<-50}  #Set to 50 for optimistic scenario
			if (trunc.dist == 30000) {limit <- 4}
			if (trunc.dist == 55000) {limit <- 5}
			if (trunc.dist == 78000) {limit <- 6}
		trc <- fits[u,limit] #Extracts results of integration from fits dataset
			if (trunc.kernel == FALSE) {trc[] <- 1}
		u<-fits[u,1:2] #Number of rainy days per year given time and modelling scenario 
		u<-cbind(u,trc)
		pop[,6:8] <- u #Put values for u and v given time into population matrix
		occp<-subset(pop, pop[,"Pres"]==1) #collect occupied sites
		occp<-cbind(occp, lambda=rpois(nrow(occp), delta)) #Adds lambda values to each waterbody
		occp <- unique.lambda(occp,delta) #Modify lambda for irrigated areas plus La Grange
		gma<-sum(occp[,"lambda"]) #Total number of propagules in current time step
		potl<-pairs[occp[,"ID"]] #collects relevant waterbodies from pairs matrix so distance from occupied points to all other points is available.
		potl<-do.call("rbind", potl) #Combines elements of the matrix into a single vector
		src.ID<-rep(occp[,"ID"], times=occp[,"n.pairs"])
		potl<-cbind(src.ID, potl) #Creates a matrix with distance from each occupied source to all available waterbodies
		lambda<-rep(occp[,"lambda"], times=occp[,"n.pairs"])
		U<-rep(occp[,"u"], times=occp[,"n.pairs"]) #expand source specific kernel parameters
		V<-rep(occp[,"v"], times=occp[,"n.pairs"])
		T<-rep(occp[,"trc"], times=occp[,"n.pairs"])
		#recruits<-lambda*dcncross(potl[,"dists"]+0.05, U, V) #Number of possible dispersers multiplied by probability of dispersing to a particular distance
		recruits<-lambda*dcncross.T(potl[,"dists"]+0.05, U, V, T) #Number of
		trun<-(potl[,"dists"]+0.05)>trunc.dist #Identify points greater than the specified truncation distance
			if (trunc.kernel == TRUE) {recruits[trun==TRUE]<-0} #Make points greater than truncation distance equal to zero
		recruits<-(pi*r^2*neigh.corr.uniqueR(pairs, r))/gma*tapply(recruits, potl[,"snk.ID"], sum) #sum densities from colonised waterbodies over all waterbodies and convert to proportion
		failures<-1-sum(recruits)
		if(failures<0) failures<-0 # catches the approximately statement (primarily happens at large r)
		recruits<-c(recruits, failures) #add on the failures
		recruit.ID<-as.integer(names(recruits)[-length(recruits)])
		recruits<-rmultinom(1, gma, recruits)[-length(recruits)]
		recruits<-cbind(recruit.ID, recruits) #List all waterbodies and the number of recruits that disperse from the currently occupied points
		recruits<-subset(recruits, recruits[,"recruits"]>2) #Exclude waterbodies with less than 2 new recruits
		#Wet season ends and toads must find permanent water
		pop[match(recruits[,"recruit.ID"], pop[,"ID"]), "Pres"]<-1
		pop[which(pop[,"Pres"]==1), "age"]<-1+pop[which(pop[,"Pres"]==1), "age"]
		#Include failure rate of managed waterbodies
		for (z in 1:length(pop[,"Man"])) {
			ifelse(runif(1) > failure.rate & pop[,"Man"][z] == 1, pop[,"Man"][z] <- 0, pop[,"Man"][z] <- pop[,"Man"][z])} #Determine which managed waterbodies fail
		#Make all waterbodies that have been managed and haven't failed absent of toads
		pop[,"Pres"][which(pop[,"Man"]==1)] <- 0
		pop[,"age"][which(pop[,"Man"]==1)] <- 0
		if (plot.spread == TRUE) { #plot simulations
			par(mfrow=c(1,1))
			plot(pop[,2:3], pch=20, col="black", cex=1.5)
			points(subset(pop,pop[,4] == 1)[,2:3], col="red",pch=20, cex=1.5)
			points(subset(pop,pop[,10] == 1)[,2:3], col="grey",pch=20, cex=1.5)
		}
    if (sum(pop[target,"Pres"])>0) {break}
	}	
	list(gen=i, popmatrix=pop)	
}

# spreads the population over gens generations or until target sites are reached
# returns number of generations
# target is a vector of rows of pop that contain targets
# differs from previous versions of spread in that popmatrix has ndays of rain rather than U.
#	u is calculated internally
spread.pilb.varndays<-function(pop, gens, pairs, target, delta, r, plot=TRUE, fits, p.extMane, extMane.val){ #pairs is a list from pdist.fast  
for (i in 1:gens){
		#if (i%%5==0) output(pop, i, K)
		occp<-subset(pop, pop[,"Pres"]==1) #collect occupied sites
		occp<-cbind(occp, lambda=rpois(nrow(occp), delta))
		gma<-sum(occp[,"lambda"])
		potl<-pairs[occp[,"ID"]] #collect relevant parts of pair list
		potl<-do.call("rbind", potl)
		src.ID<-rep(occp[,"ID"], times=occp[,"n.pairs"])
		potl<-cbind(src.ID, potl)
		lambda<-rep(occp[,"lambda"], times=occp[,"n.pairs"])
			relrain<-1+rbinom(1, 1, p.extMane)*extMane.val
			ndays<-occp[,"ndays"]*relrain
			U<-(ndays-1)/364
			U<-3*(U-U^2) + U^3
			U<-ndays+3*ndays*(1-U)
			U<-floor(U)
			U<-fits[U,1:2]
			V<-rep(U[,2], times=occp[,"n.pairs"])
			U<-rep(U[,1], times=occp[,"n.pairs"]) #expand source specific kernel parameters			
		recruits<-lambda*dcncross(potl[,"dists"]+0.05, U, V) #calculate densities attributable to each pair
		#browser()
		recruits<-(pi*r^2*neigh.corr(pairs, r))/gma*tapply(recruits, potl[,"snk.ID"], sum) #sum densities from colonised waterbodies over all waterbodies and convert to proportion
		failures<-1-sum(recruits)
		if(failures<0) failures<-0 # catches the approximately statement (primarily happens at large r)
		recruits<-c(recruits, failures) #add on the failures
		recruit.ID<-as.integer(names(recruits)[-length(recruits)])
		recruits<-rmultinom(1, gma, recruits)[-length(recruits)]
		recruits<-cbind(recruit.ID, recruits)
		recruits<-subset(recruits, recruits[,"recruits"]>2)
		pop[match(recruits[,"recruit.ID"], pop[,"ID"]), "Pres"]<-1
		pop[which(pop[,"Pres"]==1), "age"]<-1+pop[which(pop[,"Pres"]==1), "age"]
		if (plot==TRUE) plotter(pop, file.name=paste(i,".png", sep=""), gen=i)
    if (sum(pop[target,"Pres"])>0) break
	}
	list(gen=i, popmatrix=pop)	
}

# given a point in the pdist list, takes out that point and its n nearest neighbours from the
#   pairwise distance list.  Returns modified list.
knock.out.nn<-function(pdist.list, point, n.n.neighb, natural){
  if (length(point)>1) warning("Multiple point Manoval not allowed")
  temp<-pdist.list[[point]] # get ID of points to knock out
  temp<-temp[order(temp[,'dists'])[1:(n.n.neighb+1)],"snk.ID"]
  temp<-temp[!temp%in%natural] #can't knock out natural points
  lapply(pdist.list, function(x) {subset(x, !x[,1]%in%temp)}) # and Manove them from everywhere in the list
}

# given a point not in the spread table, takes out that point's n nearest (artificial) neighbours from the
#   spread table.  Returns modified spread.table and pdist.list.
knock.out.nn.xy<-function(X, Y, spread.table, n, natural){
  if (length(X)>1 | length(Y)>1) {warning("Multiple point Manoval not allowed"); return(NULL)}
  pdists<-sqrt((spread.table[,"X"]-X)^2+(spread.table[,"Y"]-Y)^2) #Calculates distance of all waterbodies from the location of interest
  top<-order(pdists)[!order(pdists)%in%natural] #orders artificial waterbodies by distance from location of interest  
  spread.table<-spread.table[-top[1:n],] #Manoves the n closest waterbodies from the spread table
  pairs.mod<-pdist.fast(spread.table[,"X"], spread.table[,"Y"], maximum=500000, space.size=500000) #Calculates distance between now spread.table (with n points closest to location A Manoved)
  spread.table[, "n.pairs"]<-do.call("c",lapply(pairs.mod,nrow)) #Modify the number of pairs that each point has in the spread.table
  spread.table[, "ID"]<-1:nrow(spread.table) #Modify the ID of the Manaining waterbodies
  list(spread.table=spread.table, pairs.mod=pairs.mod)
}

#This function Removes waterbodies at a site until a fixed budget is spent
knock.out.nn.budget<-function(X, Y, spread.table, n, natural, irrigation, dwelling, Pres, cost.WP, rmnat, rmirr, rmdwel) {
  if (length(X)>1 | length(Y)>1) {warning("Multiple point Removal not allowed"); return(NULL)}
  pdists<-sqrt((spread.table[,"X"]-X)^2+(spread.table[,"Y"]-Y)^2) #Calculates distance of all waterbodies from the location of interest
  top <- order(pdists)
  if (rmnat==FALSE) {top<-top[!top%in%natural]} 
  if (rmirr==FALSE) {top<-top[!top%in%irrigation]} 
  if (rmdwel==FALSE) {top<-top[!top%in%dwelling]}
  top<-top[!top%in%which(Pres>0)]
  cost <- 0
  Rpoints <- c()
  i <- 1
  while (cost <= n) {
		if (cost.WP[top[1]] > n) {break} #Stop if the first waterpoint costs more than the available budget
		cost <- cost + cost.WP[top[i]]
		Rpoints <- c(Rpoints,top[i])
		i <- i + 1
		if (i == (length(top)+1)) {break} #Stop if budget exceeds the cost of managing all waterpoints
	}
  Rpoints<-Rpoints[1:(length(Rpoints)-1)] #Remove last waterpoint because this always exceeds the available budget
  spread.table[,"Man"][Rpoints] <- 1
  pairs.mod<-pdist.fast(spread.table[,"X"], spread.table[,"Y"], maximum=500000, space.size=500000)
  list(spread.table=spread.table, pairs.mod=pairs.mod)
}

# given a point, identifies n points making the shortest overall pathway between points
pathway<-function(pdist.list, point, n.points=2, natural){
  if (length(point)>1) {warning("Multiple point specification not allowed"); return(NULL)}
  IDs<-point
  temp<-pdist.list[[point]] # get first matrix
  for (ii in 1:(n.points-1)){
    temp<-temp[!temp[,"snk.ID"]%in%IDs,] #Manove rows already identified in IDs
    next.id<-temp[which(temp[,"dists"]==min(temp[,"dists"]))[1], "snk.ID"]
    temp<-rbind(temp, pdist.list[[next.id]])
    IDs<-c(IDs, next.id)
  }
  if (sum(IDs%in%natural)>0) return(NULL)
  IDs
}

# given a point, takes out that point and its n sized shortest path from the
#   pairwise distance list.  Returns modified list.
knock.out.path<-function(pdist.list, point, n.points=2, natural){
  temp<-pathway(pdist.list, point, n.points=2, natural)
  if (is.null(temp)) return(NULL)
  lapply(pdist.list, function(x) {subset(x, !x[,1]%in%temp)}) # and Manove them from everywhere in the list
}

#calculates density of X2, Y2 along line given by X1, Y1   
density.line<-function(X1, Y1, X2, Y2, bw=bw.nrd(sqrt((X1[1]-X2)^2+(Y1[1]-Y2)^2)), adjust=1){
  pdists<-vector(mode="list", length=length(X1))
  pdens<-vector(mode="numeric", length=length(X1))
  for (ii in 1:length(X1)){ 
    pdists[[ii]]<-sqrt((X2-X1[ii])^2+(Y2-Y1[ii])^2)
    pdens[ii]<-sum(dnorm(x=pdists[[ii]], mean=0, sd=bw*adjust))#/length(X2)
  }
  list(x=X1, y=pdens)
}

# Approximate correction for pi*r2 calculation by the proportional overlap between waterbodies that are less than 2r distant from one another
neigh.corr<-function(pairs, r){
	#collect pairs less than 2r distant
	corrn<-function(x){
		temp<-matrix(x[x[,"dists"]<=2*r & x[,"dists"]>0,], ncol=2) #Identify if there are any points within a distance equal to 2 x detection radius from one another. First column stores ID, second column stores distance
		if (length(temp)==0) return(1) #If there are none within the radius, return a 1
		theta<-2*acos(temp[,2]/(2*r))
		out<-1-(theta-sin(theta))/(2*pi)
		out<-prod(out)
		out
	}
	unlist(lapply(pairs, corrn))
}


# Approximate correction for pi*r2 calculation by the proportional overlap between waterbodies that are less than 2r distant from one another
neigh.corr.uniqueR<-function(pairs, r){
	#collect pairs less than 2r distant
	corrn<-function(x){
		x <- cbind(x,r)
		temp <- matrix(subset(x,x[,"dists"]<=2*x[,"r"] & x[,"dists"]>0),ncol=3) 
		if (length(temp)==0) return(1) #If there are none within the radius, return a 1
		theta<-2*acos(temp[,2]/(2*temp[,3]))
		out<-1-(theta-sin(theta))/(2*pi)
		out<-prod(out)
		out
	}
	unlist(lapply(pairs, corrn))
}

corridor.smoother <- function(X, Y, plb, interval, plot.barriers){
wp2<-cbind(X,Y,plb$art_nat)
ld.points<-wp2[order(wp2[,1]),] 
lo <- loess(ld.points[,"Y"]~ld.points[,"X"],degree=2,span=0.3)
#plot(ld.points[,"X"],ld.points[,"Y"], col="blue",pch=20)
#lines(ld.points[,"X"],predict(lo), col='red', lwd=2)
smX <- ld.points[,"X"]
smY <- predict(lo)
sm.dist <- rep(0,length(smX))
for (i in 2:length(smX)) {
	sm.dist[i] <- sqrt((smX[1]-smX[i])^2+(smY[1]-smY[i])^2)}
XYinterval <- (max(sm.dist) - min(sm.dist))/interval
ld.chunks <- seq(sm.dist[1],tail(sm.dist,1),XYinterval)
position <- ld.chunks
for (i in 1:length(position)) { 
	position[i] <- which(abs(sm.dist-ld.chunks[i])==min(abs(sm.dist-ld.chunks[i])))}
corridor <- cbind(smX[position],smY[position])
colnames(corridor) <- c("X","Y")
if (plot.barriers == TRUE) {
	par(mfrow=c(1,1), mar=c(1.5,4.5,1.5,4.5)+.1, oma=c(5,6,5,3), mai=c(0.1,0.1,0,0.1))
	plot(ld.points[,"Y"] ~ ld.points[,"X"], col="black", pch=20, ylab="", xlab="", xaxt="n", yaxt="n")
	points(corridor[,"Y"] ~ corridor[,"X"], col="red", pch=20, cex=2)
}
return(corridor)
}


###########New cost model with two types of irrigation
cost.waterpoint2 <- function(plb,gens,dst,plot.cost) { #
cost.WP <- matrix(NA,nrow(plb),1)
purpose <- plb$Purpose #Irrigation, artificial, Dwelling or natural waterbody
type <- plb$Type #1 = dam, 2 = tank
Property <- plb$Property
Area <- plb$Area #Note area of irrigation is in km2, while perimater is in meters
perimeter <- plb$Perimeter/100 #Perimeter of fences in 100s of metres
disc<- c(1:gens)
for (cc in 1:nrow(plb)){
	if (purpose[cc] == "Artificial" & type[cc] == 1) #Tanks
	{tank.installation <- plb[cc,37]*1.50 + 32*100 + 8500 #Installation cost = Distance (nearest town * $1.50) + Labour (2 men, 2 days @ $100/hr) + Materials ($8500)
	tank.maintenance <- (plb[cc,37]*1.50 + 8*100 + 500)*(1-0.025)^disc #Maintenance = Distance (nearest town * $1.50) + Labour ( 1 man, 1 day @ $100/hr) + Materials ($500) 
	tank.repair <- tank.installation/50*(1-0.025)^disc #Repair = Installation cost / 50 subject to discount rate
	cost.WP[cc] <- tank.installation + sum(tank.maintenance) + sum(tank.repair)} #Installation cost (Distance + Labour + Materials)
	
	if (purpose[cc] == "Artificial" & type[cc] == 2) #Dams
	{tank.installation <- plb[cc,37]*1.50 + 32*100 + 8500 #Installation cost = Distance (nearest town * $1.50) + Labour (2 men, 2 days @ $100/hr) + Materials ($8500)
	tank.maintenance <- (plb[cc,37]*1.50 + 8*100 + 500)*(1-0.025)^disc #Maintenance = Distance (nearest town * $1.50) + Labour ( 1 man, 1 day @ $100/hr) + Materials ($500)
	tank.repair <- tank.installation/50*(1-0.025)^disc #Repair = Installation cost / 50 subject to discount rate
	cost.WP[cc] <- tank.installation + sum(tank.maintenance) + sum(tank.repair)} #Installation cost (Distance + Labour + Materials)
	
	if (purpose[cc] == "Irrigation" & !Property[cc] == "Shalemar") #Irrigation - Pay for hay production to shut down for a period
	{cost.WP[cc] <- sum(Area[cc]*100*30*275*(1-0.025)^disc)} #Convert km2 to ha and multiply by profit gained for hay per ha. Check where profit comes from. 
	
	if (purpose[cc] == "Irrigation" & Property[cc] == "Shalemar") #Irrigation - Fence year round irrigation
	{fence.installation <- plb[cc,37]*1.50 + 32*100*perimeter[cc] + 3000*perimeter[cc] #Installation = Distance (nearest town * $1.50) + Labour (Perimeter * time per unit @ $100/hr) + Materials (Cost per unit length * Perimeter)
	fence.maintenance <- (8*100*26+100)*(1-0.025)^disc #Maintenance = Labour (1 day per week @ $100 per/h for 6 months) + cost ($100)
	fence.repair <- fence.installation/10*(1-0.025)^disc #Repair = Installation cost / 10 subject to discount rate
	cost.WP[cc] <- fence.installation + sum(fence.maintenance) + sum(fence.repair)}
	if (purpose[cc] == "Dwelling") {cost.WP[cc] <- 0} #Dwelling
	if (purpose[cc] == "Natural") {cost.WP[cc] <- 0} #Don't manage natural waterpoints
	}
	
if (plot.cost == TRUE) {
	jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	plot(X,Y, col=jet.colors(12)[cut(cost.WP, 12)], pch=20, cex=1)}
	return(cost.WP)	
}