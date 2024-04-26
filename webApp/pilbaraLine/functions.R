## Data manipulation functions from Darrens Paper

# identifies the subset to keep for cost calculations, based on assumptions
modify.plb <- function(plb, future.irrigation, current.irrigation = TRUE, IPA, Dwelling){
  ss<-rep(TRUE, nrow(plb)) # full subset
  #Future irrigation
  if (future.irrigation == FALSE) {
    ss <- ss & !(plb$Purpose == "Irrigation" & plb$Type == 2)
  }
  #Current irrigation
  if (current.irrigation == FALSE) {
    ss <- ss & !(plb$Purpose == "Irrigation" & plb$Type == 1)
  }	
  #IPA areas
  if (IPA == TRUE) {
    ss <- ss & !(plb$Property == "Frazier Downs" & plb$IPA==0)
  }
  #Dwellings
  if (Dwelling == FALSE) {
    ss <- ss & !(plb$Purpose == "Dwelling")
  }
  return(ss)
}

#Sets unique detection radii for irrigation and La Grange
detection.radius <- function(plb) {
  r.samp<-rep(10^2,nrow(plb)) #Detection radius -> becomes r
  area<-sqrt((plb$Area*1000000)/pi) #Convert area from square metres to square kilometers then calculate area equivalent
  r.samp[which(plb$Status == "Current" & plb$Property == "Shamrock Gardens")] <- area[which(plb$Status == "Current" & plb$Property == "Shamrock Gardens")]
  r.samp[which(plb$Status == "Current" & plb$Property == "Shalemar")] <- area[which(plb$Status == "Current" & plb$Property == "Shalemar")]
  r.samp[which(plb$Status == "Current" & plb$Property == "Pardoo")] <- area[which(plb$Status == "Current" & plb$Property == "Pardoo")]
  r.samp[which(plb$Status == "Current" & plb$Property == "Wallal Downs")] <- area[which(plb$Status == "Current" & plb$Property == "Wallal Downs")]
  r.samp[which(plb$Name == "La Grange")] <- area[which(plb$Name == "La Grange")]
  return(r.samp)
}

# Calculates locations of Test points in the paper
test.points<-function(plb.xy, barriers=21){
  ld.points<-plb.xy[order(plb.xy[,"Y"]),] 
  lo <- loess(ld.points[,"Y"]~ld.points[,"X"],degree=2,span=0.3)
  smX <- ld.points[,"X"]
  smY <- predict(lo)
  sm.dist <- rep(0,length(smX))
  for (i in 2:length(smX)) {
    sm.dist[i] <- sqrt((smX[1]-smX[i])^2+(smY[1]-smY[i])^2)}
  XYinterval <- (max(sm.dist) - min(sm.dist))/barriers
  ld.chunks <- seq(sm.dist[1],tail(sm.dist,1),XYinterval)
  position <- ld.chunks
  for (i in 1:length(position)) { 
    position[i] <- which(abs(sm.dist-ld.chunks[i])==min(abs(sm.dist-ld.chunks[i])))}
  corridor <- cbind(smX[position],smY[position])
  colnames(corridor) <- c("X","Y")
  #Delete points at far ends of corridor
  corridor <- corridor[-c(20:22),]
  corridor <- corridor[-c(1:2),]
  rownames(corridor)<-1:nrow(corridor)
  corridor
}

