pth<-"../Paper_results/Fig2/9_IPAT_T55_S1_B21_AT_IF_B30/Results/"
fname<-"IPAT_T55_S1_B21_AT_IF_B30.RData"
load(file=paste0(pth, "times.", fname))
load(file=paste0(pth, "cost.", fname))
load(file=paste0(pth, "point.", fname))
load(file=paste0(pth, "nns.", fname))

reps <- 1000
#Proportion of simulations where toads spread to target before 100 years
prop<-tapply(times<50, list(cost, point), sum)/reps
#Number of waterpoints removed
nn <- tapply(nns, list(cost, point), mean)
#Cost of excluding nn waterpoints
cost.m<-tapply(cost, list(cost, point), mean)
#Delete barrier locations at far end of corridor
prop <- prop[,-c(1, 19:20)]
nn <- nn[,-c(1, 19:20)]

#Plot probability of reaching Pilbara vs Budget
r.obj <- prop
#palette(rainbow(nrow(corridor)))
plot(cost.m[,1]/1000000, r.obj[,1], ylim=c(0,1.1), xlab = "Budget (Million dollars)", ylab="Pr(Colonization)", col="gray", type="l", xlim=c(0,15))
for (i in 1:ncol(prop)) {
  #points(cost.m[,i]/1000000, r.obj[,i], col="grey", pch=20)
  lines(cost.m[,i]/1000000, r.obj[,i], col="grey", lwd=1)
}
text(0, 1.1, "c)",cex=1)
best<-prop[, order(prop[rownames(prop)=="9e+06",])] # Order by best performance at $9m
#best <- cbind(prop[,4], prop[,5], prop[,12], prop[,17])
points(cost.m[,1]/1000000, best[,1], col="red", pch=20)
lines(cost.m[,1]/1000000, best[,1], col="red", lwd=2)
points(cost.m[,2]/1000000, best[,2], col="green", pch=20)
lines(cost.m[,2]/1000000, best[,2], col="green", lwd=2)
points(cost.m[,3]/1000000, best[,3], col="orange", pch=20)
lines(cost.m[,3]/1000000, best[,3], col="orange", lwd=2)
points(cost.m[,4]/1000000, best[,4], col="cyan", pch=20)
lines(cost.m[,4]/1000000, best[,4], col="cyan", lwd=2)


pth<-"../Paper_results/Fig2/11_IPAT_T55_S1_B21_AT_IF_B30_noDwellings/Results/"
fname<-"IPAT_T55_S1_B21_AT_IF_B30_nodwell.RData"
load(file=paste0(pth, "times.", fname))
load(file=paste0(pth, "cost.", fname))
load(file=paste0(pth, "point.", fname))

reps <- 1000
#Proportion of simulations where toads spread to target before 100 years
prop<-tapply(times<50, list(cost, point), sum)/reps
matplot(as.numeric(rownames(prop)), y = prop)
 