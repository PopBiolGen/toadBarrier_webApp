###################################################
# from custom simulation
fName<-"FIrr0IPA0Dwell1K78S1Nat0Irr1.RData"

load(paste0("out/", fName))

temp<-do.call(what="rbind", output)

times<-unlist(temp[,"gen"])
point<-unlist(temp[,"location"])
nns<-unlist(temp[,"Removed"])
occ<-unlist(temp[,"colonised"])
cost<-unlist(temp[,"cost"])

reps <- table(cost, point)[1,1]
#Proportion of simulations where toads spread to target before 100 years
prop<-tapply(times<50, list(cost, point), sum)/reps

#Number of waterpoints removed
nn <- tapply(nns, list(cost, point), mean)

#Cost of excluding nn waterpoints
cost.m<-tapply(cost, list(cost, point), mean)

#matplot(cost.m, prop)

#matplot(nn, prop)

probOutIPA0Irr1<-list(nWb=nn, spend=cost.m, pSucc=prop)


save(probOutIPA0Irr1, file = "webApp/pilbaraLine/modelData/IPA0Irr1.RData")

rm(list=ls())


###################################################
fName<-"FIrr0IPA1Dwell1K78S1Nat0Irr1.RData"

load(paste0("out/", fName))

temp<-do.call(what="rbind", output)

times<-unlist(temp[,"gen"])
point<-unlist(temp[,"location"])
nns<-unlist(temp[,"Removed"])
occ<-unlist(temp[,"colonised"])
cost<-unlist(temp[,"cost"])

reps <- table(cost, point)[1,1]
#Proportion of simulations where toads spread to target before 100 years
prop<-tapply(times<50, list(cost, point), sum)/reps

#Number of waterpoints removed
nn <- tapply(nns, list(cost, point), mean)

#Cost of excluding nn waterpoints
cost.m<-tapply(cost, list(cost, point), mean)

#matplot(cost.m, prop)

#matplot(nn, prop)

probOutIPA1Irr1<-list(nWb=nn, spend=cost.m, pSucc=prop)

save(probOutIPA1Irr1, file = "webApp/pilbaraLine/modelData/IPA1Irr1.RData")

rm(list=ls())
