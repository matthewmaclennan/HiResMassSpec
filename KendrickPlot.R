#Kendrick plotting

kendrickPlot<-function(data,cex,xlim,ylim){
km<-as.numeric(data[,1])*14/14.01565
#Kendrick Mass defect
kmd<-round(floor(km)-km,4)
plot(kmd~km,cex=cex,xlim=xlim,ylim=ylim)
}


kendrickPoints<-function(data,cex,xlim,ylim){
km<-as.numeric(data[,1])*14/14.01565
#Kendrick Mass defect
kmd<-round(floor(km)-km,4)
points(kmd~km,cex=cex,xlim=xlim,ylim=ylim)
}


kendrickPlotNm<-function(data,cex,xlim,ylim){
km<-as.numeric(data[,1])*14/14.01565
#Kendrick Mass defect
kmd<-round(floor(as.numeric(data[,1]))-km,4)
plot(kmd~km,cex=cex,xlim=xlim,ylim=ylim)
}


kendrickPointsNm<-function(data,cex,xlim,ylim){
km<-as.numeric(data[,1])*14/14.01565
#Kendrick Mass defect
kmd<-round(floor(as.numeric(data[,1]))-km,4)
points(kmd~km,cex=cex,xlim=xlim,ylim=ylim)
}




par(col="black",pch=16)
kendrickPlot(s14posxf.20[grep("O5",s14posxf.20[,3]),],cex=.5,xlim=c(100,600),ylim=NULL)
par(col="red",pch=16)
kendrickPoints(s14posxf.20[grep("O4",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="blue",pch=16)
kendrickPoints(s14posxf.20[grep("O3",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="dark green",cex=1,pch=16)
kendrickPoints(s14posxf.20[grep("O2",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)

kendrickPoints(s14posxf.20[grep("O1",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)





par(col="black",pch=16)
kendrickPlot(mgt0.20f[grep("O5",mgt0.20f[,3]),],cex=.5,xlim=c(100,600),ylim=c(-1,-.8))
par(col="red",pch=16)
kendrickPoints(mgt0.20f[grep("O4",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="blue",pch=16)
kendrickPoints(mgt0.20f[grep("O3",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="dark green",cex=1,pch=16)
kendrickPoints(mgt0.20f[grep("O2",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)

kendrickPoints(mgt0.20f[grep("O1",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)

###############################################################################################

par(col="black",pch=16)
kendrickPlotNm(s14posxf.20[grep("O5",s14posxf.20[,3]),],cex=.5,xlim=c(100,600),ylim=c(0,0.3))
par(col="red",pch=16)
kendrickPointsNm(s14posxf.20[grep("O4",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="blue",pch=16)
kendrickPointsNm(s14posxf.20[grep("O3",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="dark green",cex=1,pch=16)
kendrickPointsNm(s14posxf.20[grep("O2",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="dark orange",pch=16)
kendrickPointsNm(s14posxf.20[grep("O1",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="purple",pch=16)
kendrickPointsNm(s14posxf.20[grep("O6",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="cyan",pch=16)
kendrickPointsNm(s14posxf.20[grep("O7",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="green",pch=16)
kendrickPointsNm(s14posxf.20[grep("O8",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)
par(col="navy blue",pch=16)
kendrickPointsNm(s14posxf.20[grep("O9",s14posxf.20[,3]),],cex=0.5,xlim=c(100,600),ylim=NULL)



par(col="black",pch=16)
kendrickPlotNm(mgt0.20f[grep("O5",mgt0.20f[,3]),],cex=.5,xlim=c(100,400),ylim=c(0,0.3))
par(col="red",pch=16)
kendrickPointsNm(mgt0.20f[grep("O4",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="blue",pch=16)
kendrickPointsNm(mgt0.20f[grep("O3",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="dark green",cex=1,pch=16)
kendrickPointsNm(mgt0.20f[grep("O2",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="dark orange",pch=16)
kendrickPointsNm(mgt0.20f[grep("O6",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="purple",pch=16)
kendrickPointsNm(mgt0.20f[grep("O7",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="cyan",pch=16)
kendrickPointsNm(mgt0.20f[grep("O8",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="green",pch=16)
kendrickPointsNm(mgt0.20f[grep("O9",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="navy blue",pch=16)
kendrickPointsNm(mgt0.20f[grep("O10",mgt0.20f[,3]),],cex=0.5,xlim=NULL,ylim=NULL)

###############################################################################################



par(col="black",pch=16)
kendrickPlotNm(s16s14f2[grep("N0O2S0",s16s14f2[,3]),],cex=.5,xlim=c(100,400),ylim=c(0,0.3))
par(col="red",pch=16)
kendrickPointsNm(s16s14f2[grep("N0O1S0",s16s14f2[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="blue",pch=16)
kendrickPointsNm(s16s14f2[grep("N0O3S1",s16s14f2[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="dark green",cex=1,pch=16)
kendrickPointsNm(s16s14f2[grep("N0O3S0",s16s14f2[,3]),],cex=0.5,xlim=NULL,ylim=NULL)
par(col="dark orange",pch=16)
kendrickPointsNm(s16s14f2[grep("N0O3S2",s16s14f2[,3]),],cex=0.5,xlim=NULL,ylim=NULL)


O2labs<-as.numeric(unlist(regmatches(s16s14f2[grep("N0O2S0",s16s14f2[,3]),3],gregexpr(" -[0-9]+",s16s14f2[grep("N0O2S0",s16s14f2[,3]),3],perl=T))))
