nagenmzhrnew<-function(Cmin,Cmax,DBE,Nn,On,Sn,massmin,massmax,charge){
C<-Cmin:Cmax
Z<-2*(0:-DBE)
N<-0:Nn
O<-0:On
S<-0:Sn
new<-as.matrix(expand.grid(Z,C,N,O,S))
colnames(new)<-c("Z-number","C","N","O","S")
H<-2*new[,"C"]+2+new[,"N"]+new[,"Z-number"]
new<-cbind(new[,1:2],H,new[,3:5])
colnames(new)<-c("Z-number","C","H","N","O","S")
new<-new[new[,"H"]>0,]
neutmass<-12.000000*new[,"C"]+1.007825*new[,"H"]+15.994915*new[,"O"]+14.003074*new[,"N"]+31.972072*new[,"S"]
new<-cbind(neutmass,new)
m2cH<-neutmass+charge*(1.007825-0.0005)
m2cHd<-(neutmass+2*charge*(1.007825-0.0005))/2
new<-cbind(new,m2cH,m2cHd)
colnames(new)<-c("Neutral mass","Z-number","C","H","N","O","S","m/z","m/z_2")
new<-new[new[,"m/z"]<massmax & new[,"m/z"]>massmin,]
new
}
