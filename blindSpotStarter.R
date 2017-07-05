
blindSpotStarter<-function(elementSymbols,periodic=TRUE,round){
#create the expand.grid matrix based on atoms in the periodic table
#need to write this code!
#factor in the full masses and resolutions
#choose("C","H","O","S") and find masses
require(FRACTION)
elementChoose<-function(periodic=TRUE,elementSymbols,round){
if(periodic){
units<-c(0,0)
for(i in 1:length(elementSymbols)){
	units<-rbind(units,c(0,as.numeric(nuclides[nuclides[,"Symbol"]==elementSymbols[i],"Mass"])))
#	nuclides[,"Symbol"]=="H(1)"|
#	nuclides[,"Symbol"]=="O(16)"|
#	nuclides[,"Symbol"]=="S(32)","Mass"]))

	}
	units

}
 else {units<-cbind(0,elementSymbols)}

units
units<-round(units,round)
}
units<-elementChoose(periodic=TRUE,elementSymbols,round)

eqs<-expand.grid(lapply(apply(units,1,list),function(x) unlist(x)))

#eqs<-as.matrix(expand.grid(c(0,12),c(0,1),c(0,16),c(0,32)))
colnames(eqs)[-1]<-elementSymbols
eqs<-apply(eqs,1,function(x) x[grep("[1-9]{1,3}",x)])
eqs<-lapply(eqs,function(x) if(length(x)>1){x})
eqs<-eqs[!sapply(eqs,is.null)]
eqs<-lapply(eqs,function(x) c(x[1],-(x[-1]))/x[1])
#eqs is now the set of difference equations.
#Now find lowest multiple

denom<-lapply(eqs,function(x) as.numeric(unlist(regmatches(fra.m(x),gregexpr("(?<=/ ).+",fra.m(x),perl=T)))))
#denom<-as.numeric(unlist(regmatches(fra.m(eqs[[11]]),gregexpr("(?<=/ )[1-9]+",fra.m(eqs[[11]]),perl=T))))
neg<-list(list())
for(i in 1:length(eqs)){neg[[i]]<-eqs[[i]][-1]*denom[[i]][-1]}
pos<-list(list())
for(i in 1:length(eqs)){pos[[i]]<-(-sum(eqs[[i]][-1]*denom[[i]][-1]))}
#assign names to elements in the vectors
for(i in 1:length(eqs)){
	names(denom[[i]])<-names(eqs[[i]])
	}
for(i in 1:length(eqs)){
	names(pos[[i]])<-names(eqs[[i]][1])
	}
final<-list(list())
for(i in 1:length(eqs)){
	final[[i]]<-c(pos[[i]],-denom[[i]][-1])
	}
final
}
