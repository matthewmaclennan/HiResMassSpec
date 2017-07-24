#formulaBANK, m/z peak list, ppm threshold, charge of data
#soc
matchppm<-function(bank,mz,thresh,Charge){
	list1<-list()
	for(i in 1:length(mz)){
		list1[[i]]<-bank[abs((round(abs(as.numeric(bank[,"m/z"])),4)-mz[i])/mz[i]*1e6)<=thresh & bank[,"Charge"]==Charge,]
	}
list1
#lapply(list1,function(x) )
}
#eoc

matchppm2<-function(bank,mz,thresh,Charge){
	list1<-list()
	for(i in 1:length(mz)){
		list1[[i]]<-bank[abs((round(abs(as.numeric(bank[,"m/z"])),4)-mz[i])/mz[i]*1e6)<=thresh,]
	}
list1
#lapply(list1,function(x) )
}

matchppm_2<-function(bank,mz,thresh,Charge){
	list1<-list()
	for(i in 1:length(mz)){
		list1[[i]]<-bank[abs((round(abs(as.numeric(bank[,"m/z_2"])),4)-mz[i])/mz[i]*1e6)<=thresh,]
	}
list1
#lapply(list1,function(x) )
}


m<-matchppm2(contam,as.numeric(dmedapxf1[,1]),1,1)
matrix(unlist(lapply(m,function(x) x[length(x)>0])),ncol=8,byrow=T)

################################################################################################
scanproc<-function(s14p,mzmax){
#Select masses whose intensity is > 0, return s14p without zero-intensity data.
mgt0<-lapply(s14p,function(x) x[x[,2]>0,])
mgt0x<-lapply(mgt0,function(x) x[x[,1]<=mzmax,])
mgt0x
}
####################################################################################



####################################################################################
formulamatch2<-function(mgt0,nagenS0N0n1,thresh,Charge){
#match formulae to within threshold ppm
#case here uses S=0, N=0, O<=8, C5-35
if(Charge==1){
s14.20<-matchppm2(nagenS0N0n1,mgt0,thresh,Charge)}
if(Charge==2){
s14.20<-matchppm_2(nagenS0N0n1,mgt0,thresh,Charge)}
#list all matched formulae
#s14.20[unlist(lapply(s14.20,function(x) length(x)>0))]
#make into matrix


#Match formula to peaks
#Generate formula
#	f1<-paste0("C",s14.20mat[,3],"H",s14.20mat[,4],"N",s14.20mat[,6],"O",s14.20mat[,5],"S",s14.20mat[,7])
#	f2<-gsub("[A-Z][a-z]?0","",f1)
mgt0.20f<-mgt0
for (i in 0:(max(unlist(lapply(s14.20,function(x) length(x)/9)))-1)){
s14.20f<-lapply(lapply(s14.20,function(x) matrix(x,ncol=9,byrow=F)),function(y) if(nrow(y)>i) y[i+1,])
s14.20f<-unlist(lapply(s14.20f, function(x) paste0("C",x[3],"H",x[4],"N",x[5],"O",x[6],"S",x[7]," \\(Z= ",x[2],"\\)")))
s14.20f<-gsub("CNAHNANNAONASNA.+\\(Z= NA.\\)","",s14.20f)
s14.20f<-gsub("CHNOS.+\\(Z= .\\)","",s14.20f)
#s14.20f<-gsub("[A-Z][a-z]?0","",s14.20f)
s14.20f

#Place next to peak data
mgt0.20f<-cbind(mgt0.20f,s14.20f)
mgt0.20f<-mgt0.20f
}
mgt0.20f
##If you wish, you can go to the plotting protocol.	
}

###########################################################################
s14posxf<-formulamatch(s14posx,nagenS0N6p1[nagenS0N6p1[,"Charge"]==1,])

#############################################
#Kendrick mass 2D (H/C,O/C)
kendrickPlot<-function(data){
km<-as.numeric(data[,1])*14/14.01565
#Kendrick Mass defect
round(km-floor(km),4)
plot(round(km-floor(km),4)~km,cex=0.2)
}
#############################################
#C number dist

#Z number dist
#################################
#Van Krevelen ratios
#################################
#H/C
vanKrevelenPlot<-function(data){
H<-as.numeric(unlist(regmatches(data[,3],gregexpr("(?<=H)[0-9]+",data[,3],perl=T))))
C<-as.numeric(unlist(regmatches(data[,3],gregexpr("(?<=C)[0-9]+",data[,3],perl=T))))
O<-as.numeric(unlist(regmatches(data[,3],gregexpr("(?<=O)[0-9]+",data[,3],perl=T))))
N<-as.numeric(unlist(regmatches(data[,3],gregexpr("(?<=N)[0-9]+",data[,3],perl=T))))
Z<-as.numeric(unlist(regmatches(data[,3],gregexpr("\\-[0-9]+",data[,3],perl=T))))
listy<-list(C=C,H=H,Z=Z,N=N,O=O)
}

HoC<-H/C
OoC<-O/C
NoC<-N/C
plot(NoC~OoC,cex=0.5)
#}
#plot O series different colours
points(2/C[O==2],H[O==2]/C[O==2],cex=0.2,col="green",pch=7)
points(4/C[O==4],H[O==4]/C[O==4],cex=0.9,col="red",pch=7)
points(6/C[O==6],H[O==6]/C[O==6],cex=0.9,col="blue",pch=7)
points(8/C[O==8],H[O==8]/C[O==8],cex=0.9,col="purple",pch=7)
#Legend?
points(2/C[O==2],N[O==2]/C[O==2],cex=0.9,col="green",pch=7)
points(4/C[O==4],N[O==4]/C[O==4],cex=0.9,col="red",pch=7)
points(6/C[O==6],N[O==6]/C[O==6],cex=0.9,col="blue",pch=7)
points(8/C[O==8],N[O==8]/C[O==8],cex=0.9,col="purple",pch=7)

####################################
#eoc
####################################

#OMIT na C number and Z number of matched formulae
#na.omit(matrix(unlist(lapply(s14.20,function(x) x[c("C","Z-number")])),ncol=2,byrow=T))

#Plotting protocols!
#Plot the scan or data
plot(mgt0.20f[,1:2],type="h",xlim=as.numeric(c(min(mgt0.20f[,1]),max(mgt0.20f[,1]))))

#####################################################################################
##################### Plot over it some atom subclasses #############################
#####################################################################################
classplot<-function(data,subclass,col,alone=F){
if(alone==F){
points(data[grep(subclass,data[,3]),1:2],type="h",col=col)
}
	else {
		plot(data[grep(subclass,data[,3]),1:2],type="h",col=col,xlim=as.numeric(c(min(mgt0.20f[,1]),max(mgt0.20f[,1]))))
	}
}
#####################################################################################

#Sum the classes
#third column refers to mass scan with cbind() formula assignments.
classsum<-function(data,class){
listy<-list()

for(i in 1:length(class)){
listy[[i]]<-sum(as.numeric(data[grep(class[i],data[,3]),2]))
}
listy
}
#######################################################
speciesplot<-function(data,sc){
tot<-sum(as.numeric(data[,2]))
unlist(classsum(data,sc))/tot
barplot(unlist(classsum(data,sc))/tot,names.arg=sc)
}
#######################################################

#######################################################
#plot C number versus DBE



#######################################################

Match positive mode ions to negative mode ions

#O2 negative data
mgt0.20f[grep("O2",mgt0.20f[,3]),]
#masses
mgt0.20f[grep("O2",mgt0.20f[,3]),1]
#plot 'em
plot(mgt0.20f[grep("O2",mgt0.20f[,3]),],type="h",col="red",xlim=c(100,500))
points(s14posxf.20[grep("N3O2",s14posxf.20[,3]),],type="h")
points(s14posxf.20[grep("N2O",s14posxf.20[,3]),],type="h",col="blue")
points(as.numeric(s14posxf.20[grep("N2O1",s14posxf.20[,3]),1])-71,s14posxf.20[grep("N2O1",s14posxf.20[,3]),2],type="h",col="blue")

plot(s14posxf.20[grep("N2O3",s14posxf.20[,3]),1:2],type="h",col="orange")
points(as.numeric(s14posxf.20[grep("N2O3",s14posxf.20[,3]),1])-156,s14posxf.20[grep("N2O3",s14posxf.20[,3]),2],type="h",col="orange")

#matching gunk
#negative ion spectrum masses
neg<-as.numeric(mgt0.20f[,1])
#arbitrary positive ion mass
#subtract from all negative points
#use differences to compare with known derivatization differences

multiClassPlot<-function(){
#Need to ad in the i j counters for loops.

pos<-as.numeric(s14posxf.20[grep("N2O",s14posxf.20[,3]),][37,1])
###############################################################################

derivshift<-function(pos,neg,derivatives,thresh){
#matches positive MS data of derivatized compounds to negative MS data of the same underivatized 
#compounds using a derivative match.
if(length(pos)>0){
listy<-list()
for(j in 1:length(derivatives)){
listy[[j]]<-list()
for(i in 1:length(pos)){
listy[[j]][[i]]<-list()
listy[[j]][[i]]<-neg[abs((abs(pos[i]-as.numeric(neg[,1]))-as.numeric(derivatives[j]))/as.numeric(derivatives[j])*1e6)<=thresh,]
}
listy
}
listy
}}
###############################################################################

derivativePlot<-function(mgt0.20f,s14posxf.20,formularule,derivatives,thresh){
#mgt0.20f is the negative ion data
require(plyr)
#s14posxf.20 is the positive ion data
matchmat<-cbind(unique(sort(mgt0.20f[nchar(mgt0.20f[,3])>0,3])),matrix(0,ncol=length(formularule),nrow=length(unique(sort(mgt0.20f[nchar(mgt0.20f[,3])>0,3])))))
#print(matchmat)
#pos is m/z values on the formularule hits

for(i in 1:(length(formularule)-1)){
print(i)
pos<-as.numeric(s14posxf.20[grep(formularule[i],s14posxf.20[,3]),1])
#print(pos)
if(length(pos)>0){
der<-derivshift(pos,mgt0.20f[,1:3],derivatives,thresh)
#mgt0.20f[,1:3] should include the columns m/z, intensity and formula in that order.
der<-lapply(unlist(der,recursive=F),function(x)matrix(x,ncol=3,byrow=F))
#omit rows of NA 
der<-lapply(der,function(x) na.omit(x))
if(length(unlist(der))>2){

#print(der)
#mat<-unlist(lapply(der,function(x) x[grep("C",x)]))
#mat<-unlist(lapply(unlist(der,recursive=F),function(x) x[,3]))
#mat<-matrix(unlist(lapply(der,function(x) if(sum(as.numeric(x[,2]))) cbind(mean(as.numeric(x[,2])),unique(x[,3])))),ncol=2,byrow=T)
#mat<-mat[grep("C",mat[,2]),]
mat<-as.matrix(ldply(der,function(x) if(sum(as.numeric(x[,2]))) {cbind(mean(as.numeric(x[,1]),trim=0),mean(as.numeric(x[,2]),trim=0),unique(x[,3]))}))
mat<-mat[grep("C",mat[,3]),]
print(mat)
##MUST DEAL WITH 1-ROW RESULTS!
if(length(mat)>3){
#check for duplicate entries and remove them
if(sum(duplicated(mat[,3],fromLast=T)>0)){
insert<-c(mean(c(as.numeric(mat[duplicated(mat[,3],fromLast=T),1]),as.numeric(mat[duplicated(mat[,3]),1]))),
	mean(c(as.numeric(mat[duplicated(mat[,3],fromLast=T),2]),as.numeric(mat[duplicated(mat[,3]),2]))),
	mat[duplicated(mat[,3]),3])
print(insert)
mat<-mat[!duplicated(mat[,3]),]
mat[duplicated(mat[,3],fromLast=T),]<-insert
print(mat)
#print(matchmat)
}
matchmat[,i+1]<-match(matchmat[,1],mat[,3],nomatch=0)
matchmat[matchmat[,i+1]>0,i+1]<-mat[,2]
print(matchmat)
} else {mat<-matrix(mat,ncol=3,byrow=T)
		matchmat[,i+1]<-match(matchmat[,1],mat[,3],nomatch=0)
		matchmat[matchmat[,i+1]>0,i+1]<-mat[,2]
		print(matchmat)}

matchmat
}
matchmat
}
matchmat
}
matchmat
#Can use the command below to name the columns. Switch the names if inappropriate.
#colnames(matchmat)<-c("Neg Ion","EDC","NHS","DMEDA","EDC_EDC","EDC_NHS","EDC_DMEDA","NHS_NHS","NHS_DMEDA","DMEDA_DMEDA")
#The next iteration of this script is to sum / average the intensities of peaks in the specific class
#in order to do class abundance.
#An important thing is to plot this information appropriately for maximum information transmission.
}

pos<-as.numeric(s14posxf.20[grep(formularule[19],s14posxf.20[,3]),1])
der<-derivshift(pos,mgt0.20f[,1:3],as.numeric(derivatives[,3]),1)
der<-lapply(unlist(der,recursive=F),function(x)matrix(x,ncol=3,byrow=F))
mat<-as.matrix(ldply(der,function(x) if(sum(as.numeric(x[,2]))) {cbind(mean(as.numeric(x[,1]),trim=0),mean(as.numeric(x[,2]),trim=0),unique(x[,3]))}))
mat<-mat[grep("C",mat[,3]),]


#Order MS data of a particular atom class (by regex) according to abundance. The purpose 
#of this was to see the maximum peak height of a particular class to determine if that type 
#of derivatization occurred often. Eventually a suite of algorithms will be used to compare 
#abundances of negative underivatized with positive derivatized.
##############################################################################
maxclass<-function(s14posxf.20,class){
s14posxf.20[grep(class,s14posxf.20[,3]),][order(as.numeric(s14posxf.20[grep(class,s14posxf.20[,3]),2])),]
}
##############################################################################
zFrFormula<-function(formula,append=T){
#Compute Z number form a chemical formula
#append=T decides to place the Z number next to the formula

}
##############################################################################

#derivatization matches
derivs<-function(s22s4DmgPxf,formulaclass,derivatives,cutoff){
c<-c()
for(i in 1:length(derivatives[,3])){
c<-cbind(c,as.numeric(s22s4DmgPxf[grep(formulaclass,s22s4DmgPxf[,3]),][as.numeric(s22s4DmgPxf[grep(formulaclass,s22s4DmgPxf[,3]),2])>cutoff,1])-as.numeric(derivatives[i,3]))
}
colnames(c)<-derivatives[,1]
c
}

#fragmentation matches
frags<-


########################################################################################################################################
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


########################################################################################################################################

#derivatives[,2] so far can be used to add to neutral precursor ("standards")
#subtract fragments (Losses)
#Losses contains neutral losses unique to derivatization reagents
#contam is the Contamination file -> check peak masses to these
#create character string descriptions of each value obtained
#NAMES by expand.grid() and paste0()

###############################################################################################
standfrag<-function(standards,derivatives,Losses){
#Create lists of standards, derivatives, and common ESI-MS losses/additions
a<-unlist(lapply(apply(as.matrix(expand.grid(standards[,2],derivatives[,2],Losses[,2])),1,function(x) paste0(x,collapse="+")),function(x) eval(parse(text=x))))
b<-paste0("[",apply(as.matrix(expand.grid(standards[,1],derivatives[,1],Losses[,1])),1,function(x) paste0(x,collapse="")),"]")
standfrag<-cbind(a,b)
colnames(standfrag)<-c("m/z","ID")
standfrag
}
#Use matchppm2 to match m/z columns in the data.
matchppm2(stan,as.numeric(s22s4DmgPxf[,1]),1,1)
#
#
#
###############################################################################################

1. WORKFLOWS for standard spectral expounding: ID all peaks.
2. WORKFLOWS for NAFC spectrum expounding: ID all peaks and classes
3. OVERLAPS: contaminations that fall into NAFC classes.
4. 

#formulamtch2 workflow
fmfull<-function(data){
data<-formulamatch2(data,newp,1,1)
data<-formulamatch2(data,newp,1,2)
data<-formulamatch2(data,stan,1,1)
data<-formulamatch2(data,contam,1,1)
data
}

fmfulln<-function(data){
data<-formulamatch2(data,newn,1,1)
data<-formulamatch2(data,newn,1,2)
#data<-formulamatch2(data,stan,1,1)
data<-formulamatch2(data,contamn,1,1)
data
}



#unlist
#library(gplots)
#venn()
library(gplots)
bob2e3<-list(list())
range<-seq(72.1020,72.1060,by=0.0001)
for(i in 1:length(range)){
bob2e3[[i]]<-venn(list(pos=(s21s2na[s21s2na[,2]<500 & s21s2na[,2]>100 & s21s2na[,22]>2e3,2]-range[i]),neg=s16s14na[s16s14na[,2]<400 & s16s14na[,22]>2e3,2]))
}
bob2e3
