library(mzR)
library(plyr)
#Function for unpacking a list of m/z & peak intensity matrices from a chromatography-mass spectrometry run and binning
#the data according to specified 2-D bins. One can bin according to m/z and time. One can also choose to bin "by a multiple"
#(mzby, tby) or choose the bin size outright (mzbins, tbins), or a combination of the two.
#or provide m/z vector. Example: mzp2<-sort(unique(round(unique(unlist(lapply(s22s4DmgP,function(x) x[x[,1]<600,1]))),4)))
mzUnpackBin2<-function(list,ret_vec,mzbins=NULL,tbins=NULL,mzby=NULL,tby=NULL,mzprovide=NULL){
require(plyr)
list <- as.list(list)
ret_vec<-as.numeric(ret_vec)

#If the length argument is specified
if(is.null(mzprovide)){
alpha_mz<-if(is.null(mzbins)){c(seq(from=(range(ldply(list,unlist)[,1])[1]),to=(range(ldply(list,unlist)[,1])[2]),by=mzby),(range(ldply(list,unlist)[,1])[2]))} else {seq(from=range(ldply(list,unlist)[,1])[1],to=range(ldply(list,unlist)[,1])[2],length=mzbins)}
alpha_t<-if(is.null(tbins)){c(seq(from=ret_vec[1],to=ret_vec[length(ret_vec)],by=tby),ret_vec[length(ret_vec)])} else {seq(from=ret_vec[1],to=ret_vec[length(ret_vec)],length=tbins)}
} else {
alpha_mz<-mzprovide
alpha_t<-if(is.null(tbins)){c(seq(from=ret_vec[1],to=ret_vec[length(ret_vec)],by=tby),ret_vec[length(ret_vec)])} else {seq(from=ret_vec[1],to=ret_vec[length(ret_vec)],length=tbins)}}

#Summation function (to start with).#
####Can add in the option of mean or specified regression function.###
##Definitions of emptymat and alpha_mz depend on "discrete" or "continuous" options

emptymat<-matrix(0,ncol=length(alpha_t),nrow=length(alpha_mz))


for (i in 2:length(alpha_t)){
for (j in 2:length(alpha_mz)){
#Instead of a match, a </=/> in the index will flag certain values.
#First is to bin intensities into the appropriate m/z bin.

new <- sum(unlist(lapply(list[c(ret_vec<=alpha_t[i] & ret_vec>=alpha_t[i-1])],function(x) x[c(x[,1]<=alpha_mz[j] & x[,1]>=alpha_mz[j-1]),2])))
emptymat[j,i] <- new

}
emptymat<-emptymat
}
emptymat<-emptymat
#Result is a list of data matrix, and bin values for m/z and time.
listf<-list(emptymat,alpha_mz,alpha_t)
}



###########################
#script#
###########################

mzmllist<-list.files("D:\\20160919 Orbitrap Negative Ion Mode\\KevinMatthew_PosIon_Sept22_2016",pattern=".mzML",full.names=T)
#mzbyvec<-c(1)
#tbinsvec<-c(1)
for(i in seq(along=mzmllist)){
#	for(j in seq(along=mzbyvec)){
#		for (k in seq(along=tbinsvec)){
			open<-openMSfile(mzmllist[i])
			peak<-peaks(open)
			peakx<-lapply(peak,function(x) x[x[,1]<600,])
			peakx<-lapply(peakx,function(x) matrix(x,ncol=2,byrow=F))
			mzp<-sort(unique(round(unique(unlist(lapply(peak,function(x) x[x[,1]<600,1]))),4)))	
			header<-header(open)
			close(open)
			binned<-mzUnpackBin2(peakx,header$ret,tby=1,mzprovide=mzp)
			write.csv(rbind(c(0,binned[[3]]),cbind(binned[[2]],binned[[1]])),file=paste0(mzmllist[i],".csv"))
			rm(peak,header)
#		}

#	}


}


mzmllist<-list.files("D:\\20160919 Orbitrap Negative Ion Mode\\20160921_KevinMatthew_PosIon_Sept21_2016",pattern=".mzML",full.names=T)
#mzbyvec<-c(1)
#tbinsvec<-c(1)
for(i in seq(along=mzmllist)){
#	for(j in seq(along=mzbyvec)){
#		for (k in seq(along=tbinsvec)){
			open<-openMSfile(mzmllist[i])
			peak<-peaks(open)
			peakx<-lapply(peak,function(x) x[x[,1]<600,])
			peakx<-lapply(peakx,function(x) matrix(x,ncol=2,byrow=F))
			mzp<-sort(unique(round(unique(unlist(lapply(peak,function(x) x[x[,1]<600,1]))),4)))	
			header<-header(open)
			close(open)
			binned<-mzUnpackBin2(peakx,header$ret,tby=1,mzprovide=mzp)
			write.csv(rbind(c(0,binned[[3]]),cbind(binned[[2]],binned[[1]])),file=paste0(mzmllist[i],".csv"))
			rm(peak,header)
#		}

#	}


}

