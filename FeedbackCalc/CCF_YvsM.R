#
#  Calculate and represent the cross covariance functions between M and Y
#
#   by Albert Ossó ,2024
#***********************************************

library('abind');
library('ncdf4');
library('zoo');
library('MASS');
library("TSA");
library("DescTools");
library("signal");
library('splus2R');

# Substitue by your own
#DIR="***"; 
#DIR2="***";


lag.l=35;  # chose number of lags
seas="JJA";

# Allocate arrays
lag=array(0,dim=c(33,lag.l*2+1));
ccf=array(0,dim=c(33,lag.l*2+1));
lag.wholebox=array(0,dim=c(33,lag.l*2+1));
ccf.wholebox=array(0,dim=c(33,lag.l*2+1));
lag.southbox=array(0,dim=c(33,lag.l*2+1));
ccf.southbox=array(0,dim=c(33,lag.l*2+1));
lag.northbox=array(0,dim=c(33,lag.l*2+1));
ccf.northbox=array(0,dim=c(33,lag.l*2+1));
acf.values=array(0,dim=c(33,35));
Area=array(0,dim=c(33));
Area.acf=array(0,dim=c(33));
lagMax=array(0,dim=c(33));


models=c("ACCESS-ESM1-5","AWI-ESM-1-1-LR","BCC-CSM2-MR","BCC-ESM1","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2","CESM2-WACCM","CMCC-CM2-HR4","CMCC-CM2-SR5","CMCC-ESM2","EC-Earth3-AerChem","EC-Earth3-CC","EC-Earth3","EC-Earth3-Veg-LR","EC-Earth3-Veg","FGOALS-g3","IITM-ESM","INM-CM4-8","INM-CM5-0","IPSL-CM5A2-INCA","IPSL-CM6A-LR","MIROC6","MPI-ESM-1-2-HAM","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","SAM0-UNICON","TaiESM1","era5");

nmodels=length(models);
efoldingMfor_lag=array(0,dim=nmodels);
efoldingCCF_lag=array(0,dim=nmodels);

boxes=c("wholebox","northbox","southbox");
for(box in boxes){  # Loop across boxes
count=1;
for(m in models){    # Loop across models

#********************************************OPEN M and Y TIMESERIES and chose wave or eddy*********************************
titlefile=sprintf("%s/Mforcing_Wave_%s_%s_1980-2014_%s_Barnesbox.txt",DIR,m,seas,box);
load(titlefile);   
Mfor=Mforcing;

# Chose for eddy forcing
#titlefile=sprintf("%s/Mforcing_%s_%s_1980-2014_%s_HIGHPASS_Barnesbox.txt",DIR,m,seas,box);
#load(titlefile);   
#Mfor=Mforcing;

titlefile=sprintf("%s/YNAO_%s_%s_1980-2014_wholebox_Barnesbox.txt",DIR,m,seas,box);
load(titlefile);   
YNAO=YNAO;
time=seq(1,length(YNAO),by=1);

#**************************Calculate the CCF***************************************
# 1. Remove linear trends
kitrend=YNAO;
kjtrend=time;
kwitrend=kitrend[!is.na(kitrend)];
kwjtrend=kjtrend[!is.na(kitrend)];
lm_fit=lm(kwitrend~kwjtrend);
predict_fit=predict(lm_fit);
YNAO_det=kwitrend-predict_fit;

kitrend=Mfor;
kjtrend=time;
kwitrend=kitrend[!is.na(kitrend)];
kwjtrend=kjtrend[!is.na(kitrend)];
lm_fit=lm(kwitrend~kwjtrend);
predict_fit=predict(lm_fit);
Mfor_det=kwitrend-predict_fit;
#**********************************************************************

CCF=ccf(Mfor_det,YNAO_det,lag=lag.l,plot=FALSE);
lag[count,]=CCF$lag;
ccf[count,]=CCF$acf;

#Calculate autocorrelation function for Mfor  
acf_result=acf(Mfor,plot=FALSE);
acf.values[count,]=acf_result$acf;
signi.acf.Mfor=qnorm((1+0.95)/2)/sqrt(sum(!is.na(Mfor)));

sig_up=2/sqrt(length(Mfor));
sig_down=-2/sqrt(length(Mfor));

count=count+1;
}

print(box);
iSf(box=="wholebox"){
lag.wholebox=lag;
ccf.wholebox=ccf;
}
if(box=="northbox"){
lag.northbox=lag;
ccf.northbox=ccf;
}
if(box=="southbox"){
lag.southbox=lag;
ccf.southbox=ccf;
}


}  #loop boxes

 count=count-1;
  
 #*******wholebox plot*********
ccf=ccf.wholebox;
lag=lag.wholebox 
q1=array(0,dim=lag.l*2+1);
q50=array(0,dim=lag.l*2+1);
q90=array(0,dim=lag.l*2+1);
 
 for(i in seq(1,lag.l*2+1)){   # loop across all lags 
 
 q1[i]=quantile(ccf[,i],probs=0.10);
 q50[i]=quantile(ccf[,i],probs=0.5);
 q90[i]=quantile(ccf[,i],probs=0.90);
  
 }
  
lagPol=c(-lag.l:lag.l,lag.l:-lag.l);
quantPol=c(q90,rev(q1)); 
 
	plotname=sprintf("%s/CCF_Mfor.wave-Ywholebox_%s_1980-2014_wholebox_HIGHPASS_Barnesbox.pdf",DIR,seas);
	pdf(plotname,width=7,height=5);
	plot(lag[32,],ccf[32,],type="l",lwd=3,col="black",ylim=c(-0.2,0.4),xlab="lag(days)",ylab="ccf");
	polygon(lagPol,quantPol,density=20,col="orangered");
	lines(lag[32,],q50,type="l",lty=1,lwd=1.5,col="orangered");
	abline(h=0,col="black",lty=1,lwd=2);
	abline(v=0,col="black",lty=1,lwd=2);
	abline(h=sig_up,col="blue",lty=2,lwd=1.5);
	abline(h=sig_down,col="blue",lty=2,lwd=1.5);
	dev.off();

 #*******southbox plot*********
ccf=ccf.southbox;
lag=lag.southbox 
q1=array(0,dim=lag.l*2+1);
q50=array(0,dim=lag.l*2+1);
q90=array(0,dim=lag.l*2+1);
 
 for(i in seq(1,lag.l*2+1)){   # loop across all lags
 
 q1[i]=quantile(ccf[,i],probs=0.10);
 q50[i]=quantile(ccf[,i],probs=0.5);
 q90[i]=quantile(ccf[,i],probs=0.90);
  
 }
  
lagPol=c(-lag.l:lag.l,lag.l:-lag.l);
quantPol=c(q90,rev(q1)); 
 
	plotname=sprintf("%s/CCF_Mfor.wave-Ywholebox_%s_1980-2014_southbox_HIGHPASS_Barnesbox3.pdf",DIR,seas);
	pdf(plotname,width=7,height=5);
	plot(lag[32,],ccf[32,],type="l",lwd=3,col="black",ylim=c(-0.2,0.4),xlab="lag(days)",ylab="ccf");
	polygon(lagPol,quantPol,density=20,col="orangered");
	lines(lag[32,],q50,type="l",lty=1,lwd=1.5,col="orangered");
	abline(h=0,col="black",lty=1,lwd=2);
	abline(v=0,col="black",lty=1,lwd=2);
	abline(h=sig_up,col="blue",lty=2,lwd=1.5);
	abline(h=sig_down,col="blue",lty=2,lwd=1.5);
	dev.off();

 #*******northbox plot*********
ccf=ccf.northbox;
lag=lag.northbox 
q1=array(0,dim=lag.l*2+1);
q50=array(0,dim=lag.l*2+1);
q90=array(0,dim=lag.l*2+1);
 
 for(i in seq(1,lag.l*2+1)){   # loop across all lags
 
 q1[i]=quantile(ccf[,i],probs=0.10);
 q50[i]=quantile(ccf[,i],probs=0.5);
 q90[i]=quantile(ccf[,i],probs=0.90);
  
 }
  
lagPol=c(-lag.l:lag.l,lag.l:-lag.l);
quantPol=c(q90,rev(q1)); 
 
	plotname=sprintf("%s/CCF_Mfor.wave-Ywholebox_%s_1980-2014_northbox_HIGHPASS_Barnesbox3.pdf",DIR,seas);
	pdf(plotname,width=7,height=5);
	plot(lag[32,],ccf[32,],type="l",lwd=3,col="black",ylim=c(-0.2,0.4),xlab="lag(days)",ylab="ccf");
	polygon(lagPol,quantPol,density=20,col="orangered");
	lines(lag[32,],q50,type="l",lty=1,lwd=1.5,col="orangered");
	abline(h=0,col="black",lty=1,lwd=2);
	abline(v=0,col="black",lty=1,lwd=2);
	abline(h=sig_up,col="blue",lty=2,lwd=1.5);
	abline(h=sig_down,col="blue",lty=2,lwd=1.5);
	dev.off();


