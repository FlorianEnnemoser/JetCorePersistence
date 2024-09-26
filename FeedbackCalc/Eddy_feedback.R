#
#  Calculate the eddy feedback strenght following Simpson et al. 2013 and Barnes and Hartmann 2010
#  Check the correlation between Y and NAOI persitence and plots
#  Boxplot of Y and NAOI persistence
#  Calculate autocorrelations and plot
#  Calculate the scatter plot between eddy feedback and Y persistence and plot 
#       by Albert Ossó, 2024
#***************************************************
library('abind');
library('ncdf4');
library('zoo');
library('MASS');
library("TSA");
library("DescTools");
library("dplyr");

# substitude by your own directories
#DIR="****";
#DIR2="****";

num.mod=33;
nlags=35;
seas="JJA";
#box="northbox";
box="southbox";
#box="wholebox";

# Allocate arrays*********************
efoldingNAO=array(0,dim=num.mod);
efoldingYNAO=array(0,dim=num.mod);
efoldingMfor=array(0,dim=num.mod);

q1=array(0,dim=nlags);
q50=array(0,dim=nlags);
q99=array(0,dim=nlags);

b=array(0,dim=c(num.mod,nlags));

acf.values.NAO=array(0,dim=c(num.mod,nlags));
q1.NAO=array(0,dim=nlags);
q50.NAO=array(0,dim=nlags);
q99.NAO=array(0,dim=nlags);

acf.values.YNAO=array(0,dim=c(num.mod,nlags));
q1.YNAO=array(0,dim=nlags);
q50.YNAO=array(0,dim=nlags);
q99.YNAO=array(0,dim=nlags);

acf.values.Mfor=array(0,dim=c(num.mod,nlags));
q1.Mfor=array(0,dim=nlags);
q50.Mfor=array(0,dim=nlags);
q99.Mfor=array(0,dim=nlags);

#***********
models=c("ACCESS-ESM1-5","AWI-ESM-1-1-LR","BCC-CSM2-MR","BCC-ESM1","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2","CESM2-WACCM","CMCC-CM2-HR4","CMCC-CM2-SR5","CMCC-ESM2","EC-Earth3-AerChem","EC-Earth3-CC","EC-Earth3","EC-Earth3-Veg-LR","EC-Earth3-Veg","IITM-ESM","INM-CM4-8","INM-CM5-0","IPSL-CM5A2-INCA","IPSL-CM6A-LR","MIROC6","MPI-ESM-1-2-HAM","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","SAM0-UNICON","TaiESM1","era5");

nmodels=length(models);
count=1;

for(m in models){  # Loop for all models

#********************************************OPEN M and Y TIMESERIES*********************************

titlefile=sprintf("%s/Mforcing_%s_%s_1980-2014_%s_HIGHPASS_Barnesbox.txt",DIR,m,seas,box);
load(titlefile);   
Mfor=Mforcing;

#titlefile=sprintf("%s/Mforcing_Wave_%s_%s_1980-2014_%s_Barnesbox.txt",DIR,m,seas,box);
#load(titlefile);   
#Mfor.wave=Mforcing;

#Mfor=Mfor.wave;

titlefile=sprintf("%s/YNAO_%s_%s_1980-2014_wholebox_Barnesbox.txt",DIR,m,seas,box);
load(titlefile);   
YNAO=YNAO;

titlefile=sprintf("%s/NAOI_%s_%s_1980-2014_SNAO_Barnesbox.txt",DIR,seas,m);
load(titlefile);   
NAO=NAOI;


#*******************YNAO,SNAO and Mfor autocorrelation and persistence****************************************************
acf_result=acf(YNAO,lag.max=35,plot=FALSE);
acf.values.YNAO[count,]=acf_result$acf;
acf.lags=acf_result$lag;
signi.YNAO=qnorm((1+0.95)/2)/sqrt(sum(!is.na(YNAO)));

# Find the threshold value for e-folding
threshold=1/exp(1);

# Find the first lag where ACF falls below the thresold
efoldingYNAO[count]=which(acf.values.YNAO[count,]<=threshold*acf.values.YNAO[1])[1];

acf_result=acf(NAO,lag.max=35,plot=FALSE);
acf.values.NAO[count,]=acf_result$acf;
acf.lags=acf_result$lag;
signi.NAO=qnorm((1+0.95)/2)/sqrt(sum(!is.na(NAO)));
efoldingNAO[count]=which(acf.values.NAO[count,]<=threshold*acf.values.NAO[1])[1];

acf_result=acf(Mfor,lag.max=35,plot=FALSE);
acf.values.Mfor[count,]=acf_result$acf;
acf.lags=acf_result$lag;
signi.Mfor=qnorm((1+0.95)/2)/sqrt(sum(!is.na(Mfor)));

efoldingMfor[count]=which(acf.values.Mfor[count,]<=threshold*acf.values.Mfor[1])[1];

#**************************eddy-feedback parameter***************************
#Calculate lagged regression coefficients as in Simpson et al.
for(i in seq(1,nlags,by=1)){
nlag=i;
YNAO_lead=lead(YNAO,nlag);

YNAO_lead=YNAO_lead[!is.na(YNAO_lead)];
Mfor_align=Mfor[(nlag+1):length(Mfor)];

lm.model=lm(Mfor_align~YNAO_lead);
b[count,i]=lm.model$coefficients[2];

}
#*******************************************************************
count=count+1;
}

efoldingYNAO[is.na(efoldingYNAO)]=nlags;
efoldingNAO[is.na(efoldingNAO)]=nlags;
efoldingMfor[is.na(efoldingMfor)]=nlags;

#*****************Autocorrelation Plots***********************
lags=seq(1,35);

for(i in seq(1,35,by=1)){   # loop for all lags
q1.NAO[i]=quantile(acf.values.NAO[c(1:(num.mod-1)),i],probs=0.01);
q50.NAO[i]=quantile(acf.values.NAO[c(1:(num.mod-1)),i],probs=0.50);
q99.NAO[i]=quantile(acf.values.NAO[c(1:(num.mod-1)),i],probs=0.99);
}

for(i in seq(1,35,by=1)){   # loop for all lags
q1.YNAO[i]=quantile(acf.values.YNAO[c(1:(num.mod-1)),i],probs=0.01);
q50.YNAO[i]=quantile(acf.values.YNAO[c(1:(num.mod-1)),i],probs=0.50);
q99.YNAO[i]=quantile(acf.values.YNAO[c(1:(num.mod-1)),i],probs=0.99);
}

for(i in seq(1,35,by=1)){   # loop for all lags
q1.Mfor[i]=quantile(acf.values.Mfor[c(1:(num.mod-1)),i],probs=0.01);
q50.Mfor[i]=quantile(acf.values.Mfor[c(1:(num.mod-1)),i],probs=0.50);
q99.Mfor[i]=quantile(acf.values.Mfor[c(1:(num.mod-1)),i],probs=0.99);
}

lagPol=c(1:35,35:1);
quantPol.NAO=c(q99.NAO,rev(q1.NAO)); 
quantPol.YNAO=c(q99.YNAO,rev(q1.YNAO)); 
quantPol.Mfor=c(q99.Mfor,rev(q1.Mfor)); 

#*****************Autocorrelation plots*******************
#plotname=sprintf("%s/Autocorrelation_%s_YNAOvsMfor_%s_Barnesbox.pdf",DIR,seas,box);
#	plottitle=sprintf("r=%3.2f pvalue=%3.2f",cor,pvalue);
#	pdf(plotname,width=7,height=5);
#	par(las=1,cex.main=1.0);    #ylim=c(224,232),

#	plot(acf.values.YNAO[32,]~lags,type="l",col="blue",lwd=2,ylim=c(-0.2,1));
#	lines(lags,q99.YNAO,lty=2,col="blue");
#	lines(lags,q1.YNAO,lty=2,col="blue");

#	lines(lags,acf.values.Mfor[32,],type="l",col="red");
#	lines(lags,q99.Mfor,lty=2,col="red");
#	lines(lags,q1.Mfor,lty=2,col="red");

#	abline(h=signi.Mfor,col='red',lty=2);
#	abline(h=signi.YNAO,col='blue',lty=2);
#	abline(h=0,col='black',lty=1);

#	dev.off();


#plotname=sprintf("%s/Autocorrelation_%s_YNAOvsNAO_%s_Barnesbox.pdf",DIR,seas,box);
##	plottitle=sprintf("r=%3.2f pvalue=%3.2f",cor,pvalue);
#	pdf(plotname,width=7,height=5);
#	par(las=1,cex.main=1.0);    #ylim=c(224,232),

#	plot(acf.values.YNAO[32,]~lags,type="l",col="blue",lwd=2,ylim=c(-0.2,1));
#	lines(lags,q99.YNAO,lty=2,col="blue");
#	lines(lags,q1.YNAO,lty=2,col="blue");

#	lines(lags,acf.values.NAO[32,],type="l",col="red",lwd=2);
#	lines(lags,q99.NAO,lty=2,col="red");
#	lines(lags,q1.NAO,lty=2,col="red");

#	abline(h=signi.NAO,col='red',lty=2);
#	abline(h=signi.YNAO,col='blue',lty=2);
#	abline(h=0,col='black',lty=1);

#	dev.off();


#**************YNAO vs NAOI to check model validity*****************

	corre=cor.test(efoldingNAO,efoldingYNAO);
	pvalue=corre$p.value;
	cor=corre$estimate;
	lm.model=lm(efoldingNAO~efoldingYNAO);

	ymax=max(efoldingNAO)+2;
	ymin=min(efoldingNAO)-2;

plotname=sprintf("%s/efolding_YNAOwholevsNAO_%s_%s_Barnesbox.pdf",DIR,seas,box);
	plottitle=sprintf("r=%3.2f pvalue=%5.4f",cor,pvalue);
	pdf(plotname,width=7,height=5);
	par(las=1,cex.main=1.0);    #ylim=c(224,232),
	plot(jitter(efoldingNAO[1:31])~jitter(efoldingYNAO[1:31]),ylim=c(ymin,ymax),type="p",main=plottitle);
	points(efoldingNAO[32]~efoldingYNAO[32],pch=19,col="red",cex=2);
      axis(1,at=seq(3,6,by=1),labels=c(3,4,5,6));
	abline(lm.model,col='red',lwd=2);
	dev.off();


bands.M=rep(c(1,2),each=(num.mod-1));
efolding.M=c(efoldingNAO[1:(num.mod-1)],efoldingYNAO[1:(num.mod-1)]);

efolding.NAO.mean=mean(efoldingNAO[1:(num.mod-1)]);
efolding.YNAO.mean=mean(efoldingYNAO[1:(num.mod-1)]);

bands.M.era=c(1,2);
efolding.M.era=c(efoldingNAO[num.mod],efoldingYNAO[num.mod]);


plotname=sprintf("%s/boxplot_efolding_YNAOwholevsNAO_%s_%s_Barnesbox.pdf",DIR,seas,box);
	plottitle=sprintf("r=%3.2f pvalue=%3.2f",cor,pvalue);
	pdf(plotname,width=7,height=5);
	par(las=1,cex.main=1.0);    #ylim=c(224,232),
	boxplot(jitter(efolding.M)~bands.M,col="white",boxwex=0.5);
	points(jitter(efolding.M.era)~bands.M.era,pch=19,col="red",cex=2);
	abline(h=efolding.YNAO.mean,col="red");    
	abline(h=efolding.NAO.mean,col="green");    
	dev.off();



#**************************Feedback strength as a function of lag plot**************
for(i in seq(1,nlags,by=1)){
q1[i]=quantile(b[c(1:(num.mod-1)),i],probs=0.01);
q50[i]=quantile(b[c(1:(num.mod-1)),i],probs=0.50);
q99[i]=quantile(b[c(1:(num.mod-1)),i],probs=0.99);
}

ymax=max(q99)+0.01;
ymin=min(q1)-0.01;
lag=seq(1,nlags,by=1);

plotname=sprintf("%s/feedback_lag1to35_%s_%s_Barnesbox.pdf",DIR,seas,box);
	pdf(plotname,width=7,height=5);
	par(las=1,cex.main=1.0);    #ylim=c(224,232),
        plot(b[32,]~lag,type="l",lty=1,ylim=c(ymin,ymax));
	lines(lag,q99,type="l",lty=2,col="red");
	lines(lag,q1,type="l",lty=2,col="blue");

	dev.off();


#*************************************Persistence vs average feedback strength*********

bmean=rowMeans(b[,c(5:35)],na.rm=TRUE);  # average the feedback paramater from 5 to 35 lag
ymin=min(efoldingYNAO)-2;
ymax=max(efoldingYNAO)+2;
xmin=min(bmean)-0.1;
xmax=max(bmean)+0.1;

	corre=cor.test(bmean,efoldingYNAO);
	pvalue_bo=corre$p.value;
	cor=corre$estimate;
	lm.model=lm(efoldingYNAO~bmean);

plotname=sprintf("%s/efoldingYNAOwhole_feeback.bmean1:35_%s_%s_HIGHPASS_Barnesbox.pdf",DIR,seas,box);
	plottitle=sprintf("r=%3.2f pvalue=%3.2f",cor,pvalue_bo);
	pdf(plotname,width=7,height=5);
	 plot(jitter(bmean[1:(num.mod-1)]),jitter(efoldingYNAO[1:(num.mod-1)]),
	    main=plottitle,ylim=c(ymin,ymax),xlim=c(xmin,xmax),type="p",col="black");
	points(jitter(bmean[num.mod]),jitter(efoldingYNAO[num.mod]),type="p",pch=16,cex=1.5,col="green");
	abline(lm.model,col='red',lwd=2);
	abline(h=0);
	dev.off();



