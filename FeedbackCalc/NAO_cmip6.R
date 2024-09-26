#********************************************
#  Calculate the First EOF and PC of the NATL SLP using eof.mca.R and cov4gappy.R
#
#   by  Albert Ossó, 2024
#*************************************************
library('abind');
library('ncdf4');
library('MASS');
source("eof.mca.R");
source("cov4gappy.R");
source("trend_AO.R");  

# substitue by your own directories
#DIR="******";
#DIR2="******";  
#DIR3="******";
#DIR4="******";

seas="JJA";

models=c("ACCESS-ESM1-5","AWI-ESM-1-1-LR","BCC-CSM2-MR","BCC-ESM1","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2","CESM2-WACCM","CMCC-CM2-HR4","CMCC-CM2-SR5","CMCC-ESM2","EC-Earth3-AerChem","EC-Earth3-CC","EC-Earth3","EC-Earth3-Veg-LR","EC-Earth3-Veg","FGOALS-g3","IITM-ESM","INM-CM4-8","INM-CM5-0","IPSL-CM5A2-INCA","IPSL-CM6A-LR","MIROC6","MPI-ESM-1-2-HAM","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","SAM0-UNICON","TaiESM1","era5");

nmodels=length(models);


for(m in models){   # Loop across models. 

file.day=sprintf("%s/slpanom_%s_day_%s.nc",DIR2,seas,m);   # daily slp data 
file.mon=sprintf("%s/slpanom_%s_monmean_%s.nc",DIR2,seas,m); # monthly slp data

#*******Open and read file*********************
title=file.mon;
slp.nc=nc_open(title,write=FALSE);
slp=ncvar_get(slp.nc,"psl"); 
slp=slp/100;
lat=ncvar_get(slp.nc,"lat"); 
lon=ncvar_get(slp.nc,"lon");  
time=ncvar_get(slp.nc,"time");  

slp_e=slp[1:73,,];    
slp_w=slp[74:144,,];
slp=abind(slp_w,slp_e,along=1);

durat=length(time);
lon=seq(-177.5,180,2.5);
#*************************************************

#******* Select NATL area to calculate the SNAO**********
slp_area=slp[36:84,11:29,];    
lat_area=lat[11:29];      
lon_area=lon[36:84];      # Barnes SNAO 25 to 70 and 90 to 30  
#*********************************************

slp_area=matrix(slp_area,durat,length(lon_area)*length(lat_area),byrow=TRUE); #Arrange data in S mode  

#*********************cos latitude weighting*********************************************************

vpesos=cos(lat_area*pi/180.0);  		# weigth vector
vpesos=rep(vpesos,each=length(lon_area)); 	# map of weights
vpesos=rep(vpesos,each=durat); 
mpesos=matrix(vpesos,durat,length(lon_area)*length(lat_area)); #weight matrix

slp_area=slp_area*mpesos;  	
rm(vpesos,mpesos);

#********Calculate the EOFS using the R function eof.mca************

mca=eof.mca(slp_area,nu=1,nv=1,centered=TRUE,scaled=FALSE);

#****** Open and read daily SLP timeseries************
title=fileout.day;
slp.nc=nc_open(title,write=FALSE);
slp.day=ncvar_get(slp.nc,"psl"); 
slp.day=slp.day/100;
lat=ncvar_get(slp.nc,"lat"); 
lon=ncvar_get(slp.nc,"lon");  
time.day=ncvar_get(slp.nc,"time");  
durat=length(time.day);

slp_e=slp.day[1:73,,];     
slp_w=slp.day[74:144,,];

slp.day=abind(slp_w,slp_e,along=1);
lon=seq(-177.5,180,2.5);
#************************************************
#*********Select NATL area************** 
slp.day_area=slp.day[36:84,11:29,];    
lat_area=lat[11:29];      
lon_area=lon[36:84];      # Barnes SNAO 25 to 70 and 90 to 30  
#*****************************************************

slp.dayM=matrix(slp.day_area,durat,length(lon_area)*length(lat_area),byrow=TRUE);  #Arrange in S mode 
slp.dayM=replace(slp.dayM,is.na(slp.dayM),0);   # Remove NANs

EOF1V=as.vector(mca$u);
EOF1V=replace(EOF1V,is.na(EOF1V),0);

#******Project daily slp anomalies onto the monthly NAO EOF1 ************
NAOI=slp.dayM%*%EOF1V;    
sdM=sd(NAOI);
NAOI=(NAOI-mean(NAOI))/sd(NAOI);  # standardize time series 
NAOI=NAOI[,1];
titlefile=sprintf("%s/NAOI_%s_%s_1980-2014_SNAO_Barnesbox.txt",DIR4,seas,m);
save(NAOI,file=titlefile);
#*********************************************************

#*******Extract PCs********************************
PCslp=mca$A[,1];
PC1Var=mca$expl_var;
sdPC1=sd(PCslp);
PCslp=as.vector(PCslp);
PCslp_norm=PCslp/sd(PCslp);
titlefile=sprintf("%s/PC1.mon_slp_%s_%s_1980-2014_NAO_Sepbox.txt",DIR4,m,seas);
save(sdPC1,PC1Var,PCslp_norm,file=titlefile);


#*********************************Regression analysis for representation of SNAO***********************

	print("Regression analysis");
	
	pvalue_slp=matrix(0,length(lon),length(lat));
	reg_slp=matrix(0,length(lon),length(lat));
	regslp_sig=matrix(0,length(lon),length(lat));


  for(i in seq(1,length(lon),by=1)){
	    for(j in seq(1,length(lat),by=1)){
			
		kitrend=slp[i,j,];
		kjtrend=PCslp_norm;
		trend=trend_AO(kitrend,kjtrend);
		reg_slp[i,j]=trend$b;		
		pvalue_slp[i,j]=trend$sig;
		logi_sig=is.na(pvalue_slp[i,j]);
		
	if(logi_sig=="TRUE"){
	regslp_sig[i,j]=NA;	
	}else{
	
	if(pvalue_slp[i,j]==1){
	regslp_sig[i,j]=reg_slp[i,j];
	}else{
	regslp_sig[i,j]=NA;
	}
	}
		
	     }
	}


#******************************

londim=ncdim_def("lon","degrees_east",lon); 
latdim=ncdim_def("lat","degrees_north",lat); 

 # define variables
fillvalue=1e32;

dlname="reg_slp";
reg_slp.def=ncvar_def("regslp","nounits",list(londim,latdim),fillvalue,dlname,prec="single");

dlname="reg_slpsig";
reg_slpsig.def=ncvar_def("regslpsig","nounits",list(londim,latdim),fillvalue,dlname,prec="single");

dlname="pvalue_slp";
pvalue_slp.def=ncvar_def("pvalslp","nounits",list(londim,latdim),fillvalue,dlname,prec="single");


# create netCDF file and put arrays  
nctitle1=sprintf("%s/reg_slp_PC1slpMon_%s_%s_80-2014_NAO_Barnesbox.nc",DIR,m,seas);
ncfname=nctitle1;
ncout=nc_create(ncfname,list(reg_slp.def,reg_slpsig.def,pvalue_slp.def));

# put variables
ncvar_put(ncout,reg_slp.def,reg_slp);
ncvar_put(ncout,reg_slpsig.def,regslp_sig);
ncvar_put(ncout,pvalue_slp.def,pvalue_slp);

# put additional attributes into dimension and data variables
ncatt_put(ncout,"p","axis","Z"); #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y");
ncatt_put(ncout,"lon","axis","X");

# close the file, writing data to disk
nc_close(ncout);

}





