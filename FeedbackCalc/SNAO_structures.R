#*******************************
# Project the vorticity terms onto the NAOI timeseries for preparing Fig.2.
# By Albert Ossó 2024.

#*******************************************

library('abind');
library('ncdf4');
library('zoo');
library('MASS');
library("TSA");

DIR="/nas/home/aos/Desktop/JetPersistence";
DIR2="/data/users/aos/00_BARNS_HARTMANN_PAPER_save";

level="250"; # Chose vertical level
seas="JJA";  # Chose season  
#**********Chose parameter********************* 

#term="VORTICITY";
#term.path="vorticity";
#var.name="rel_vort";

#term="EDDY";
#term.path="eddy";
#var.name="eddy";     # Change filein and final file for synoptic EDDY

#term="STRECHING"
#term.path="streching";
#var.name="streching";

term="WAVE"
term.path="wave";
var.name="wave";
#**************************************************
models=c("ACCESS-ESM1-5","AWI-ESM-1-1-LR","BCC-CSM2-MR","BCC-ESM1","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2","CESM2-WACCM","CMCC-CM2-HR4","CMCC-CM2-SR5","CMCC-ESM2","EC-Earth3-AerChem","EC-Earth3-CC","EC-Earth3","EC-Earth3-Veg-LR","EC-Earth3-Veg","FGOALS-g3","IITM-ESM","INM-CM4-8","INM-CM5-0","IPSL-CM5A2-INCA","IPSL-CM6A-LR","MIROC6","MPI-ESM-1-2-HAM","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","SAM0-UNICON","TaiESM1","era5");

nmodels=length(models);
count=1;

for(m in models){   # Loop across models

#********************************************OPEN NAOI TIMESERIES*********************************
titlefile=sprintf("/nas/home/aos/Desktop/JetPersistence/newData/NAOI_%s_%s_1980-2014_SNAO_Barnesbox.txt",seas,m);
load(titlefile);   
NAOi=NAOI;
time=seq(1,length(NAOi),by=1);


#*****************************OPEN DYN VARIABLE*******************************************************

# Switch for HIGHPASS or not HIGHPASS********
filein=sprintf("%s/%s/%s_%s.nc",DIR2,term.path,term,m);
#filein=sprintf("%s/%s/%s_%s_HIGHPASS.nc",DIR2,term.path,term,m);
#*******************************CALCULATE ANOMALIES USING CDO OPERATORS****************
fileout=sprintf("%s/merda/vort_%s.nc",DIR2,m);
command=sprintf("cdo seldate,1980-1-1,2014-11-30 %s %s",filein,fileout);
system(command);

filein=fileout;
fileout=sprintf("%s/merda/VORTICITY_250_%s.nc",DIR2,m);
command=sprintf("cdo sellevel,%s %s %s",level,filein,fileout);
system(command);

filein=fileout;
fileout1=sprintf("%s/merda/VORTICITY_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo selseas,%s %s %s",seas,filein,fileout1);
system(command);

filein=fileout1;
fileout2=sprintf("%s/merda/vort_seasm_%s.nc",DIR2,m);
command=sprintf("cdo yearmean %s %s",filein,fileout2);
system(command);

filein=fileout2;
fileout3=sprintf("%s/merda/vortclim_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo timmean %s %s",filein,fileout3);
system(command);

fileout=sprintf("%s/merda/vortanomm_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo -b 64 sub %s %s %s",fileout1,fileout3,fileout);
system(command);

#*************************READ file***********
title=fileout;
u.nc=nc_open(title,write=FALSE);
vort=ncvar_get(u.nc,var.name);  
lat=ncvar_get(u.nc,"latitude");    # change for era5 mod file
lon=ncvar_get(u.nc,"longitude");
time.data=ncvar_get(u.nc,"time");
time.vort=seq(1,length(time.data),by=1);

command=sprintf("rm -f %s/merda/*.nc",DIR2);    
system(command);

#**********************Calculate regression***********************************************
	reg=matrix(0,length(lon),length(lat));
	pvalue=matrix(0,length(lon),length(lat));

   for(i in seq(1,length(lon),by=1)){
	    for(j in seq(1,length(lat),by=1)){
	
		kitrend=vort[i,j,];
		kjtrend=NAOi;
		kwitrend=kitrend[!is.na(kitrend)];
		kwjtrend=kjtrend[!is.na(kitrend)];
		if(length(kwitrend)>=0.5*length(kitrend)){

		regre=lm(kwitrend~kwjtrend);     #Calculem la regressió lineal.
		reg[i,j]=regre$coefficients[2]; 
		errors=summary.lm(regre);
		pvalue[i,j]=errors$coefficients[2,4];
		}else{
		reg[i,j]=NA;
		pvalue[i,j]=NA;		
		}
	     }
	}

#*************************************NETCDF TO PLOT***************************
 londim=ncdim_def("lon","degrees_east",lon); 
 latdim=ncdim_def("lat","degrees_north",lat); 

 # define variables
	fillvalue=1e32;

dlname="reg";
trend.def=ncvar_def("reg","",list(londim,latdim),fillvalue,dlname,prec="single");

dlname="pvalue";
pvalue.def=ncvar_def("pval","",list(londim,latdim),fillvalue,dlname,prec="single");

# create netCDF file and put arrays
nctitle1=sprintf("%s/newData/reg_%s-%s-NAO_%s_%s_Barnesbox.nc",DIR,term,level,seas,m);
ncfname=nctitle1;
ncout=nc_create(ncfname,list(trend.def,pvalue.def));

# put variables
ncvar_put(ncout,trend.def,reg);
ncvar_put(ncout,pvalue.def,pvalue);

# put additional attributes into dimension and data variables
#ncatt_put(ncout,"p","axis","Z"); #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y");
ncatt_put(ncout,"lon","axis","X");

# close the file, writing data to disk
nc_close(ncout);

}














