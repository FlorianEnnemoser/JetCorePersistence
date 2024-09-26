#
#  Calculate the Y time series and the forcing time series M
#   by Albert Ossó, 2024
#**************************************
library('abind');
library('ncdf4');
library('zoo');
library('MASS');
library("TSA");


# Substitue by your own
#DIR="***"; 
#DIR2="***";

#Select SNAO lobe
northbox=1;
box="northbox";

southbox=0;
#box="southbox";

wholebox=0;
#box="wholebox";

# if waves =0 then we calculate the eddy forcing otherwise the linear wave forcing
waves=1;

seas="JJA";   # select season

models=c("ACCESS-ESM1-5","AWI-ESM-1-1-LR","BCC-CSM2-MR","BCC-ESM1","CanESM5","CESM2-FV2","CESM2","CESM2-WACCM-FV2","CESM2-WACCM","CMCC-CM2-HR4","CMCC-CM2-SR5","CMCC-ESM2","EC-Earth3-AerChem","EC-Earth3-CC","EC-Earth3","EC-Earth3-Veg-LR","EC-Earth3-Veg","IITM-ESM","INM-CM4-8","INM-CM5-0","IPSL-CM5A2-INCA","IPSL-CM6A-LR","MIROC6","MPI-ESM-1-2-HAM","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","SAM0-UNICON","TaiESM1","era5");

nmodels=length(models);
count=1;

for(m in models){   # Loop across models

#***************************************CALCULATE FORCING M*********************************
# Open low level NAO vorticity pattern
title=sprintf("%s/reg_VORTICITY-850-NAO_%s_%s_Barnesbox.nc",DIR,seas,m);
u.nc=nc_open(title,write=FALSE);
vort.low=ncvar_get(u.nc,"reg");  
lat_vort=ncvar_get(u.nc,"lat");
lon_vort=ncvar_get(u.nc,"lon");
vort.low=replace(vort.low,is.na(vort.low),0);

#*********Select Area*************
if(wholebox!=0){        
lat.S=45;                        
lat.N=70;			 
lon.W=-60;        
lon.E=0;           
}
  
if(northbox!=0){
lat.S=57.5;         
lat.N=70;              
lon.W=-60;          
lon.E=0;           
}

if(southbox!=0){
lat.S=45;
lat.N=57.5;       
lon.W=-60;        
lon.E=0;
}
                    
# select coordinates to get right lat and lon
lat1=which(lat_vort==lat.S);
lat2=which(lat_vort==lat.N);	
lon1=which(lon_vort==lon.W);
lon2=which(lon_vort==lon.E);	   
vort.low=vort.low[lon1:lon2,lat1:lat2];     
lat=lat_vort[lat1:lat2];
lon=lon_vort[lon1:lon2];

#*****************************OPEN DYN VARIABLE*******************************************************

if(waves!=0){
term="WAVE";
term.path="wave";
var.name="wave";
filein=sprintf("%s/%s/%s_%s.nc",DIR2,term.path,term,m);

}else{
term="EDDY";
term.path="eddy";
var.name="eddy";
# select highpass or not
#filein=sprintf("%s/%s/%s_%s.nc",DIR2,term.path,term,m);
filein=sprintf("%s/%s/%s_%s_HIGHPASS.nc",DIR2,term.path,term,m);
}

#****************CALCULATE ANOMALIES USING CDO****************
fileout=sprintf("%s/vort_%s.nc",DIR2,m);
command=sprintf("cdo seldate,1980-1-1,2014-11-30 %s %s",filein,fileout);
system(command);

filein=fileout;
fileout=sprintf("%s/VORTICITY_250_%s.nc",DIR2,m);
command=sprintf("cdo sellevel,250 %s %s",filein,fileout);
system(command);

filein=fileout;
fileout1=sprintf("%s/VORTICITY_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo selseas,%s %s %s",seas,filein,fileout1);
system(command);

filein=fileout1;
fileout2=sprintf("%s/vort_seasm_%s.nc",DIR2,m);
command=sprintf("cdo yearmean %s %s",filein,fileout2);
system(command);

filein=fileout2;
fileout3=sprintf("%s/vortclim_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo timmean %s %s",filein,fileout3);
system(command);

fileout=sprintf("%s/vortanom_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo -b 64 sub %s %s %s",fileout1,fileout3,fileout);
system(command);
#****************************

#**********************READ FIELD************************
title=fileout;
u.nc=nc_open(title,write=FALSE);
eddy=ncvar_get(u.nc,var.name);  
lat_eddy=ncvar_get(u.nc,"latitude");
lon_eddy=ncvar_get(u.nc,"longitude");
time.data=ncvar_get(u.nc,"time");
time=seq(1,length(time.data),by=1);
durat=length(time);

#*************SELECT COORDINATES***********************
lat1=which(lat_eddy==lat.S);
lat2=which(lat_eddy==lat.N);	
lon1=which(lon_eddy==lon.W);
lon2=which(lon_eddy==lon.E);

eddy=eddy[lon1:lon2,lat1:lat2,];
lat=lat_eddy[lat1:lat2];
lon=lon_eddy[lon1:lon2];


# ARRANGE eddy and vort into matrices
eddyM=matrix(eddy,durat,length(lon)*length(lat),byrow=TRUE);  # Arrange in S mode (Every row is a map at time t1, colus are timeseries)
vortV=as.vector(vort.low);
vortV_noNAN=replaceNAN(vortV);

# Calculate the M time series
Mforcing=eddyM%*%vortV_noNAN;    # We project the eddy forcing onto low level NAO vorticity field
sdM=sd(Mforcing);
Mforcing=Mforcing/sdM;
Mforcing=Mforcing[,1];


# WRITE THE TIMESERIES IN A .txt file.
if(waves!=0){
titlefile=sprintf("%s/Mforcing_Wave_%s_%s_1980-2014_%s_Barnesbox.txt",DIR2,m,seas,box);
save(Mforcing,file=titlefile);

}else{
titlefile=sprintf("%s/Mforcing_%s_%s_1980-2014_%s_HIGHPASS_Barnesbox.txt",DIR2,m,seas,box);
save(Mforcing,file=titlefile);
}

#************************************CALCULATE Y TIMESERIES*************************************
term="VORTICITY";
term.path="vorticity";
var.name="rel_vort";
#****************CALCULATE ANOMALIES USING CDO****************
filein=sprintf("%s/%s/%s_%s.nc",DIR2,term.path,term,m);
fileout=sprintf("%s/trash/vort_%s.nc",DIR2,m);
command=sprintf("cdo seldate,1980-1-1,2014-11-30 %s %s",filein,fileout);
system(command);

filein=fileout;
fileout=sprintf("%s/trash/VORTICITY_250_%s.nc",DIR2,m);
command=sprintf("cdo sellevel,250 %s %s",filein,fileout);
system(command);

filein=fileout;
fileout1=sprintf("%s/trash/VORTICITY_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo selseas,%s %s %s",seas,filein,fileout1);
system(command);

filein=fileout1;
fileout2=sprintf("%s/trash/vort_seasm_%s.nc",DIR2,m);
command=sprintf("cdo yearmean %s %s",filein,fileout2);
system(command);

filein=fileout2;
fileout3=sprintf("%s/trash/vortclim_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo timmean %s %s",filein,fileout3);
system(command);

fileout=sprintf("%s/trash/vortanom_%s_%s.nc",DIR2,seas,m);
command=sprintf("cdo -b 64 sub %s %s %s",fileout1,fileout3,fileout);
system(command);

#******************read data*********************
title=fileout;
u.nc=nc_open(title,write=FALSE);
vort.high=ncvar_get(u.nc,var.name);  
lat_vort=ncvar_get(u.nc,"latitude");
lon_vort=ncvar_get(u.nc,"longitude");
time.data=ncvar_get(u.nc,"time");
time=seq(1,length(time.data),by=1);
durat=length(time);

#*******select coordinates****************************
lat1=which(lat_vort==lat.S);
lat2=which(lat_vort==lat.N);	
lon1=which(lon_vort==lon.W);
lon2=which(lon_vort==lon.E);	  

vort.high=vort.high[lon1:lon2,lat1:lat2,];      
lat=lat_vort[lat1:lat2];
lon=lon_vort[lon1:lon2];

# ARRANGE upper vort into a matrix
vortHighM=matrix(vort.high,durat,length(lon)*length(lat),byrow=TRUE);  # Arrange in S mode 
vortV=as.vector(vort.low);
vortV_noNAN=replaceNAN(vortV);
YNAO=vortHighM%*%vortV_noNAN;    # We project the upper vorticity onto low level NAO vorticity field
sdM=sd(YNAO);
YNAO=YNAO/sdM;
YNAO=YNAO[,1];
command=sprintf("rm -f %s/trash/*.nc",DIR2);    
system(command);

# WRITE THE TIMESERIES IN A .txt file.
titlefile=sprintf("%s/YNAO_%s_%s_1980-2014_%s_Barnesbox.txt",DIR,m,seas,box);
save(YNAO,file=titlefile);

}


