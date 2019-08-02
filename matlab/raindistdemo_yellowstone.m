%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This demo is meant to walk you through calculating the distribution of 
%%% rain from an input data set of daily rainfall accumulation, as described
%%% in: 
%%% Pendergrass, A.G. and D.L. Hartmann, 2014: Two modes of change of the 
%%%   distribution of rain. Journal of Climate, 27, 8357-8371. 
%%%   doi:10.1175/JCLI-D-14-00182.1.  
%%% and the shift and increase modes of response of the rainfall distribution
%%% to warming, occuring across ENSO events or global warming simulations. 
%%% The response to warming is described in: 
%%% Pendergrass, A.G. and D.L. Hartmann, 2014: Changes in the distribution 
%%%   of rain frequency and intensity in response to global warming. 
%%%   Journal of Climate, 27, 8372-8383. doi:10.1175/JCLI-D-14-00183.1. 

%%% Please cite one or both of these papers if you use or alter these scripts. 

%%% THIS script will get you going quickly with CESM large
%%% ensemble data if you happen to have access to the NCAR supercomputer, 
%%% yellowstone https://www2.cisl.ucar.edu/resources/yellowstone 

%%% But if you dont, thats ok, any gridded precipiation dataset will do. 
%%% See the companion script for more information. 

%%% 14 January 2016, Angeline Pendergrass, NCAR, Boulder CO. apgrass@uw.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  input data
%%% first 10 years

%%%% CESM LE 
%/glade/p/cesm0005/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/PRECT
%b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h1.PRECT.20060101-20801231.nc  b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h1.PRECT.20810101-21001231.nc
%%%

%lat
%lon

% 10 years of daily p data in mm/d. [lat,lon,days].  global mean is 2.6-2.9 mm/d.
%pdata1 % first epoch
%pdata2 % second (assumed warmer) epoch

% global mean surface (air, if you want) temperature change
%dt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filedir='/glade/p/cesm0005/CESM-CAM5-BGC-LE/atm/proc/tseries/daily/PRECT/';
file1='b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h1.PRECT.20060101-20801231.nc';
lat=ncread([filedir file1],'lat');
lon=ncread([filedir file1],'lon');

pr=ncread([filedir file1],'PRECT');
pr=pr(:,:,1:3650);

L=2.5e6; % w/m2. latent heat of vaporization of water
wm2tommd=1./L*3600*24; % conversion from w/m2 to mm/d

pdata1=permute(pr,[2 1 3])*1000*L*wm2tommd;

file2='b.e11.BRCP85C5CNBDRD.f09_g16.001.cam.h1.PRECT.20810101-21001231.nc';
pr=ncread([filedir file2],'PRECT');
pr=pr(:,:,3651:7300);
pdata2=permute(pr,[2 1 3])*1000*L*wm2tommd;
clear pr


tfiledir='/glade/p/work/apgrass/code/lensptimeseries/';
tfile='TREFHT.001.rcp85.nc';
tas=ncread([tfiledir tfile],'TREFHT');
years=1920:2100;
years(87:96)
years(172:181)

dt=mean(tas(172:181)-tas(87:96));

%%% dt=4.0470;  %%% in case the files in my directories arent readable
%%% anymore

%save demopdata.mat pdata1 pdata2 lat lon dt

[ppdf1,pamt1,ppdf2,pamt2,bincrates]=makeraindist(pdata1,pdata2,lat,lon);


%save raindistdemodata.mat ppdf1 pamt1 ppdf2 pamt2 bincrates dt

makeshiftincplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates);

