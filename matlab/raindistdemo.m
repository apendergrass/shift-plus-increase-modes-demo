
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

%%% A companion demo script will get you going quickly with CESM large
%%% ensemble data if you happen to have access to the NCAR supercomputer, 
%%% yellowstone https://www2.cisl.ucar.edu/resources/yellowstone 

%%% But if you dont, thats ok, any gridded precipiation dataset will do. 
%%% For example, you can use daily rainfall data from a CMIP5 simulation 
%%% (GFDL-ESM2G is shown below, for the 1pctCO2 scenario)
%%% which you might be able to download, for example from PCMDI 
%%% (https://pcmdi.llnl.gov/projects/cmip5/)

%%% Or you can use data from GPCP 1dd 
%%% https://climatedataguide.ucar.edu/climate-data/gpcp-daily-global-precipitation-climatology-project
%%% or TRMM 3B42
%%% https://climatedataguide.ucar.edu/climate-data/trmm-tropical-rainfall-measuring-mission
%%% though looking at changes in those is different, and you'll have to get
%%% corresponding temperature data somewhere, if you're into that. 

%%% 14 January 2016, Angeline Pendergrass, NCAR, Boulder CO. apgrass@uw.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% SAMPLE INPUT FILE NAMES. You'll want to update these. If you actually
%%% use the GFDL simulation, you'll have to concatenate data from 2
%%% files... I'll assume you can figure that out. 

%%% first 10 years from GFDL ESM2G
filedir='/directory/with/data/';
file1='pr_day_GFDL-ESM2G_1pctCO2_r1i1p1_00010101-00051231.nc';
% file1='pr_day_GFDL-ESM2G_1pctCO2_r1i1p1_00060101-00101231.nc';

% 10 years immediately prior to co2 doubling 
file2='pr_day_GFDL-ESM2G_1pctCO2_r1i1p1_00610101-00651231.nc';
% file2='pr_day_GFDL-ESM2G_1pctCO2_r1i1p1_00660101-00701231.nc';

% To get the warming, you'll have to do a little data processing yourself. 
% This assumes an input file with global-annual mean surface temperature
% over the entire simulation 
tfiledir='/where/your/temperature/data/is/';
tfile='tas_global-annual-mean.nc';

%%% Otherwise, just get your daily precip data, in mm/d, and in an array
%%% with dimension order [lat,lon,days]. 
%%% As a sanity check on units, global, annual-mean rainfally should be around 2.6-2.9 mm/d.

%pdata1 % first epoch
%pdata2 % second (assumed warmer) epoch

% global mean surface (or near-surface air) temperature change
%dt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat=ncread([filedir file1],'lat');
lon=ncread([filedir file1],'lon');

pr=ncread([filedir file1],'pr'); 
%%% potentially deal with concatenation first. 
pr=pr(:,:,1:3650);

L=2.5e6; % w/m2. latent heat of vaporization of water
wm2tommd=1./L*3600*24; % conversion from w/m2 to mm/d

pdata1=permute(pr,[2 1 3])*1000*L*wm2tommd;

pr=ncread([filedir file2],'pr');
pr=pr(:,:,3651:7300);
pdata2=permute(pr,[2 1 3])*1000*L*wm2tommd;
clear pr

tas=ncread([tfiledir tfile],'tas');
years=1:70;
years(1:10)
years(61:70)

dt=mean(tas(61:70)-tas(1:10));

%save demopdata.mat pdata1 pdata2 lat lon dt

[ppdf1,pamt1,ppdf2,pamt2,bincrates]=makeraindist(pdata1,pdata2,lat,lon);


%save raindistdemodata.mat ppdf1 pamt1 ppdf2 pamt2 bincrates dt

makeshiftincplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates);

