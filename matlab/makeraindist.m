function [ppdf1d,pamt1d,ppdf1d2,pamt1d2,bincrates]=makeraindist(pdata1,pdata2,lat,lon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given two epochs of rainfall data, calculate the distributions of rain 
% amount and rain frequency for each, and then area-weight them, following 
% Pendergrass and Hartmann (2014), Two modes of change of the distribution 
% of rain.  
%
% To fit the shift and increase modes to the change in rain distributions 
% calculated here, and then plot the changes as well as the shift and 
% increase modes, see the companion script makeshiftincplots.m 
%
% Angie Pendergrass, 7 August 2015
%
% Demo:
%  load demopdata.mat pdata1 pdata2 lat lon dt
%  [ppdf1,pamt1,ppdf2,pamt2,bincrates]=makeraindist(pdata1,pdata2,lat,lon);
%  save raindistdemodata.mat ppdf1 pamt1 ppdf2 pamt2 bincrates dt
%  makeshiftincplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs 
%
% pdata1: first epoch of daily rain accumulation in mm/d, with shape [lat,lon,days]
%   suggest: 10 years of global data
%   units hint: global mean is 2.6-2.9 mm/d
% pdata2: second epoch of rain data 
% lat,lon: assumed to be regular
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs 
% 
% ppdf1d: area-weighted rain frequency distribution calculated from pdata1 
% pamt1d: area-weighted rain amount distirbution calculated from data1 
% ppdf1d2: area-weighted rain frequency distribution from pdata2 
% pamt1d2: area-weighted rain amount distribution from pdata2 
% bincrates: rain rates at bin center, which are needed for fitting and
%   plotting. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check to see that input data have the right shape. 
sp1=size(pdata1);
sp2=size(pdata2);
if (sp1(1)~=length(lat))&&(sp1(2)~=length(lon))
    disp('pdata1 should be [lat,lon,days]')
end

if (sp2(1)~=length(lat))&&(sp2(2)~=length(lon))
    disp('pdata2 should be [lat,lon,days]')
end

lat=lat(:); lon=lon(:); % make sure lat and lon are column vectors 

% define some constants
L=2.5e6; % W/m2. latent heat of vaporization of water
wm2tommd=1./L*3600*24; % conversion from W/m2 to mm/d

%%% step 1: set up the bins.
pmax=max(max(pdata1(:)),max(pdata2(:))); % find the highest rain accumulation in the dataset

% define initial set of logarithmically-spaced bins
firstupperp=1500; % choose an arbitrary upper bound for initial distribution, in w/m2
minp=1; % arbitrary lower bound for raining threshold, in W/m2
nbins=100;

%%% Some notes: Here, an initial upper threshold and a lower threshold, are 
%%% specified. It might be better to specify the minimum threshold and the
%%% bin spacing, which are around 7% for firstupperp=1500 and minp=1. The 
%%% goals are: 
%%%    - to capture as much of the distribution as possible and 
%%%    - to balance sampling against resolution. 
%%% Capturing the upper end is easy: just extend the bins to include the 
%%% heaviest precipitation event in the dataset (this is done below). The 
%%% lower end is harder: it can go all the way to machine epsilon, and 
%%% there is no obvious reasonable threshold for "rain" over a large 
%%% spatial scale. The value I chose here captures 99.971% of rainfall in a 
%%% CESM run.

binrlog=linspace(log(minp),log(firstupperp),nbins);
dbinlog=diff(binrlog);
binllog=binrlog-dbinlog(1);
binr=exp(binrlog)./L*3600*24;
binl=exp(binllog)./L*3600*24;

dbin=dbinlog(1);
binrlogex=binrlog;
binrend=exp(binrlogex(end));

% extend the bins until the maximum precip anywhere in the dataset falls
% within the bins
while pmax>binr(end)
    binrlogex(end+1)=binrlogex(end)+dbin;
    binrend=exp(binrlogex(end));
    binrlog=binrlogex;
    binr=exp(binrlog)./L*3600*24;
end
binllog=binrlog-dbinlog(1);
binl=exp(binllog)./L*3600*24; % rain rate at bin left edges; this will be used to make distributions
bincrates=[0 (binl+binr)]/2; % rain rate at bin centers; we'll use this for plotting.

% calculate the distributions 
disp('epoch 1...')
[ppdfmap,pamtmap]=makedists(pdata1,binl);

disp('epoch 2...')
[ppdfmap2,pamtmap2]=makedists(pdata2,binl);

% area-weight the distributions 
weight=repmat(cosd(lat),[1 length(lon)]); % cosine-weight. this is correct for rectilinear grids, and only a little off for nearly-rectilinear grids. 
weight=weight./(nansum(nansum(weight))); % normalize the matrix. 
weightp=repmat(permute(weight,[1 2 3]),[1 1 size(ppdfmap,3)]); % replicate the matrix for each bin. 
ppdf1d=squeeze(nansum(nansum(ppdfmap.*weightp,1),2)); % area-weighted rain frequency
pamt1d=squeeze(nansum(nansum(pamtmap.*weightp,1),2)); % area-weighted rain amount 

weightp=repmat(permute(weight,[1 2 3]),[1 1 size(ppdfmap2,3)]); % repeat area-weighting for epoch 2. 
ppdf1d2=squeeze(nansum(nansum(ppdfmap2.*weightp,1),2)); 
pamt1d2=squeeze(nansum(nansum(pamtmap2.*weightp,1),2));

%%% This is a fudge.  Set the precip between the lower p threshold and zero
%%% to zero. 
pamt1d(1)=0;
pamt1d2(1)=0;

    function [thisppdfmap,thispamtmap]=makedists(prdctl,binl)
        % calculate the distributions of rain amount and rain frequency
        % from daily rainfall data and bin left edges. 
        nlat=size(prdctl,1);
        nlon=size(prdctl,2);
        nd=size(prdctl,3);
        
        % use matlab's histc function, which returns histogram bins and 
        % also indices pointing to where the bins came from in the original
        % dataset 
        [n,bin]=histc(prdctl,[0 binl(1:(end-1)) inf],3);  
        ndmat=repmat(nansum(n,3),[1 1 length(binl)+1]); % make a matrix of the total number of days at each gridpoint. 
        thisppdfmap=n./ndmat; % normalize the histogram counts by the total number of days to get frequencies
        
        locmat=repmat(permute(reshape(1:(nlat*nlon),[nlat nlon]),[1 2 3]),[1 1 nd]); % a matrix of indices of the rainfall data 
        locmat=locmat(:); % make the location matrix into a vector
        bin=bin(:); % make the bin-number matrix into a vector
        prdctlvec=double(prdctl(:)); % data type of rain data must be double  
        idx=find(bin>0); % use only non-missing data
        testpamtmap=full(sparse(bin(idx),locmat(idx),prdctlvec(idx))); % make a sparse matrix of the rain falling on days in each bin   
        testpamtmap=(cat(1,testpamtmap,zeros(length(binl)+1-size(testpamtmap,1),nlat*nlon))); % exploit sparse matrix algebra to sum the rain from all the days in each bin 
        thispamtmap=reshape(permute(testpamtmap,[3 2 1]),[nlat nlon length(binl)+1])./ndmat; % reshape back to lat/lon space and normalize by the total number of days at each gridpoint. 
    end
end
