function [shift,inc,err]=makeshiftincplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given two distributions of rain amount and rain frequency and a change in 
% temperature, 
%   - fit the shift and increase modes to the change, and then 
%   - make plots showing 
%       - the change in rain amount and frequency distributions, 
%       - the change in rain rate as a function of percentile, and 
%       - the shift and increase modes, 
% following Pendergrass and Hartmann (2014), Two modes of change of the 
% distribution of rain.  
%
% To calculate rain amount and frequency distributions from rainfall data,
% see the companion script makraindist.m 
%
% Angie Pendergrass, 7 August 2015
%           update  2 October 2015: improved accuracy of dp/dlnr 
%
% Demo:
%  load raindistdemodata.mat ppdf1 pamt1 ppdf2 pamt2 bincrates dt
%  [shift,inc,err]=makeshiftincplots(ppdf1,pamt1,ppdf2,pamt2,dt,bincrates);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs 
%
% ppdf1d: area-weighted rain frequency distribution calculated from pdata1 
% pamt1d: area-weighted rain amount distirbution calculated from data1 
% ppdf1d2: area-weighted rain frequency distribution from pdata2 
% pamt1d2: area-weighted rain amount distribution from pdata2 
% dt: change in temperature. Should be a scaler, eg change in global mean
%   surface-air temperature. Set to 1 if you don't want to normlaize by
%   temperature.
% bincrates: rain rates at bin center, which are needed for fitting and
%   plotting. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs 
% 
% figure 1: 9-panel plot showing changes in rain amount, frequency and rain
% rate as a function of percentile, as well as the shift and increase modes
% shift: shift mode (%/K)
% inc: increase mode (%/K)
% err: error of the fit (absolute value of difference between change in 
%      rain amount and shift-plus-increase, %) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations 

% calculate rain distribution changes for 1 K warming
pamt21k=pamt1+(pamt2-pamt1)./dt;
ppdf21k=ppdf1+(ppdf2-ppdf1)./dt;

% solve for shift and increase modes, and get shifted and increased
% rain amount distributions 
[shift,inc,pamtshift,pamtinc,pamtshiftinc,err]=findshiftinc(pamt1,pamt21k,bincrates);

% calculate rain frequency distributions from shifted and increased rain
% amount 
ppdfshift=pdffrompamt(pamt1,ppdf1,pamtshift,bincrates);
ppdfinc=pdffrompamt(pamt1,ppdf1,pamtinc,bincrates);
ppdfshiftinc=pdffrompamt(pamt1,ppdf1,pamtshiftinc,bincrates);

% calculate the change in rain rate as a function of percentile of distribution from rain
% frequency distributions
[dfreq]=dratefunpercentile(ppdf1,ppdf21k,bincrates);
[dfreqshift]=dratefunpercentile(ppdf1,ppdfshift,bincrates);
[dfreqinc]=dratefunpercentile(ppdf1,ppdfinc,bincrates);
[dfreqshiftinc]=dratefunpercentile(ppdf1,ppdfshiftinc,bincrates);


    % calculating functions
    function [shift,inc,pamtshift,pamtinc,pamtshiftinc,err]=findshiftinc(pamt1,pamt21k,bincrates)
        % solve for shift and increase modes to optimally fit the change in rain amount distribution
        %, return shifted and increased rain amount

        % setup
        bincrates=bincrates(:);
        pamt1=pamt1(:);
        pamt21k=pamt21k(:);
        dpamt=pamt21k-pamt1; 
        db=(bincrates(3)-bincrates(2))./bincrates(2);        
        
        % calculate dp/dln(r)
        % dpdlnr=[0; diff(pamt1(:))./db]; % old method: forward euler (1st order accurate)
        nb=length(pamt1);
        dpdlnrc=[0;pamt1(3:nb)-pamt1(1:(nb-2));0;]./(2*db); % centered difference (2nd order accurate; now used near endpoints)
        prm2=pamt1(1:(nb-4));
        prm1=pamt1(2:(nb-3));
        prp1=pamt1(4:(nb-1));
        prp2=pamt1(5:nb);
        dpdlnr=[0;0;8*(prp1-prm1)-(prp2-prm2);0;0;]./(12*db); % 4th order accurate 
        dpdlnr([2 nb-1])=dpdlnrc([2 nb-1]); % fill in with second order accurate near endpoints

        % calculate the terms in the matrices that will need to be solved
        % (equation 10 in Pendergrass and Hartmann 2014, Two Modes of Change of the Distribution of Rain) 
        spr2=nansum(pamt1.^2);
        sdpdlnr2=nansum(dpdlnr.^2);
        spdpdlnr=nansum(pamt1.*dpdlnr);
        spdpm=nansum(pamt1.*dpamt);
        sdpdlnrdpm=nansum(dpdlnr.*dpamt);
        
        % set up the matrices
        B=[spdpm; -sdpdlnrdpm];
        A=[spr2 -spdpdlnr; -spdpdlnr sdpdlnr2];
        x=A\B; % solve for shift and increase modes 
        
        inc=100*x(1); % increase mode (%) 
        shift=100*x(2); % shift mode (%)
                
        % shifted and increased rain amount distributions 
        pamtshift=-x(2).*(dpdlnr)+pamt1;
        pamtinc=x(1).*pamt1+pamt1;
        pamtshiftinc=x(1).*pamt1-x(2).*(dpdlnr)+pamt1;

        % error of the fit (%)
        err=100*nansum(abs(pamtshiftinc(3:end)-pamt21k(3:end)))./nansum(abs(dpamt(3:end)));
    end

    function [newppdf]=pdffrompamt(pamt,ppdf,newpamt,bincrates)
        % calculate rain frequency distributions from shifted and increased rain amount 
        bincrates=bincrates(:);
        newpamt(isnan(newpamt))=0; % get rid of NaNs 

        % make sure there is no rain in the zero rain-rate bin. 
        pamt(1)=0; 
        newpamt(1)=0; 
        
        tr1=pamt./bincrates; % calculate the "time raining in each bin" in the original rain amount distribution. [mm / (mm/d)] Note that this is a fuzzy concept. 
        totaltime=sum(tr1(2:end))./(1-ppdf(1)); % calculate the equivalent of the total amount of time. 
        
        tr=newpamt./bincrates; % calculate the "time raining in each bin" in the new rain amount distribution
        newppdf=tr./totaltime; % normalize this by the total amount of time. 

        fadj=ppdf./(tr1./totaltime); % ratio between the initial rain frequency distribution and the frequency of normalized "time raining in each bin." Will be used as an adjustment to the new rain frequency distribution. 
        fadj(isinf(fadj))=1;  
        fadj(isnan(fadj))=1;
        newppdf=newppdf.*fadj; % adjust the new rain frequency distribution. 

        newppdf(1)=1-sum(newppdf(2:end)); % the important part: the new dry day frequency 
    end

    function [dfreq]=dratefunpercentile(ppdf1,ppdf2,bincrates)
        % calculate the change in rain rate as a function of percentile of distribution from rain frequency distributions

        bincrates=bincrates(:);
        
        % calculate the rain rate as a function of percentile of distribution from a rain frequency distribution
        [prrates1]=ratefunpercentile(ppdf1,bincrates);
        [prrates2]=ratefunpercentile(ppdf2,bincrates);
        dfreq=100*(prrates2-prrates1)./prrates1;
        
        function [prrates]=ratefunpercentile(ppdf,bincrates)
            % calculate the rain rate as a function of percentile of distribution from a rain frequency distribution

            pbinm=(1-10.^(-.1:-.1:-4)); % percentile bins 
            cdf=cumsum(ppdf);         
            % find the monotonically increasing part of the cdf for interpolation
            startind=2; % initially, start at bin 2; bin 1 is dry frequency and should be much larger than bin 2
            lgz=find(diff(cdf)<=0); % initially, end where the cdf stops going up. there is no more rain in the higher bins.  
            if isempty(lgz); % in this case, there is rain even in the last bin
                lgz=length(cdf); 
            else   % otherwise, we need to find the smooth part by force.  
                if ~isempty(find(lgz<20)); % in this case, something wonky is happening at light rain rates.  
                    startind=1+max(lgz(find(lgz<20))); % try again, skipping the light rain rates. 
                end
                lgz=lgz(find(lgz>20,1)); % try again to find the smooth part, making sure it's not at the beginning. 
                if isempty(lgz);
                    lgz=length(cdf); % it may still be the case that all the bins have rain, despite something weird happening at light rain rates. 
                end
            end
            mrrates=interp1(cdf(startind:lgz),startind:lgz,pbinm);% interpolate percentiles onto bin indices from the monotonically increasing part of the distribution
            prrates=interp1(1:length(bincrates),bincrates,mrrates);% interpolate bin indices onto rain rates from the monotonically increasing part of the distribution
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preliminary items

% colors 
shiftcolor=[27 158 119]/256; 
inccolor=[117 112 179]/256; 
spicolor=[231 41 138]/256; 

% x-axis label data
xtickrates=[0 0.1:.1:.9 1:1:9 10:10:90 100:100:1000]; % rain rates in mm/d for x axis 
xticks=interp1(bincrates,1:length(bincrates),xtickrates); % bin numbers associated with nice number rain rate 
xtickratesp={0 0.1 [] [] [] [] [] [] [] [] 1 [] [] [] [] [] [] [] [] 10 [] [] [] [] [] [] [] [] 100 [] [] [] [] [] [] [] [] 1000}; % rain rates with white space in between 
pbinm=(1-10.^(-.1:-.1:-4)); % percentile bins for calculations 
xtickpercent=[30:10:90 91:99 99.1:.1:99.9 99.91:.01:99.99]/100; % percentiles for x-axes of percentile plots 
xticks99=interp1(pbinm,1:length(pbinm),xtickpercent); % bin numbers of percentiles for plotting

% bin width
db=(bincrates(3)-bincrates(2))./bincrates(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);clf
set(gcf,'units','centimeters','paperpositionmode','auto');
set(gcf,'position',[1 1 16.9 16.9]); % 8.4 16.9 23.7 are the lengths

% rain amount 
subplot(3,3,1) % shift 
dpamtplot((pamt21k-pamt1)*1./db,'k',(pamtshift-pamt1)*1./db,shiftcolor,xticks,xtickratesp);
l=legend('Model','Shift','location','northwest');set(l,'box','off')
ylabel('\DeltaRain amount (mm/d/K)');% If you wanted to be precise, the y-axis label is in mm/d/K/\Delta bin
title(['Shift = ' num2str(shift,2)  '%'])

subplot(3,3,2) % increase
dpamtplot((pamt21k-pamt1)*1./db,'k',(pamtinc-pamt1)*1./db,inccolor,xticks,xtickratesp);
l=legend('Model','Increase','location','northwest');set(l,'box','off')
title(['Increase = ' num2str(inc,2)  '%'])

subplot(3,3,3) % shift+increase
dpamtplot((pamt21k-pamt1)*1./db,'k',(pamtshiftinc-pamt1)*1./db,spicolor,xticks,xtickratesp);
l=legend('Model','Shift+inc','location','northwest');set(l,'box','off')
title(['Shift + increase'])

% rain frequency distribution
subplot(3,3,4)
dppdfplot((ppdf21k-ppdf1),'k',(ppdfshift-ppdf1),shiftcolor,xticks,xtickratesp);
ylabel('\DeltaRain amount (mm/d/K)');% If you wanted to be precise, the y-axis label is in mm/d/K/\Delta bin

subplot(3,3,5)
dppdfplot((ppdf21k-ppdf1),'k',(ppdfinc-ppdf1),inccolor,xticks,xtickratesp);

subplot(3,3,6)
dppdfplot((ppdf21k-ppdf1),'k',(ppdfshiftinc-ppdf1),spicolor,xticks,xtickratesp);

% rain rate as a function of percentile
subplot(3,3,7)
drpercentileplot(dfreq,'k',dfreqshift,shiftcolor,pbinm,xticks99)
ylabel('\DeltaRain rate (%/K)')

subplot(3,3,8)
drpercentileplot(dfreq,'k',dfreqinc,inccolor,pbinm,xticks99)

subplot(3,3,9)
drpercentileplot(dfreq,'k',dfreqshiftinc,spicolor,pbinm,xticks99)


    % plotting functions 
    function []=dpamtplot(dpamt1,color1,dpamt2,color2,xticks,xtickratesp)
        % plot the change in rain amount distribution.  two curves in two
        % different colors.  
        nb=length(dpamt1);
        p=plot(1:nb,dpamt1,'-k',1:nb,dpamt2,'-r');
        set(p(1),'color',color1)
        set(p(2),'color',color2)
        set(p,'linewidth',1.5)
        
        l=line([0 130],[0 0]);set(l,'color',[1 1 1]*.5) % add a zero line 
        ylim([-.4 1]*.103) % choose y-axes that worked for CESM 
        xlim([4 130]) % exclude the zero bin and also the last few, where little should be happening. 
        set(gca,'xtick',xticks,'xticklabel',xtickratesp) % make x-axes show rain rate in mm/d
        
        %  ylabel('\DeltaRain amount (mm/d/K)');% If you wanted to be precise, the y-axis label is in mm/d/K/\Delta bin
        xlabel('Rain rate (mm/d)');
    end

    function []=dppdfplot(dppdf1,color1,dppdf2,color2,xticks,xtickratesp)
        % plot the change in rain frequency distribution.  two curves in two
        % different colors.  
        
        nb=length(dppdf1);

        p=plot(1:nb,100*dppdf1,'-k',1:nb,100*dppdf2,'-r');
        set(p(1),'color',color1)
        set(p(2),'color',color2)
        set(p,'linewidth',1.5)
        
        l=line([0 130],[0 0]);set(l,'color',[1 1 1]*.5) % add a zero line 
        xlim([4 130]); % exclude the zero bin (which will throw off y-axis) and also the last few, where little should be happening. 
        set(gca,'xtick',xticks,'xticklabel',xtickratesp) % make x-axes show rain rate in mm/d
        
        ylim([-.04 .04]) % choose y-axes that worked for CESM 
        
        % add the changes in dry-day frequency at the top left of the plot
        a=axis;
        t1=text(a(1),a(4)-1*.01,[' ' num2str(100*dppdf1(1),2)],'horizontalalignment','left','verticalalignment','top','color',color1,'fontsize',10);
        t4=text(a(1),a(4),[' ' num2str(100*dppdf2(1),2)],'horizontalalignment','left','verticalalignment','top','color',color2,'fontsize',10);
        
        %  ylabel('\DeltaFrequency (%/K)'); % If you wanted to be precise, the y-axis label is in mm/d/K/\Delta bin
        xlabel('Rain rate (mm/d)');
                
    end

    function []=drpercentileplot(dfreq1,color1,dfreq2,color2,pbinm,xticks99)
        % plot the change in rain rate as a function of percentile.  two curves in two
        % different colors.  
        
        p=plot(1:length(pbinm),dfreq1,'-k',1:length(pbinm),dfreq2,'-m');
        set(p(1),'color',color1)
        set(p(2),'color',color2)
        set(p,'linewidth',1.5)
        
        ylim([-10 20]) % choose y-axes that worked for CESM 
        
        % add minor ticks. without forcing them to have labels.  Code from
        % matlab central. 
        grid on
        xg = xticks99;
        xx = reshape([xg;xg;NaN(1,length(xg))],1,length(xg)*3);
        a=axis;
        mth=(a(4)-a(3))/100;
        yg = a(3)+[0 mth];
        yy = repmat([yg NaN],1,length(xg));
        yg2 = a(4)+[0 -mth];
        yy2 = repmat([yg2 NaN],1,length(xg));
        hold on
        h_minorgrid = plot(xx,yy,'k',xx,yy2,'k');
        hold off

        set(gca,'xtick',10:10:length(pbinm),'xticklabel',{90 99 99.9 99.99}) % make x-axes show percentile
        
        l=line([2 length(pbinm)],[0 0]);set(l,'color',[1 1 1]*.5) % add zero line 
        xlim([1 length(pbinm)])
        
        % ylabel('\DeltaRain rate (%/K)')
        xlabel('Percentile')
    end

end
