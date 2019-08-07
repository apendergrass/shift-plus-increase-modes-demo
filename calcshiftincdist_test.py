
###########################################################################
# Given two distributions of rain amount and rain frequency and a change in 
# temperature, 
#   - fit the shift and increase modes to the change, and then 
# following Pendergrass and Hartmann (2014), Two modes of change of the 
# distribution of rain.  
#
# To calculate rain amount and frequency distributions from rainfall data,
# see the companion script makraindist.m 
#
# Angie Pendergrass, 7 August 2015
#           update  2 October 2015: improved accuracy of dp/dlnr 


# Translated to python, R Chadwick July 2019
#
###########################################################################
# Inputs 
#
# ppdf1: area-weighted rain frequency distribution calculated from pdata1 
# pamt1: area-weighted rain amount distirbution calculated from pdata1 
# ppdf2: area-weighted rain frequency distribution from pdata2 
# pamt2: area-weighted rain amount distribution from pdata2 
# dt: change in temperature. Should be a scaler, eg change in global mean
#   surface-air temperature. Set to 1 if you don't want to normalize by
#   temperature.
# bincrates: rain rates at bin center, which are needed for fitting and
#   plotting. 
# 
###########################################################################
# Outputs 
# 
# shift: shift mode (#/K)
# inc: increase mode (#/K)
# err: error of the fit (absolute value of difference between change in 
#      rain amount and shift-plus-increase, #) 
# pamtshift: area-weighted change in rain amount distribution from shift mode
# pamtinc: area-weighted change in rain amount distribution from increase mode 
# pamtshiftinc: area-weighted change in rain amount distribution from both modes 
# ppdfshift; area-weighted change in rain frequency distribution from shift mode
# ppdfinc: area-weighted change in rain frequency distribution from increase mode 
# ppdfshiftinc: area-weighted change in rain frequency distribution from both modes 
# dfreq: Change in frequency of each rain-rate percentile
# dfreqshift: Change in frequency of each rain-rate percentile from shift mode
# dfreqinc: Change in frequency of each rain-rate percentile from increase mode
# dfreqshiftinc:Change in frequency of each rain-rate percentile from both modes   
#
###########################################################################

import numpy as np
import numpy.linalg as lin
import iris # xarray or similar could also be used to import the data

# calculating functions
def findshiftinc(pamt1,pamt21k,bincrates):
    # solve for shift and increase modes to optimally fit the change in rain amount distribution
    #, return shifted and increased rain amount

    dpamt=pamt21k-pamt1 
    db=(bincrates[2]-bincrates[1])/bincrates[1]        
    
    # calculate dp/dln(r)
    nb=len(pamt1);
    dpdlnrc= np.concatenate((np.array([0]),pamt1[2:nb]-pamt1[0:(nb-2)],np.array([0])))/(2*db) # centered difference (2nd order accurate; now used near endpoints)
    prm2=pamt1[:(nb-4)]
    prm1=pamt1[1:(nb-3)]
    prp1=pamt1[3:(nb-1)]
    prp2=pamt1[4:]
    dpdlnr= np.concatenate((np.array([0]),np.array([0]),8*(prp1-prm1)-(prp2-prm2),np.array([0]),np.array([0])))/(12*db) # 4th order accurate 
    dpdlnr[1]=dpdlnrc[1] # fill in with second order accurate near endpoints
    dpdlnr[nb-2]=dpdlnrc[nb-2] 
    
    # calculate the terms in the matrices that will need to be solved
    # (equation 10 in Pendergrass and Hartmann 2014, Two Modes of Change of the Distribution of Rain) 
    spr2=np.nansum(pamt1**2)
    sdpdlnr2=np.nansum(dpdlnr**2)
    spdpdlnr=np.nansum(pamt1*dpdlnr)
    spdpm=np.nansum(pamt1*dpamt)
    sdpdlnrdpm=np.nansum(dpdlnr*dpamt)
    
    # set up the matrices
    B=np.array([[spdpm],[-sdpdlnrdpm]])
    A=np.array([[spr2,-spdpdlnr],[-spdpdlnr,sdpdlnr2]])
    x = lin.solve(A,B)  # solve for shift and increase modes 
    
    inc=100*x[0] # increase mode (%)
    shift=100*x[1] # shift mode (%)
            
    # shifted and increased rain amount distributions 
    pamtshift=-x[1]*dpdlnr+pamt1
    pamtinc=x[0]*pamt1+pamt1
    pamtshiftinc=x[0]*pamt1-x[1]*dpdlnr+pamt1

    # error of the fit (#)
    err=100*np.nansum(np.abs(pamtshiftinc[2:]-pamt21k[2:]))/np.nansum(np.abs(dpamt[2:]))
    
    return (shift,inc,pamtshift,pamtinc,pamtshiftinc,err)



def pdffrompamt(pamt,ppdf,newpamt,bincrates):
    # calculate rain frequency distributions from shifted and increased rain amount
    # This is an indirect estimate of frequency, so an adjustment is applied. 
    
    newpamt[np.isnan(newpamt)]=0 # get rid of NaNs 

    # make sure there is no rain in the zero rain-rate bin. 
    pamt[0]=0; 
    newpamt[0]=0 
    
    tr1=pamt/bincrates # calculate the "time raining in each bin" in the original rain amount distribution. [mm / (mm/d)] Note that this is a fuzzy concept. 
    totaltime=np.sum(tr1[1:])/(1-ppdf[0]) # calculate the equivalent of the total amount of time. 
    
    tr=newpamt/bincrates # calculate the "time raining in each bin" in the new rain amount distribution
    newppdf=tr/totaltime # normalize this by the total amount of time. 

    fadj=ppdf/(tr1/totaltime) # ratio between the initial rain frequency distribution and the frequency of normalized "time raining in each bin." 
                              #  Will be used as an adjustment to the new rain frequency distribution. 
    fadj[np.isinf(fadj)]=1  
    fadj[np.isnan(fadj)]=1
    newppdf=newppdf*fadj # adjust the new rain frequency distribution. 

    newppdf[0]=1-sum(newppdf[1:]) # the important part: the new dry day frequency 
    
    return newppdf
    
    

def dratefunpercentile(ppdf1,ppdf2,bincrates):
    # calculate the change in rain rate as a function of percentile of distribution from rain frequency distributions
 
    prrates1=ratefunpercentile(ppdf1,bincrates)
    prrates2=ratefunpercentile(ppdf2,bincrates)
    dfreq=100*(prrates2-prrates1)/prrates1
    
    return dfreq


    
def ratefunpercentile(ppdf,bincrates):
    # calculate the rain rate as a function of percentile of distribution from a rain frequency distribution

    pbinm=1-10**(np.arange(-.1,-4,-.1)) # percentile bins 
    
    cdf=np.cumsum(ppdf)   
          
    # find the monotonically increasing part of the cdf for interpolation
    startind=1 # initially, start at bin 2; bin 1 is dry frequency and should be much larger than bin 2
    lgz=np.where(np.diff(cdf)<=0)[0] # initially, end where the cdf stops going up. there is no more rain in the higher bins.  
    if not lgz.any(): # in this case, there is rain even in the last bin
        lgz=len(cdf); 
    else:   # otherwise, we need to find the smooth part by force.  
        if (lgz<20).any(): # in this case, something wonky is happening at light rain rates.  
            startind=np.max(lgz[np.where(lgz<20)[0]]) # try again, skipping the light rain rates. 
            
        lgz=lgz[np.where(lgz>20)[0]][0] # try again to find the smooth part, making sure it's not at the beginning. 
        if not lgz.any():
            lgz=len(cdf); # it may still be the case that all the bins have rain, despite something weird happening at light rain rates. 
    
    mrrates=np.interp(pbinm,cdf[startind:lgz],np.arange(startind,lgz)) # interpolate percentiles onto bin indices from the monotonically increasing part of the distribution
    prrates=np.interp(mrrates,np.arange(0,len(bincrates)),bincrates) # interpolate bin indices onto rain rates from the monotonically increasing part of the distribution
        
    return prrates


pdfdir = 'your_pdf_path/' # Change this to the path containing test_spi.nc


ppdf1 = iris.load_cube(pdfdir+'test_spi.nc','ppdf1').data.data
ppdf2 = iris.load_cube(pdfdir+'test_spi.nc','ppdf2').data.data       
pamt1 = iris.load_cube(pdfdir+'test_spi.nc','pamt1').data.data
pamt2 = iris.load_cube(pdfdir+'test_spi.nc','pamt2').data.data

## Load change in global mean tas
## Set this to 1K for the test version
dt = 1.0

## Load bin centre values (same for all models)
bincrates = iris.load_cube(pdfdir+'test_spi.nc','ppdf1').coord('bincrates').points

# calculate rain distribution changes for 1 K warming
pamt21k=pamt1+(pamt2-pamt1)/dt
ppdf21k=ppdf1+(ppdf2-ppdf1)/dt

# solve for shift and increase modes, and get shifted and increased
# rain amount distributions 
modes = findshiftinc(pamt1,pamt21k,bincrates)
shift = modes[0]
inc = modes[1]
pamtshift = modes[2]
pamtinc = modes[3]
pamtshiftinc = modes[4]
err = modes[5]

np.save(pdfdir+'pr_modes_test.npy',modes)

print('Shift '+shift)
print('Increase '+inc)
print('Error '+err)

# calculate rain frequency distributions from shifted and increased rain
# amount 
ppdfshift=pdffrompamt(pamt1,ppdf1,pamtshift,bincrates)
ppdfinc=pdffrompamt(pamt1,ppdf1,pamtinc,bincrates)
ppdfshiftinc=pdffrompamt(pamt1,ppdf1,pamtshiftinc,bincrates)

np.save(pdfdir+'pr_pdfshift_test.npy',ppdfshift)
np.save(pdfdir+'pr_pdfinc_test.npy',ppdfinc)
np.save(pdfdir+'pr_pdfshiftinc_test.npy',ppdfshiftinc)

# calculate the change in rain rate as a function of percentile of distribution from rain
# frequency distributions
dfreq=dratefunpercentile(ppdf1,ppdf21k,bincrates)
dfreqshift=dratefunpercentile(ppdf1,ppdfshift,bincrates)
dfreqinc=dratefunpercentile(ppdf1,ppdfinc,bincrates)
dfreqshiftinc=dratefunpercentile(ppdf1,ppdfshiftinc,bincrates)

np.save(pdfdir+'pr_dfreq_test.npy',dfreq)
np.save(pdfdir+'pr_dfreqshift_test.npy',dfreqshift)
np.save(pdfdir+'pr_dfreqinc_test.npy',dfreqinc)
np.save(pdfdir+'pr_dfreqshiftinc_test.npy',dfreqshiftinc)


