
# Plot histograms of change in CMIP5 daily precipitation data and fitted shift and increase modes
# Uses output from calcshiftincdist_test.pro, and a file test_spi.nc
# Angie Pendergrass, 7 August 2015
#           update  2 October 2015: improved accuracy of dp/dlnr 


# Translated to python, R Chadwick July 2019
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
# bincrates: rain rates at bin center
# pamtshift: area-weighted change in rain amount distribution from shift mode
# pamtinc: area-weighted change in rain amount distribution from increase mode 
# pamtshiftinc: area-weighted change in rain amount distribution from both modes 
# ppdfshift; area-weighted change in rain frequency distribution from shift mode
# ppdfinc: area-weighted change in rain frequency distribution from increase mode 
# ppdfshiftinc: area-weighted change in rain frequency distribution from both modes    


import numpy as np
import matplotlib.pyplot as plt
import iris # xarray or similar could also be used to import the data

pdfdir = 'your_pdf_path/' # Change this to the path containing input to and output from calcshiftincdist_test.py
plotdir = 'your_plot_path/' # Change this to the path where you want to save your plots

ppdf1 = iris.load_cube(pdfdir+'test_spi.nc','ppdf1').data.data
ppdf2 = iris.load_cube(pdfdir+'test_spi.nc','ppdf2').data.data       
pamt1 = iris.load_cube(pdfdir+'test_spi.nc','pamt1').data.data
pamt2 = iris.load_cube(pdfdir+'test_spi.nc','pamt2').data.data

## Load change in global mean tas
## Set this to 1K for the test version
dt = 1.0

## Load bin centre values (same for all models)
bincrates = iris.load_cube(pdfdir+'test_spi.nc','ppdf1').coord('bincrates').points
        
         
#### Normalize to get changes per degree warming. 
dppdfnorm = (ppdf2-ppdf1)/dt
dpamtnorm = (pamt2-pamt1)/dt

# Calculate number of dry days
dry1 = ppdf1[0]*100
dry2 = ppdf2[0]*100
ddry = (dry2-dry1)/dt

modes = np.load(pdfdir+'pr_modes_test.npy')
pamtshift = modes[2]
pamtinc = modes[3]
pamtshiftinc = modes[4]
  
ppdfshift = np.load(pdfdir+'pr_pdfshift_test.npy')
ppdfinc = np.load(pdfdir+'pr_pdfinc_test.npy')
ppdfshiftinc = np.load(pdfdir+'pr_pdfshiftinc_test.npy')          

dpamtshift = pamtshift - pamt1          
dpamtinc = pamtinc - pamt1 
dpamtshiftinc = pamtshiftinc - pamt1          
dppdfshift = ppdfshift - ppdf1          
dppdfinc = ppdfinc - ppdf1 
dppdfshiftinc = ppdfshiftinc - ppdf1
        
# rain rates in mm/d for x axis ticks and labeling 
otn=np.linspace(1,9,9)
xtickrates=np.append(0,otn*.1)
xtickrates=np.append(xtickrates,otn)
xtickrates=np.append(xtickrates,otn*10)
xtickrates=np.append(xtickrates,otn*100)
xticks=np.interp(xtickrates,bincrates,range(0,len(bincrates))); #% bin numbers associated with nice number rain rate
xticks,indices=np.unique(xticks,return_index=True)

### Bin width - needed to normalize the rain amount distribution
db=(bincrates[2]-bincrates[1])/bincrates[1];
binwidth = db*100
print('binwidth (%)',binwidth)

### Plot 
plt.figure(figsize=(4,6))
plt.clf()
ax=plt.subplot(211)
plt.plot(range(0,len(pamt1)),dpamtnorm/db, 'k')
plt.plot(range(0,len(pamt1)),dpamtshiftinc/db, 'r')
plt.plot(range(0,len(pamt1)),dpamtinc/db, 'b')
plt.plot(range(0,len(pamt1)),dpamtshift/db, 'g')
                        
plt.plot((0,len(pamt1)),(0,0),'0.5')
plt.setp(ax,xticks=xticks,xticklabels=[''])
    #plt.xlabel('Rain rate (mm/d)')
plt.ylim((-.08,.12))
plt.gca().set_xlim([1, len(bincrates)])
plt.title('Rain amount change (mm/day/K)')

ax=plt.subplot(212)
plt.plot(range(0,len(ppdf1)),dppdfnorm*100/db, 'k')  
plt.plot(range(0,len(ppdf1)),dppdfshiftinc*100/db, 'r') 
plt.plot(range(0,len(pamt1)),dppdfinc*100/db, 'b')
plt.plot(range(0,len(pamt1)),dppdfshift*100/db, 'g')  
                         
plt.plot((0,len(ppdf1)),(0,0),'0.5')
    ### Annotate with the dry day frequency
t=plt.text(4,.4, "{:.1f}".format(ddry)+'%')
plt.setp(t,va='top',ha='left')
plt.setp(ax,xticks=xticks,xticklabels=['0','0.1','','','','','','','','','1','','','','','','','','','10','','','','','','','','','100','','','','','','','','','1000'])
plt.ylim((-0.6,0.6))
plt.gca().set_xlim([1, len(bincrates)])
plt.xlabel('Rain rate (mm/day)')
plt.title('Rain frequency change (%/K)')
#plt.show()
filename=plotdir+'raindistshiftinc_test.pdf'
plt.savefig(filename)
print("wrote "+filename)
plt.close()
        
        




