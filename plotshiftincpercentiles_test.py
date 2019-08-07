
# Plot change in percentiles of CMIP5 daily precipitation data and fitted modes
# Uses output from calcshiftincdist_test.pro
# Angie Pendergrass, 7 August 2015
#           update  2 October 2015: improved accuracy of dp/dlnr 


# Translated to python, R Chadwick July 2019
###########################################################################
# Inputs 
#
# dfreq: Change in frequency of each rain-rate percentile
# dfreqshift: Change in frequency of each rain-rate percentile from shift mode
# dfreqinc: Change in frequency of each rain-rate percentile from increase mode
# dfreqshiftinc:Change in frequency of each rain-rate percentile from both modes

import numpy as np
import matplotlib.pyplot as plt

pdfdir = 'your_pdf_path/' # Change this to the path containing input to and output from calcshiftincdist_test.py
plotdir = 'your_plot_path/' # Change this to the path where you want to save your plots

pbinm=1-10**(np.arange(-.1,-4,-.1)) # percentile bins 
xtickpercent=np.concatenate((np.arange(30,100,10),np.arange(91,100), np.arange(99.1,99.9,.1), np.arange(99.91,99.99,.01)))/100 # percentiles for x-axes of percentile plots 
xticks99=np.interp(xtickpercent,pbinm,np.arange(len(pbinm))) # bin numbers of percentiles for plotting

dfreq = np.load(pdfdir+'pr_dfreq_test.npy')
dfreqshift = np.load(pdfdir+'pr_dfreqshift_test.npy')
dfreqinc = np.load(pdfdir+'pr_dfreqinc_test.npy')
dfreqshiftinc = np.load(pdfdir+'pr_dfreqshiftinc_test.npy')
        
### Plot 
plt.figure()
ax=plt.gca()

plt.plot(range(0,len(pbinm)),dfreq,'k')
plt.plot(range(0,len(pbinm)),dfreqshiftinc,'r')
                                 
plt.setp(ax,xticks=xticks99,xticklabels=['','','','','','','90','','', '','','','','','','99' ,'','','','','','','','','99.9','','','','','','','',''])

plt.ylim((-20,20))
plt.xlabel('Percentile')
plt.ylabel('Rain rate change (%/K)')        
filename=plotdir+'prpercentiles_test.pdf'
plt.savefig(filename)
print("wrote "+filename)
plt.close()
        
        




