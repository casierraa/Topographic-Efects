
# -*- coding: utf-8 -*-

from IPython import get_ipython
get_ipython().magic('reset -sf') 

"""
Created on Sat Nov 15 16:25:42 2014
Modified on Feb 10. Plots for each eta

@author: Cesar Sierra - Mario Saenz 

Title: FTtotal
Porpuse: Progam to show Tranfer function by whole points on space and a frequency
The error is calcuculate

**** Files for load *****

load('FT_Amp.txt')        from Forder Canonics Forms
load('FT_Imag.txt')       from Forder Canonics Forms
load('FT_Imag.txt')       from Forder Canonics Forms
"""

# %% paths
#'/Volumes/MacintoshHD/
pathsystem = 'C:/Users/caugu/'
pathsystem = '/Volumes/MacintoshHD/'
pathdoc = pathsystem + 'Dropbox/01_Eafit/01_Doctorate/31_Papers/'
pathbins = pathsystem + 'Dropbox/01_Eafit/01_Doctorate/31_Papers/10_Codes/06_Scripts/'
pathproy = '03_BandMethod/'
pathinter = '01_Canonics/'

#from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import Functions as fn
from time import time				# Python Module by control of time
import os
os.sys.path.append(pathbins)
import Plot as ploter

start_time = time()								#	Control of time

plt.close("all")
print("start Proccess")  

#%% define some data
case = ['15_degree','20_degree','25_degree','30_degree','35_degree','40_degree']
nd = [1,10,25,50]
id = ['01_','02_','03_','04_'] #['01_','02_','03_','04_']
curva = ['_curvo_','_recto_']
s = 0.01                            # s _  : Step for interpolating
d = 0.10                            # d _  : distance betteen straiht and curve approach in [km]                
form = ['/01_Canyon/', '/02_Trapecio/','/03_Hill/','/04_Colina/']
wave = ['SH', 'SV']

file = 'FTX_A.txt'

color=[]; color.append('black'); color.append('blue');color.append('green'); color.append('magenta'); color.append('red')
label=[]; label.append('Lw = 1'); label.append('Lw = 10'); label.append('Lw = 25'); label.append('Lw = 50');
color=[]; color.append('#68affc'); color.append('#1f9383');color.append('#4ad9e1'); color.append('#1f3e9e'); color.append('dimgray')
angle = []
eta_ld = np.zeros([4,6])    # allocate [ld,degree]

for idform in range(0,1): #['/01_Canyon/', '/02_Trapecio/','/03_Hill/','/04_Colina/']
    for idcase in range(6): #['15_degree','20_degree','25_degree','30_degree','35_degree','40_degree']
        pathcanyon1 = '0' + str(idcase + 1) + '_Canonic' + case[idcase][0:2] + form[idform]    
        folderld = pathdoc + pathproy + pathinter + 'Results_nd'  + form[idform] 
        # ***** Load Files *****

        for idwave in range(1):#     ['SH', 'SV']  
                                    
            Nplotx = 1.1; Nploty = 1.1
            Size = [16.0/2.54/Nplotx, 23/2.54/Nploty]
            figgs = plt.figure(figsize=(Size[0], Size[1]))
            marg =[0.06 * Nploty , 0.95, 0.08 * Nplotx, 0.94, 0.30]
            plotgs = (gridspec.GridSpec(nrows=3, ncols=2, height_ratios=[1, 1, 1],
                     bottom=marg[0],top=marg[1],left=marg[2],right=marg[3],hspace=marg[4]))
            
            #%% print difference curve and straigth approaches for each Ld
            Nplotx = 3.0; Nploty = 3.4
            Size = [23.0/2.54/Nplotx,16/2.54/Nploty]
            marg =[0.08 * Nploty , 0.95, 0.06 * Nplotx, 0.94]
            Limx = [0, 5, 6]; Limy = [0, 0.3, 4]
            Format = ['','$\eta$','$\delta$','$%.1f$','$%.1f$','linear','linear','']
            figdelta, plotdelta = plt.subplots(1,1,figsize=(Size[0], Size[1]),dpi=200)
            figdelta.subplots_adjust(bottom=marg[0],top=marg[1],left=marg[2],right=marg[3])                          
            
            deltalim = 0.10
            plotdelta.axhline(y=deltalim, color=color[4], linewidth=2.0, linestyle='--')
            angle = np.append(angle,float(case[idcase][0:2]))
                   
            for ld in range(4):#  [1,10,25,50]

                #%% Read the amplitud of transfer functions SH and SV horizontal
                """
                Frequency 				     : array frequency
                XR(SH,SV,P); XC(SH,SV,P)      : points on the straight canoyn  and the circular one
                Amp(XRSH)                     : displacement value
                FTAmp : Transfer functiona    : SV,SH are the kind of waves
                """  
                geome = wave[idwave] + curva[0] + 'd' + str(nd[ld]) + '/'                    
                pathcanyon = pathcanyon1 + id[ld] + 'd' + str(nd[ld]) + '/' + geome; print(pathcanyon)   
                Folder = pathdoc + pathproy + pathinter + pathcanyon + 'Output/'
            
                FTAmpC = np.loadtxt(Folder + file)
                geome = wave[idwave] + curva[1] + 'd' + str(nd[ld]) + '/'                    
                pathcanyon = pathcanyon1 + id[ld] + 'd' + str(nd[ld]) + '/' + geome; print(pathcanyon) 
                Folder = pathdoc + pathproy + pathinter + pathcanyon + 'Output/'
 
                FTAmpR = np.loadtxt(Folder + file); 

                Frequency = FTAmpR[0, 2:]; eta = Frequency * d * 100.0    # The frequency is scaleeta = 
                XR = FTAmpR[1: -1, 0];   YR = FTAmpR[1: -1, 1]; AmpR = FTAmpR[1: -1,2: ]
                XC = FTAmpC[1: -1, 0];   YC = FTAmpC[1: -1, 1]; AmpC = FTAmpC[1: -1,2: ]
                              
                #%% Compute differences between curve and straigth approaches for each frequency
                Totfreq = len(Frequency)
                               
                xmr,xmc,XR,YR,XC,YC,AmpR,AmpC = fn.limredim(XR,YR,XC,YC,AmpR,AmpC)
                                              
                #for Nf in range(0,len(Frequency)):  
                eval = 0
                delta = np.zeros([1,Totfreq])
                for Nf in range(Totfreq):
                    AmpC1 = AmpC[:,Nf]; AmpR1 = AmpR[:,Nf]                    
                    AmpC1 = np.append(AmpC1[1:][::-1],AmpC1)
                    AmpR1 = np.append(AmpR1[1:][::-1],AmpR1)
                    
                    if idcase == 0: Amp = 1.4
                    elif idcase == 1: Amp = 1.1
                    else: Amp = 1.0 
                    
                    delta[:,Nf] = (fn.delta(XR,XC,AmpR1,AmpC1,xmc) * Amp)
                     
                    if ld == 0: row = 0; col = 0
                    if ld == 1: row = 1; col = 0
                    if ld == 2: row = 0; col = 1
                    if ld == 3: row = 1; col = 1
                    if delta[0,Nf] >= 0.10 and eval == 0: 
                        eta_ld[ld,idcase] = eta[Nf]; eval = 1
                        subplot = figgs.add_subplot(plotgs[row, col])
                        subplot.plot(XR/xmr,AmpR1, color='black', linewidth=1.5, linestyle='--', label = label[ld])
                        subplot.plot(XC/xmr,AmpC1, color=color[ld], linewidth=1.0,linestyle='-', label = label[ld])
                        Format1 = [label[ld],'$%.1f$','$%.1f$']
                        fn.formatplotgs(subplot, row, col, Format1)                       
                    
                    elif delta[0,Totfreq-1] < 0.10 and eval == 0:
                        eta_ld[ld,idcase] = 5.0
                     
                plotdelta.plot(eta[0:Totfreq], delta[0, :],color[ld], linewidth=1.0, linestyle='-', label = label[ld])                
                
                subplot = figgs.add_subplot(plotgs[2:, :])
                subplot.plot(eta[0:Totfreq], delta[0, :],color[ld], linewidth=1.0, linestyle='-', label = label[ld]) 
                subplot.axhline(y=deltalim, color=color[4], linewidth=2.0, linestyle='--')
                fn.formatplotgs(subplot, 2, 2, Format1)
                subplot.legend(loc='upper left',ncol = 2, prop={'size':7},framealpha=0.2)
            # ***** save plot by .pdf *****                      
            Format[7] = folderld + 'delta_nd' + '_' + wave[idwave] + '_' + case[idcase][0:2]          
            plotdelta.legend(loc='upper left',ncol=2, prop={'size':6},framealpha=0.2)
            ploter.formatplot(figdelta, plotdelta, Limx, Limy, Format); #plt.close(figdelta)
            
            figgs.savefig(Format[7] + '_FT.pdf', dpi=100); #plt.close(figgs)
            
    #%% Create figures with 
    Nplotx = 2.0; Nploty = 3.4
    Size = [16.5/2.54/Nplotx,23/2.54/Nploty]
    marg =[0.05 * Nploty , 0.95, 0.08 * Nplotx, 0.92]
    Limx = [15.0, 40.0, 6]; Limy = [1, 6, 6]
    Format = ['','$\\alpha$','$\eta_{max}$','$%.1f$','$%.1f$','linear','linear','']
    figld, plotld = plt.subplots(1,1,figsize=(Size[0], Size[1]),dpi=200)
    figld.subplots_adjust(bottom=marg[0],top=marg[1],left=marg[2],right=marg[3])                          
    
    for ld in range(4):             
        plotld.plot(angle, eta_ld[ld], color[3],marker='.', linewidth=0.0, label = label[ld])
    
    #plotld.legend(loc='upper left',ncol = 2, prop={'size':6})
    plotld.axhspan(2.0, 4.0, facecolor='0.5', alpha=0.5)
    Format[7] = folderld + 'delta_nd' + '_' + wave[idwave]

    ploter.formatplot(figld, plotld, Limx, Limy, Format); #plt.close(figld)


#%% end process
elapsed_time = time() - start_time

print("End Proccess")
print("Elapsed time: %.10f seconds." % elapsed_time)
