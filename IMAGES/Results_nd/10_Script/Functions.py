# -*- coding: utf-8 -*-
"""
"""

# ***** References *****
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import numpy as np

 
def Setpath(Path,L):

    """ Determinate the path for result. """
    if L == 1:
        Pathin = Path + '01_Plot_d1/'
        Pathout = Pathin
    elif L == 10:
        Pathin = Path + '02_Plot_d10/'
        Pathout = Pathin
    elif L == 25:
        Pathin = Path + '03_Plot_d25/'
        Pathout = Pathin
    elif L == 50:
        Pathin = Path + '04_Plot_d50/'
        Pathout = Pathin
    elif L == 100:
        Pathin = Path + '05_Plot_d100/'
        Pathout = Pathin

    return Pathin,Pathout
#%%
def fun_inter (Vx,Vy,Vxn):

    """ from a vector Vx and vector Vy, determinate a vector Vyn interpolate Vyn 
    
    """
    Vyn = InterpolatedUnivariateSpline(Vx, Vy)(Vxn)  # para suavizar con spline
    return (Vyn)

#%%
def Recal_amplitude (xmin,xmax,xr,xc,ampr,ampc):

    """ Compute a interpolated amplitude """
    
    s = 0.02
    xn=np.arange(xmin,xmax,s)

    amprn= fun_inter(xr,ampr,xn)             # Interpolated Amplitude for the canyon
    ampcn= fun_inter(xc,ampc,xn)             # Interpolated Amplitude for the  curve canyon
    #
    return (xn,amprn,ampcn)

#%%
def FuncError(Amprn,Ampcn):

    """Funtion to compute errors between two signals
        
        Amprn, Ampcn  : Amplitudes by two apprach
        ErrorProm     : mean error
    
    """

    Error = (Amprn - Ampcn)**2
    ErrorProm = np.sqrt(sum(Error) / len(Amprn))
    return ErrorProm

#%%
def delta(Xr,Xc,Ampr,Ampc,xm):
    
    """Function for computing erros wave for canyons 
        ErrorProm : mean error between two approac
    
    """
    
    #%% Define limits for compute
    #ErrorProm = np.zeros(3)
    
    # ***** Compute spectrum for total surface *****
    xmin1 = -2*xm; xmax1 = 2 * xm
    (Xn,Amprn,Ampcn) = Recal_amplitude(xmin1,xmax1,Xr,Xc,Ampr,Ampc)
    deltaProm = FuncError(Amprn,Ampcn)

#    # ***** Compute spectrum on horizontal line *****
#    xmin1 = -xm; xmax1 = xm
#    (Xn1,Amprn1,Ampcn1) = Recal_amplitude(xmin1,xmax1,Xr,Xc,Ampr,Ampc)
#    ErrorProm[1] = FuncError(Amprn1,Ampcn1)
#
#    # ***** Compute spectrum on slope *****
#    xmin1 = -xm; xmax1 = xm
#    (Xn2,Amprn2,Ampcn2) = Recal_amplitude(xmin1,xmax1,Xr,Xc,Ampr,Ampc)
#    ErrorProm[2] = FuncError(Amprn2,Ampcn2)

    return deltaProm

#%% 
def limredim(XR,YR,XC,YC,AmpR,AmpC):
    
    """Function determinate limits of compute
        X,Y : coordenates curve canyon.
    
    """
    
    #%% Define limits for compute
    rows1=0
    for cont in range(len(XC)): 
        if XC[cont] < 0.0:
            rows1 = rows1 + 1
    
    XC = XC[rows1:]; YC= YC[rows1:];  AmpC = AmpC[rows1:]
    
    #%% define limits for plotting and fitting
    rows1=0
    for cont in range(len(XC)):
        if YC[cont] != 0.0:
            rows1 = rows1 + 1
    xmc = XC[rows1]
    
    # Define mirror for geometry
    XC = np.append(-XC[1:][::-1],XC); YC = np.append(-YC[1:][::-1],YC)
   
    rows1=0
    for cont in range(len(XR)): 
        if XR[cont] < 0.0:
            rows1 = rows1 + 1
            
    XR = XR[rows1:]; YR = YR[rows1:];  AmpR = AmpR[rows1:]

    rows1=0
    for cont in range(len(XR)):
        if YR[cont] != 0.0:
            rows1 = rows1 + 1
    xmr = XR[rows1]

    XR = np.append(-XR[1:][::-1],XR); YR = np.append(-YR[1:][::-1],YR)

    
    return xmr,xmc,XR,YR,XC,YC,AmpR,AmpC
#%% 
def PlotFT(ploter,xmax,XR,XC,AmpR,AmpC):
    
    import Plot as ploter				# 	Python Module by control of time

    """ Function for plotting Functions Transfers  
    
    """

    Nplotx = 3.0; Nploty = 3.4
    Size = [23.0/2.54/Nplotx,16/2.54/Nploty]
    marg =[0.05 * Nploty , 0.95, 0.10 * Nplotx, 0.90]
    Format = ['','$x/a$','$\mid FT \mid$','$%.1f$','$%.1f$','linear','linear','']
         
    # **** Create figure *****
    ylim1 = 0.0; ylim2 = 4.0
    Limx = [-3, 3, 5]; Limy = [ylim1, ylim2, 5]

    figFT, plotFT = plt.subplots(1,1,figsize=(Size[0], Size[1]),dpi=200)
    figFT.subplots_adjust(bottom=marg[0],top=marg[1],left=marg[2],right=marg[3])        
    plotFT.plot(XR/xmax,AmpR, color='blue', linewidth=1.0,linestyle='-')
    plotFT.plot(XC/xmax,AmpC, color='black', linewidth=1.5, linestyle='--')
    ploter.formatplot(figFT, plotFT, Limx, Limy, Format)
        
    #Name = Pathout + 'FT_' + Type + 'SH' + Kind.geo[1:3] + '_L' + str(Kind.LD[NLD]) + '_' + str(int(Frequency[NFig] * 10))
    #plt.savefig(Name + '.pdf', dpi=200); 
    #plt.close(figFT)
        
 #%% funcion
def formatplotgs(plot, row, col, Format):
    """
     Given nr signals each one of nt samples
     plot each signal and its FAS.
     
     cont, i , j, rows, cols: control loop
	 Limx,Limy    :  array by 3 with information about the grid for plotting: [xi,xr,dx]
     Format  : array of 8 formats for title and axis: [title,labelx,labely,fmx,fmy,scalex,scaley,name]

    """
    
    plt.rc('font',family='Times New Roman')
    plot.grid(True, which="both",ls="-")
    font = 'Times New Roman'
    sizetick = 8.0; sizelabel = 9; pad =4; sizetitle = 10
    
    
    ylim1 = 0.0; ylim2 = 4.0
    Limx = [-3, 3, 5]; Limy = [ylim1, ylim2, 5]
    xticks = np.linspace(Limx[0], Limx[1], Limx[2], endpoint=True)    
    plot.set_xticks(xticks)
    plot.set_xlim(Limx[0], Limx[1])
    plot.set_xticklabels([Format[1] % x for x in xticks])
    yticks = np.linspace(Limy[0], Limy[1], Limy[2], endpoint=True)        
    plot.set_yticks(yticks)        
    plot.set_yticklabels([Format[2] % y for y in yticks])
    plot.set_ylim(Limy[0], Limy[1])
         
    plot.tick_params(labelbottom='off')  
    plot.tick_params(labelleft='off') 
    # 
    if col == 0 : 
        plot.tick_params(labelleft ='on') 
        plot.set_ylabel('$\mid FT \mid$',size=sizelabel, weight='normal',labelpad=pad)            
    if row == 0 or row == 1: 
        plot.tick_params(labelbottom='on')
        plot.set_title(Format[0],fontname=font,size=sizetitle,weight='normal')
  
    if row == 1: 
        plot.set_xlabel('$x/a$',size=sizelabel, weight='normal',labelpad=pad +1)  
    elif row == 2: 
        plot.tick_params(labelleft ='on') 
        plot.tick_params(labelbottom='on')
        plot.set_ylabel('$\delta$',size=sizelabel, weight='normal',labelpad=pad)  
        plot.set_xlabel('$\eta$',size=sizelabel, weight='normal',labelpad=pad)

        ylim1 = 0.0; ylim2 = 0.4
        Limx = [0, 10, 6]; Limy = [ylim1, ylim2, 5]
        xticks = np.linspace(Limx[0], Limx[1], Limx[2], endpoint=True)    
        plot.set_xticks(xticks)
        plot.set_xlim(Limx[0], Limx[1])
        plot.set_xticklabels([Format[1] % x for x in xticks])
        yticks = np.linspace(Limy[0], Limy[1], Limy[2], endpoint=True)        
        plot.set_yticks(yticks)        
        plot.set_yticklabels([Format[2] % y for y in yticks])
        plot.set_ylim(Limy[0], Limy[1])
        
    plot.tick_params(labelsize=sizetick)
    return