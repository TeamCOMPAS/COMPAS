#!/usr/bin/env python3


"""
Two scripts made for the figures in the introduction
of cosmic integration. Hardcoded and unedited 
"""
import numpy as np
import matplotlib.pyplot as plt

def layoutAxes(ax, nameX='', nameY='', \
               labelSizeMajor = 10, fontsize = 25, second=False, labelpad=None):
    """
    Tiny code to do the layout for axes in matplotlib
    """
    tickLengthMajor = 10
    tickLengthMinor = 5
    tickWidthMajor  = 1.5
    tickWidthMinor  = 1.5
    
    #rc('axes', linewidth=2)
    #label1 always refers to first axis not the twin 
    if not second:
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
    if second:
        for tick in ax.xaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label2.set_fontsize(fontsize)
            #tick.label1.set_fontweight('bold')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.2)
    ax.tick_params(length=tickLengthMajor, width=tickWidthMajor, which='major')
    ax.tick_params(length=tickLengthMinor, width=tickWidthMinor, which='minor')
    ax.set_xlabel(nameX, fontsize=fontsize,labelpad=labelpad)#,fontweight='bold')
    ax.set_ylabel(nameY, fontsize=fontsize,labelpad=labelpad)#, fontweight='bold')    

    return ax


#SFR prescription for illustration
def SFR_Madau(z): #[Msun yr-1 Gpc-3]
    return 0.015* ((1+z)**2.7) / ( 1 + ((1+z)/2.9)**5.6) * 1e9 





def drawConcentricFigure():
    #Having some fun drawing
    minRedshift   = 0
    maxRedshift   = 2.
    nRedshiftBins = 10.
    redshiftEdges = np.linspace(0,maxRedshift,nRedshiftBins+1) #The bin edges in redshift
    redshifts = 0.5*(redshiftEdges[:-1] + redshiftEdges[1:])      #Central value of each redshift bin
    angles = np.linspace(0,2*np.pi, 1e4)

    fig, axes = plt.subplots(1,1, figsize=(15,3.75))

    for nrBin, zBin in enumerate(redshiftEdges[:-1]):
        #Draw circles of at every redshift edge
        if nrBin == 0:
            label='shell Edges'
        else:
            label = None
        x = redshiftEdges[nrBin] * np.cos(angles)
        y = redshiftEdges[nrBin] * np.sin(angles)
        axes.plot(x,y,c='k', lw=2, linestyle=':', label=label)
        
        #for one shell draw a line thick enough
        #so it looks like the shell filled
        if nrBin == 5:
            x = redshifts[nrBin] * np.cos(angles)
            y = redshifts[nrBin] * np.sin(angles)
            axes.plot(x,y,c='b', lw=35, alpha=0.2, label='volume Shell')
        
    #redshift point at center/halfway point shell which we use for luminosity distance
    axes.scatter(redshifts, np.zeros(len(redshifts)), label='redshift at center shell')

    #draw arrow from detector to redshift for luminosity Distance
    axes.arrow(0, 0.05, redshifts[5]-0.05, 0.0, head_width=0.04, fc='k')
    axes.annotate('luminosity distance to shell', xy=(0.1, 0.1), xytext=(0.1, 0.1), fontsize=15)
    #Extra text
    axes.arrow(1.9, 0.3, 0., -0.22, head_width=0.04, fc='k')
    axes.annotate('calculate mergers at \n each of these points', xy=(1.9, 0.3), xytext=(1.9, 0.3), fontsize=15)
    #Central point to show location detector
    axes.scatter(0,0, c='r', label='detector', s=100)

    axes.get_yaxis().set_visible(False)
    axes.legend(loc=2, prop={'size':18})
    axes.set_ylim(-0.5, 0.5)
    nameX = r'$\rm redshift \ z$'
    nameY = r'$ $'
    axes= layoutAxes(axes, nameX=nameX, nameY=nameY)

    plt.show()



def drawingConcept(x, y):
    fig,axes = plt.subplots(1,1, figsize=(30,9))
    redshifts = np.linspace(0,3,100)
    axes.plot(redshifts, SFR_Madau(redshifts)/float(1e8), c='k', linestyle=':', label='SFR Madau et al.', lw=2.)
    axes.set_ylabel('SFR*1e8 [Msun yr-1 Gpc-3]')
    axes.set_xlabel('redshift z')
    #Extra text
    axes.arrow(0.0+x, 0.3, 0.5, 0.0, head_width=0.05, fc='k')
    axes.arrow(0.1+x, 0.6, -0.07, -0.23, head_width=0.02, fc=None)
    axes.arrow(0.63+x, 0.35+y, -0.05, 0.1, head_width=0.02, fc=None)
    axes.annotate('delay time', xy=(0.1, 0.31), xytext=(0.1+x, 0.31), fontsize=25)
    axes.annotate('binary merges here', xy=(0.1, 0.6), xytext=(0.1+x, 0.6), fontsize=25)
    axes.annotate('binary born here with this SFR', xy=(0.63, 0.35), xytext=(0.63+x, 0.35+y), fontsize=25)
    nameX = r'$\rm redshift \ z$'
    nameY = r'$\rm SFR*1e8 [Msun yr-1 Gpc-3]$'
    axes= layoutAxes(axes, nameX=nameX, nameY=nameY)
    axes.legend(loc=2, prop={'size':25})
    plt.show()
