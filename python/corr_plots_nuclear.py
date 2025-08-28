#load libraries 
import time,os
from matplotlib.offsetbox import AnchoredText
import sys,pandas as pd, matplotlib , matplotlib.pyplot as plt, matplotlib.lines , numpy as np,cupy as cp, math, pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit
import matplotlib.ticker as mticker
import matplotlib.ticker as ticker

viridis = cm.get_cmap('viridis', 12)
inferno = cm.get_cmap('inferno', 12)
autumn = cm.get_cmap('autumn', 12)
bwr = cm.get_cmap('bwr', 20)

import matplotlib.colors as colors

def offset(a, amount=np.pi/2):
    return a+2*np.pi*(a<-np.pi+amount)-2*np.pi*(a>=np.pi+amount)

def summary_plot(df,df_mixed,df_triggers=None,bins=10,projyrange=(0,2.2), bins1d=14, text="", plotFourier=3,ylimrat=(0,0.02),
                include_1d=True,zlims=None,zlimrat=None,
            surfrat = True, bottomRowOnly=False, topRowOnly=False,logzrat=False):
    deta_range=projyrange
    
    #elevation and azimuth for surface plots
    elev = 70
    azim = -112.5#-135
    
    
    if include_1d:
        fig = plt.figure(figsize=(20,8))
        ax_r1 = fig.add_subplot(gs[0,5])
        ax_r2 = fig.add_subplot(gs[1,5],sharex = ax_r1)
        ax_r1.tick_params(axis='x',direction='inout',bottom=True)
        ax_r2.tick_params(axis='x',direction='inout',top=True)
        gs = fig.add_gridspec(ncols=6, nrows=2, width_ratios=[1,1,1,1,0.5,1])
    else:
        if not bottomRowOnly and not topRowOnly:
            fig = plt.figure(figsize=(16,8))
            gs = fig.add_gridspec(ncols=5, nrows=2, width_ratios=[1,1,1,1,0.15])
        elif topRowOnly:
            fig = plt.figure(figsize=(16,4))
            gs = fig.add_gridspec(ncols=5, nrows=1, width_ratios=[1,1,1,1,0.15])
        elif bottomRowOnly:
            fig = plt.figure(figsize=(12,4))
            gs = fig.add_gridspec(ncols=4, nrows=1, width_ratios=[1,1,1,0.15])
    
    axt_d = None
    axb_c = None
    
    yc_d = None
    dyc_d = None
    
    hists_2d = []
    
    for i in range(7):
        dphi_range = (-np.pi/2,3*np.pi/2)
        

        if df_triggers is None:
            denom = len(df[i])*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins
        else:
            denom = len(df_triggers[i])*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins

        #hist1,xedges, yedges = np.histogram2d([2,2,2,2],[0,0,0,0], bins=bins,range=[deta_range, dphi_range])
        hist1, xedges, yedges = np.histogram2d(df[i].diff_rap_cm, offset(df[i].diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
        hist1 = np.divide(hist1, denom)



        #print(xedges)
        hist2, xedges, yedges = np.histogram2d(df_mixed[i].diff_rap_cm, offset(df_mixed[i].diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
        for j in range(bins):
            if xedges[j+1]>0:
                break
        for j in range(bins):
            if yedges[j+1]>0:
                break;

        hist2 = np.divide(hist2,len(df_mixed[i])*(xedges[1]-xedges[0])*(yedges[1]-yedges[0]))

        hist3 = np.divide(hist1,hist2)
        hists_2d.append(hist3)

        xpos, ypos = np.meshgrid(np.add(xedges[:-1],xedges[1:])/2, np.add(yedges[:-1],yedges[1:])/2)
        zpos = 0
        #help(hist1.transpose)

        for j in range(len(xedges)-1):
            if(xedges[j]>projyrange[0]):
                break
        #print("i=",i)
        def plot2d(ax, hist,cm = viridis,vmin=None,vmax=None, surf = True,log=False):
            if surf:
                if not log:
                    return ax.plot_surface(xpos, ypos, hist.transpose(), cmap=cm,edgecolor='k',vmin=vmin,vmax=vmax)
                else :
                    def log_tick_formatter(val, pos=None):
                        return f"$10^{{{int(val)}}}$"
                        #return "%.1f"%10**val
                    
                    ret = ax.plot_surface(xpos, ypos, np.log10(hist.transpose()), cmap=cm,edgecolor='k',vmin=np.log10(vmin),vmax=np.log10(vmax))
                    ax.set_zlim(np.log10(vmin),np.log10(vmax))
                    ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
                    ax.zaxis.set_major_locator(ticker.MultipleLocator(1))
                    return ret;
            
            else:
                return ax.pcolor(xpos,ypos,hist.transpose(), cmap=cm,vmin=vmin,vmax=vmax, norm=colors.LogNorm(vmin=vmin, vmax=vmax) if log else None)
        if i<=3:
            if bottomRowOnly:
                continue
            ax1 = fig.add_subplot(gs[0,i], projection='3d', sharez = axt_d)
            ax1.set_title('D C Fe Pb'.split()[i],fontsize='xx-large')
            if i == 0:
                axt_d = ax1

            
            ax1.view_init(azim=azim,elev = elev)
            ax1.set_xlabel("$\\Delta y$",fontsize='large')
            ax1.set_ylabel("$\\Delta\\phi$ [rad]",fontsize='large')
            ax1.set_zlabel("             $C(\\Delta\\phi,\\Delta y)$\n",fontsize='large',linespacing=3.2)
            ax1.zaxis.set_rotate_label(False)
            if not zlims is None:
                ax1.set_zlim(*zlims)
            
            p2d = plot2d(ax1,hist3, vmin=zlims[0] if not zlims is None else None, vmax=zlims[1] if not zlims is None else None)
            if i == 3:
                fig.colorbar(p2d, orientation='vertical',cax=fig.add_subplot(gs[0,-1]))
        if i >= 4 :
            if topRowOnly:
                break
            #use surface plot for the ratio?
            if surfrat:
                ax = fig.add_subplot(gs[-1,i-3-int(bottomRowOnly)], projection='3d', sharez = axb_c)
                ax.view_init(azim=azim,elev = elev)
                ax.zaxis.set_rotate_label(False)
                if not zlimrat is None:
                    
                    ax.set_zlim(*zlimrat)
                    #if logzrat:
                    #    ax.set_zscale('log')
            else :
                ax = fig.add_subplot(gs[-1,i-3-int(bottomRowOnly)])
            if i == 4:
                axb_c = ax
        
            #surf = ax3.plot_surface(xpos, ypos, hist3.transpose(), cmap=viridis,edgecolor='k')
            p2d = plot2d(ax,np.divide(hists_2d[i-3],hist3),cm=bwr,vmin=zlimrat[0] if not zlimrat is None else None,vmax=zlimrat[1] if not zlimrat is None else None, surf=surfrat,log=logzrat)
            if i == 6:
                def log_tick_formatter(val, pos=None):
                    return f"$10^{{{int(val)}}}$"
                cax = fig.add_subplot(gs[-1,-1])
                
                cbar = fig.colorbar(p2d, orientation='vertical',cax=cax)
                if logzrat:
                    cbar.ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
                    cbar.ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.set_xlabel("$\\Delta y$",fontsize='large')
            ax.set_ylabel("$\\Delta\\phi$ [rad]",fontsize='large')
            a = "D C Fe Pb".split()[i-3]
            ax.set_title(f"$\\frac{{C_{{\\mathrm{{{a}}}}}(\\Delta\\phi,\\Delta y)}}{{C_D(\\Delta\\phi,\\Delta y)}}$", fontsize='x-large')
            #ax.set_title("$C(\\Delta\\phi,\\Delta y)$",rotation=0)
            ax.set_xlim(*deta_range)
            ax.set_ylim(*dphi_range)
            
            
        
        if not include_1d:
            continue
        
        if df_triggers is None:
            denom = len(df[i])*2*np.pi/bins1d
        else:
            denom = len(df_triggers[i])*2*np.pi/bins1d
        y, x = np.histogram(offset(df[i].query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange).diff_phi_cm), bins=bins1d, range=(-np.pi/2,3*np.pi/2))
        x = np.add(x[1:],x[:-1])/2
        dy = np.sqrt(y)
        y=np.divide(y,denom)
        dy=np.divide(dy,denom)
        
        mixed_highy = df_mixed[i].query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
        denom = len(mixed_highy)/bins1d
        ym, _ = np.histogram(offset(mixed_highy.diff_phi_cm), bins=bins1d, range=(-np.pi/2,3*np.pi/2))
        dym = np.sqrt(ym)
        ym=np.divide(ym,denom)
        dym=np.divide(dym,denom)
        
        dyc = y/ym*np.hypot(dy/y, dym/ym)
        yc= y/ym
        
        if i == 0:
            dyc_d = dyc
            yc_d = yc
        color = 'k tab:blue tab:orange tab:green'.split()[i]
        label = 'D C Fe Pb'.split()[i]
        fillstyle='full'
        ax_r1.errorbar(x,yc,dyc,marker='o',color=color,linestyle='',label=label,fillstyle=fillstyle)
    
        if plotFourier:
            n = plotFourier
            b = [2*sum(yc*np.cos(k*x))/len(x) for k in range(0,n+1)]
            b[0]/=2
            s = 0
            xs = np.linspace(-np.pi/2,3*np.pi/2,50)
            for k in range(0,n+1):
                s = b[k]*np.cos(k*xs)+s
            ax_r1.plot(xs,s,linestyle='-',color=color)
         
        if i != 0:
            dyc_rat = yc/yc_d*np.hypot(dyc/yc, dyc_d/yc_d)
            yc_rat = yc/yc_d
            ax_r2.errorbar(x,yc_rat,dyc_rat,marker='o',color=color,linestyle='',label=label,fillstyle=fillstyle)

            if plotFourier:
                n = plotFourier
                b = [2*sum(yc_rat*np.cos(k*x))/len(x) for k in range(0,n+1)]
                b[0]/=2
                s = 0
                xs = np.linspace(-np.pi/2,3*np.pi/2,50)
                for k in range(0,n+1):
                    s = b[k]*np.cos(k*xs)+s
                ax_r2.plot(xs,s,linestyle='-',color=color)
    if include_1d:
        ax_r1.legend(loc='upper left',framealpha=0)
        ax_r1.set_xlabel("$\\Delta\\phi$ [rad]",fontsize='large')
        ax_r1.set_ylim(*ylimrat)
        ax_r2.set_xlabel("$\\Delta\\phi$ [rad]",fontsize='large')
        ax_r1.set_ylabel("$C(\\Delta\\phi)$",fontsize='large')
        ax_r2.set_ylabel("$C_A(\\Delta\\phi)/C_D(\\Delta\\phi)$",fontsize='large')
        ax_r1.set_xlim(-np.pi/2,3*np.pi/2)
        ax_r2.axhline(1,color='0.7',linestyle='--')
    fig.subplots_adjust(wspace=0)
    
       
    if text and not topRowOnly and not bottomRowOnly:
        fig.text(0.15,0.25,text,bbox=dict(facecolor='white', alpha=1),fontsize='large')


