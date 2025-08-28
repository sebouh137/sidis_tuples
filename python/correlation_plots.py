#load libraries 
import time,os
from matplotlib.offsetbox import AnchoredText
import sys,pandas as pd, matplotlib , matplotlib.pyplot as plt, matplotlib.lines , numpy as np,cupy as cp, math, pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit

viridis = cm.get_cmap('viridis', 12)
inferno = cm.get_cmap('inferno', 12)
autumn = cm.get_cmap('autumn', 12)

def offset(a, amount=np.pi/2):
    return a+2*np.pi*(a<-np.pi+amount)-2*np.pi*(a>=np.pi+amount)
def corr1d(df,df_mixed,bins=20,style='normal',area=None,
                          fig=None,projyrange=(1.5,2.5),bins1d=50,
           label=None,color=None,showFourier=True,Ntriggers=None):
    dphi_range = (-np.pi/2,3*np.pi/2)
    
    bins = bins1d
    #now for 1d projections:
    #ax4 = fig.add_subplot(234)
    
    if Ntriggers == None:
        denom = len(df)*2*np.pi/bins
    else :
        denom = Ntriggers/bins
    y, x = np.histogram(offset(df.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange).diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
    x = np.add(x[1:],x[:-1])/2
    dy = np.sqrt(y)
    y=np.divide(y,denom)
    dy=np.divide(dy,denom)
    
    
    def plot1d(x,y,dy,ax,n=5):
        if style=='normal':
            ax.errorbar(x,y,dy,marker='o',color='k'  if color==None else color,linestyle='', label = label)
            ax.axhline(0,color='0.7',linestyle=':')
            if showFourier:
                #a = [2*sum(y*np.sin(i*x))/len(x) for i in range(0,n+1)]
                #remove sin term
                #a = [0 for i in range(0,n+1)]
                b = [2*sum(y*np.cos(i*x))/len(x) for i in range(0,n+1)]
                b[0]/=2
                s = 0
                for i in range(0,n+1):
                    s = b[i]*np.cos(i*x)+s
                ax.plot(x,s,linestyle='-',color='k',label='all' if label == None else None)
                for i in range(1,n+1):
                    ax.plot(x,b[i]*np.cos(i*x)+b[0]*(i!=0),linestyle='--', label="n=%s"%i if label == None else None)
            ax.set_xlim(-np.pi/2,np.pi*3/2)
        elif style=='fill':
            ax.bar(x,y,x[1]-x[0],color='#ff7777',alpha=0.5)
    
    if df_mixed is None:
        plot1d(x,y,dy,plt.gca())
        return
    
    mixed_highy = df_mixed.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
    
    denom = len(mixed_highy)/bins1d
    ym, _ = np.histogram(offset(mixed_highy.diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
    ym = ym+(ym==0)*(max(ym)*.000001)
    dym = np.sqrt(ym)
    ym=np.divide(ym,denom)
    dym=np.divide(dym,denom)
    
    dyc = y/ym*np.hypot(dy/y, dym/ym)
    yc= y/ym
    if(area != None):
        scale = area/(sum(yc)*(x[1]-x[0]))
        yc = yc*scale
        dyc =dyc*scale
    
    plot1d(x,yc,dyc,plt.gca())
    plt.gca().set_xlabel("$\\Delta\\phi$ [rad]")
    return x,yc, dyc

def dphi_deta_plot_3_proj(df,df_mixed,deta_range=(-1.5,2.5),bins=20,df_triggers=None,
                          fig=None,projyrange=(1.5,2.5),bins1d=50,minbincontentindenom=0.003, 
                          normalizeMixedAt00=False,text="",zlims=None,ylims=None,fourierOrder=3):
    dphi_range = (-np.pi/2,3*np.pi/2)
    if fig == None:
        fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(231, projection='3d')
    
    elev = 70
    azim = -112.5#-135
    ax1.view_init(azim=azim,elev = elev)
    
    if df_triggers is None:
        denom = len(df)*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins
    else:
        denom = len(df_triggers)*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins
    
    #hist1,xedges, yedges = np.histogram2d([2,2,2,2],[0,0,0,0], bins=bins,range=[deta_range, dphi_range])
    hist1, xedges, yedges = np.histogram2d(df.diff_rap_cm, offset(df.diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
    hist1 = np.divide(hist1, denom)
    
    
    
    print(xedges)
    hist2, xedges, yedges = np.histogram2d(df_mixed.diff_rap_cm, offset(df_mixed.diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
    for i in range(bins):
        if xedges[i+1]>0:
            break
    for j in range(bins):
        if yedges[j+1]>0:
            break;
            
    radius = .3
    if normalizeMixedAt00:
        M00 = len(df_mixed.query('sqrt(diff_rap_cm**2+diff_phi_cm**2)<%s'%radius))/(np.pi*radius**2)
        #print("M(0,0) = ", M00)
        
        hist2 = np.divide(hist2,M00*(xedges[1]-xedges[0])*(yedges[1]-yedges[0]))
    else :
        hist2 = np.divide(hist2,len(df_mixed)*(xedges[1]-xedges[0])*(yedges[1]-yedges[0]))
    
    
    
    
    hist3 = np.divide(hist1,np.add(hist2,0.0001))
    
    
    
    xpos, ypos = np.meshgrid(np.add(xedges[:-1],xedges[1:])/2, np.add(yedges[:-1],yedges[1:])/2)
    zpos = 0
    #help(hist1.transpose)
    
    for i in range(len(xedges)-1):
        if(xedges[i]>projyrange[0]):
            break
    print("i=",i)
    def plot2d(ax, hist):
        if projyrange != deta_range:
            ax.plot_surface(xpos.transpose()[:i+1].transpose(), ypos.transpose()[:i+1].transpose(), hist[:i+1].transpose(), cmap=viridis,edgecolor='k')
            ax.plot_surface(xpos.transpose()[i:].transpose(), ypos.transpose()[i:].transpose(), hist[i:].transpose(), cmap=autumn,edgecolor='k')
        else :
            ax.plot_surface(xpos, ypos, hist.transpose(), cmap=viridis,edgecolor='k')
    plot2d(ax1,hist1)
    ax1.set_xlabel("$\\Delta y$")
    ax1.set_ylabel("$\\Delta\\phi$ [rad]")
    ax1.zaxis.set_rotate_label(False)
    ax1.set_title("$S(\\Delta\\phi,\\Delta y)$",rotation=0)
    ax1.set_xlim(*deta_range)
    ax1.set_ylim(*dphi_range)
    
    
    ax2 = fig.add_subplot(232, projection='3d')
    ax2.view_init(azim=azim,elev = elev)
    #surf = ax2.plot_surface(xpos, ypos, hist2.transpose(), cmap=viridis,edgecolor='k')
    plot2d(ax2,hist2)
    ax2.set_xlabel("$\\Delta y$")
    ax2.set_ylabel("$\\Delta\\phi$ [rad]")
    ax2.zaxis.set_rotate_label(False)
    ax2.set_title("$M(\\Delta\\phi,\\Delta y)$",rotation=0)
    ax2.set_xlim(*deta_range)
    ax2.set_ylim(*dphi_range)
    
    
    
    ax3 = fig.add_subplot(233, projection='3d')
    ax3.view_init(azim=azim,elev = elev)
    #surf = ax3.plot_surface(xpos, ypos, hist3.transpose(), cmap=viridis,edgecolor='k')
    plot2d(ax3,hist3)
    ax3.set_xlabel("$\\Delta y$")
    ax3.set_ylabel("$\\Delta\\phi$ [rad]")
    ax3.zaxis.set_rotate_label(False)
    ax3.set_title("$C(\\Delta\\phi,\\Delta y)$",rotation=0)
    ax3.set_xlim(*deta_range)
    ax3.set_ylim(*dphi_range)
    
    if not zlims is None:
        ax1.set_zlim(*(zlims[0]))
        ax2.set_zlim(*(zlims[1]))
        ax3.set_zlim(*(zlims[2]))
    
    
    bins = bins1d
    #now for 1d projections:
    ax4 = fig.add_subplot(234)
    
    if df_triggers is None:
        denom = len(df)*2*np.pi/bins
    else:
        denom = len(df_triggers)*2*np.pi/bins
    y, x = np.histogram(offset(df.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange).diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
    x = np.add(x[1:],x[:-1])/2
    dy = np.sqrt(y)
    y=np.divide(y,denom)
    dy=np.divide(dy,denom)
    
    def plot1d(x,y,dy,ax,n=5):
        ax.errorbar(x,y,dy,marker='o',color='k',linestyle='')
        ax.axhline(0,color='0.7',linestyle=':')
        #a = [2*sum(y*np.sin(i*x))/len(x) for i in range(0,n+1)]
        #remove sin term
        #a = [0 for i in range(0,n+1)]
        b = [2*sum(y*np.cos(i*x))/len(x) for i in range(0,n+1)]
        b[0]/=2
        s = 0
        xx = np.linspace(-np.pi/2,np.pi*3/2,30)
        for i in range(0,n+1):
            s = b[i]*np.cos(i*xx)+s
        ax.plot(xx,s,linestyle='-',color='k',label='all')
        for i in range(1,n+1):
            ax.plot(xx,b[i]*np.cos(i*xx)+b[0]*(i!=0),linestyle='--', label="n=%s"%i)
        ax.set_xlim(-np.pi/2,np.pi*3/2)
    plot1d(x,y,dy,ax4,n=fourierOrder)
    
    if projyrange != deta_range:
        yrangestr = "[$%s<\Delta y<%s$]" % projyrange
    else:
        yrangestr = ''
    ax4.set_xlabel("$\\Delta\\phi$ [rad]")
    ax4.set_title("$S(\\Delta\\phi)$ " + yrangestr,rotation=0)
    
    
    mixed_highy = df_mixed.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
    if(normalizeMixedAt00):
        denom = M00*2*np.pi/bins1d
    else :
        denom = len(mixed_highy)/bins1d
    ym, _ = np.histogram(offset(mixed_highy.diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
    dym = np.sqrt(ym)
    ym=np.divide(ym,denom)
    dym=np.divide(dym,denom)
    
    
    ax5 = fig.add_subplot(235)
    ax5.set_title("$M(\\Delta\\phi)$ " + yrangestr,rotation=0)
    plot1d(x,ym,dym,ax5,n=fourierOrder)
    ax5.set_xlabel("$\\Delta\\phi$ [rad]")
    
    
    ax6 = fig.add_subplot(236, sharey = ax4)
    ax6.set_title("$C(\\Delta\\phi)$ " + yrangestr,rotation=0)
    dyc = y/ym*np.hypot(dy/y, dym/ym)
    yc= y/ym
    plot1d(x,yc,dyc,ax6,n=fourierOrder)
    ax6.set_xlabel("$\\Delta\\phi$ [rad]")
    
    if not ylims is None:
        ax4.set_ylim(*(ylims[0]))
        ax5.set_ylim(*(ylims[1]))
        ax6.set_ylim(*(ylims[2]))
    
    fig.text(-0.037,0.65,text,bbox=dict(facecolor='white', alpha=1),fontsize='large')
    return fig,[ax1,ax2,ax3,ax4, ax5, ax6]


def corr1d_ratio(df,df_mixed,df_D,df_mixed_D,bins=20,style='normal',area=None,
                          fig=None,projyrange=(1.5,2.5),minbincontentindenom=3,
           label=None,color=None,showFourier=True,Ntriggers=None,NtriggersD=None,offsetX=0,alpha=1):
    dphi_range = (-np.pi/2,3*np.pi/2)
    
    #now for 1d projections:
    #ax4 = fig.add_subplot(234)
    ycs = []
    dycs = []
    for data,mixed,Ntrig in (df,df_mixed,Ntriggers),(df_D,df_mixed_D,NtriggersD):
        if Ntriggers == None:
            denom = len(data)*2*np.pi/bins
        else :
            denom = Ntriggers/bins
        y, x = np.histogram(abs(data.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange).diff_phi_cm), bins=bins, range=(0,np.pi))
        x = np.add(x[1:],x[:-1])/2
        dy = np.sqrt(y)
        y=np.divide(y,denom)
        dy=np.divide(dy,denom)

        if mixed is None:
            ycs.append(y)
            dycs.append(dy)
            continue

        mixed_highy = mixed.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)

        denom = len(mixed_highy)/bins
        ym, _ = np.histogram(abs(mixed_highy.diff_phi_cm), bins=bins, range=(0,np.pi))
        ym = ym+(ym==0)*(max(ym)*.000001)
        dym = np.sqrt(ym)
        ym=np.divide(ym,denom)
        dym=np.divide(dym,denom)

        dyc = y/ym*np.hypot(dy/y, dym/ym)
        yc= y/ym
        
        ycs.append(yc)
        dycs.append(dyc)
        #print(yc)
        #print(dycs)
    
    yc = ycs[0]/ycs[1]
    dyc = ycs[0]/ycs[1]*np.hypot(dycs[0]/ycs[0],dycs[1]/ycs[1])
    
    #print(yc)
    #print(dyc)
    
    plt.errorbar(x+offsetX,yc,dyc,color=color,linestyle='',marker='o',label=label,alpha=alpha)
    plt.gca().set_xlabel("$\\Delta\\phi$ [rad]")
    plt.gca().set_ylim(0)
    return x,yc, dyc

#create plots for uncertainties on R_2h using the mixed events
def plot_uncertainty_from_mixed(mixed_A,mixed_D,bins=20,
           label=None,color=None,offsetX=0,alpha=1, fillstyle=None,xvar='abs(diff_phi_cm)',range=(0,np.pi)):
    
    #now for 1d projections:
    #ax4 = fig.add_subplot(234)
    ys = []
    dys = []
    for mixed in (mixed_A, mixed_D):
        if mixed is None:
            continue
        y, x = np.histogram(mixed.eval(xvar), bins=bins, range=range)
        #print(y,x)
        x = np.add(x[1:],x[:-1])/2
        dy = np.sqrt(y)
        ys.append(y)
        dys.append(dy)
    
    #skip the denominator
    if len(ys) == 1:
        denom = len(mixed_A)/bins
        dy = np.divide(dy,denom)
        y = np.divide(y,denom)
        plt.errorbar(x+offsetX,y,dy,color=color,linestyle='',marker='o',label=label,alpha=alpha,fillstyle=fillstyle)
        #plt.gca().set_xlabel("$\\Delta\\phi$ [rad]")
        plt.gca().set_ylim(0)
        return x,y, dy
    
    y = ys[0]/ys[1]
    dy = ys[0]/ys[1]*np.hypot(dys[0]/ys[0],dys[1]/ys[1])
    denom = np.average(y)
    y/=denom
    dy/=denom
    
    plt.errorbar(x+offsetX,y,dy,color=color,linestyle='',marker='o',label=label,alpha=alpha,fillstyle=fillstyle)
    #plt.gca().set_xlabel("$\\Delta\\phi$ [rad]")
    plt.gca().set_ylim(0)
    
    return x,y, dy


#2d N/D plots
def NoverD_2x3(df_A,df_D,deta_range=(-1.5,2.5),bins=20,df_triggers_A=None,df_triggers_D=None,
                              fig=None,projyrange=(1.5,2.5),bins1d=50,minbincontentindenom=0.003,
               text="",label_A='',ylim=None,ylimrat=None,zlim=None,zlimrat=None,label_other=None,
              df_mix_A=None,df_mix_D=None):
    dphi_range = (-np.pi/2,3*np.pi/2)
    if fig == None:
        fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(231, projection='3d')
    
    elev = 70
    azim = -112.5#-135
    ax1.view_init(azim=azim,elev = elev)
    
    if df_triggers_A is None:
        denom = len(df_A)*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins
    else:
        denom = len(df_triggers_A)*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins
    
    #hist1,xedges, yedges = np.histogram2d([2,2,2,2],[0,0,0,0], bins=bins,range=[deta_range, dphi_range])
    hist1, xedges, yedges = np.histogram2d(df_A.diff_rap_cm, offset(df_A.diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
    hist1 = np.divide(hist1, denom)
    
    if not df_mix_A is None:
        hist1_mix, _,_ = np.histogram2d(df_mix_A.diff_rap_cm, offset(df_mix_A.diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
        hist1_mix = np.divide(hist1_mix, len(df_mix_A)/bins*(deta_range[1]-deta_range[0])/bins)
        hist1 = np.divide(hist1,hist1_mix)
    
    print(xedges)
    hist2, xedges, yedges = np.histogram2d(df_D.diff_rap_cm, offset(df_D.diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
    for i in range(bins):
        if xedges[i+1]>0:
            break
    for j in range(bins):
        if yedges[j+1]>0:
            break;
            
    radius = .3
    if df_triggers_D is None:
        denom = len(df_D)*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins
        hist2 = np.divide(hist2,denom)
    else :
        denom = len(df_triggers_D)*2*np.pi/bins*(deta_range[1]-deta_range[0])/bins
        hist2 = np.divide(hist2,denom)
    
    if not df_mix_D is None:
        hist2_mix, _,_ = np.histogram2d(df_mix_D.diff_rap_cm, offset(df_mix_D.diff_phi_cm), bins=bins, range=[deta_range, dphi_range])
        hist2_mix = np.divide(hist2_mix, len(df_mix_D)/bins*(deta_range[1]-deta_range[0])/bins)
        hist2 = np.divide(hist2,hist2_mix)
    
    
    hist3 = np.divide(hist1,np.add(hist2,0.000001))
    
    
    
    xpos, ypos = np.meshgrid(np.add(xedges[:-1],xedges[1:])/2, np.add(yedges[:-1],yedges[1:])/2)
    zpos = 0
    #help(hist1.transpose)
    
    for i in range(len(xedges)-1):
        if(xedges[i]>projyrange[0]):
            break
    print("i=",i)
    def plot2d(ax, hist):
        if deta_range != projyrange:
            ax.plot_surface(xpos.transpose()[:i+1].transpose(), ypos.transpose()[:i+1].transpose(), hist[:i+1].transpose(), cmap=viridis,edgecolor='k')
            ax.plot_surface(xpos.transpose()[i:].transpose(), ypos.transpose()[i:].transpose(), hist[i:].transpose(), cmap=autumn,edgecolor='k')
        else :
            ax.plot_surface(xpos.transpose()[:].transpose(), ypos.transpose()[:].transpose(), hist[:].transpose(), cmap=viridis,edgecolor='k')
    plot2d(ax1,hist1)
    ax1.set_xlabel("$\\Delta y$")
    ax1.set_ylabel("$\\Delta\\phi$ [rad]")
    ax1.zaxis.set_rotate_label(False)
    ax1.set_title(label_A,rotation=0)
    ax1.set_xlim(*deta_range)
    ax1.set_ylim(*dphi_range)
    
    
    ax2 = fig.add_subplot(232, projection='3d',sharez = ax1)
    ax2.view_init(azim=azim,elev = elev)
    #surf = ax2.plot_surface(xpos, ypos, hist2.transpose(), cmap=viridis,edgecolor='k')
    plot2d(ax2,hist2)
    ax2.set_xlabel("$\\Delta y$")
    ax2.set_ylabel("$\\Delta\\phi$ [rad]")
    ax2.zaxis.set_rotate_label(False)
    ax2.set_title("D",rotation=0)
    ax2.set_xlim(*deta_range)
    ax2.set_ylim(*dphi_range)
    if(zlim == None):
        zlim=ylim
    ax2.set_zlim(*zlim)
    
    
    ax3 = fig.add_subplot(233, projection='3d')
    ax3.view_init(azim=azim,elev = elev)
    #surf = ax3.plot_surface(xpos, ypos, hist3.transpose(), cmap=viridis,edgecolor='k')
    plot2d(ax3,hist3)
    ax3.set_xlabel("$\\Delta y$")
    ax3.set_ylabel("$\\Delta\\phi$ [rad]")
    ax3.zaxis.set_rotate_label(False)
    ax3.set_title(f"{label_A}/D",rotation=0)
    ax3.set_xlim(*deta_range)
    ax3.set_ylim(*dphi_range)
    if(zlimrat == None):
        zlimrat=ylimrat
    ax3.set_zlim(*zlimrat)
    
    
    bins = bins1d
    #now for 1d projections:
    ax4 = fig.add_subplot(234)
    
    df_A = df_A.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
    if df_triggers_A is None:
        denom = len(df_A)*2*np.pi/bins
    else:
        denom = len(df_triggers_A)*2*np.pi/bins
    y, x = np.histogram(offset(df_A.diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
    x = np.add(x[1:],x[:-1])/2
    dy = np.sqrt(y)
    y=np.divide(y,denom)
    dy=np.divide(dy,denom)
    
    if not df_mix_A is None:
        df_mix_A = df_mix_A.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
        y_mix,_ = np.histogram(offset(df_mix_A.diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
        y_mix = np.divide(y_mix,len(df_mix_A)/bins1d)
        y = y/y_mix
    
    
    
    def plot1d(x,y,dy,ax,n=5):
        ax.errorbar(x,y,dy,marker='o',color='k',linestyle='')
        ax.axhline(0,color='0.7',linestyle=':')
        #a = [2*sum(y*np.sin(i*x))/len(x) for i in range(0,n+1)]
        #remove sin term
        #a = [0 for i in range(0,n+1)]
        b = [2*sum(y*np.cos(i*x))/len(x) for i in range(0,n+1)]
        b[0]/=2
        s = 0
        for i in range(0,n+1):
            s = b[i]*np.cos(i*x)+s
        ax.plot(x,s,linestyle='-',color='k',label='all')
        for i in range(1,n+1):
            ax.plot(x,b[i]*np.cos(i*x)+b[0]*(i!=0),linestyle='--', label="n=%s"%i)
        ax.set_xlim(-np.pi/2,np.pi*3/2)
    plot1d(x,y,dy,ax4)
    
    if deta_range != projyrange:
        yrangestr = "[$%s<\Delta y<%s$]" % projyrange
    else:
        yrangestr = ''
    ax4.set_xlabel("$\\Delta\\phi$ [rad]")
    ax4.set_title(f"{label_A} " + yrangestr,rotation=0)
    
    
    df_D = df_D.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
    if df_triggers_D is None:
        denom = len(df_D)*2*np.pi/bins1d
    else :
        denom = len(df_triggers_D)*2*np.pi/bins1d
    yd, _ = np.histogram(offset(df_D.diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
    dyd = np.sqrt(yd)
    yd=np.divide(yd,denom)
    dyd=np.divide(dyd,denom)
    
    
    if not df_mix_D is None:
        df_mix_D = df_mix_D.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
        y_mix,_ = np.histogram(offset(df_mix_D.diff_phi_cm), bins=bins, range=(-np.pi/2,3*np.pi/2))
        y_mix = np.divide(y_mix,len(df_mix_D)/bins1d)
        yd = yd/y_mix
    
    
    ax5 = fig.add_subplot(235, sharey = ax4)
    ax5.set_title("D " + yrangestr,rotation=0)
    plot1d(x,yd,dyd,ax5)
    ax5.set_xlabel("$\\Delta\\phi$ [rad]")
    ax5.set_ylim(*ylim)
    
    
    ax6 = fig.add_subplot(236)
    ax6.set_title(f"{label_A}/D " + yrangestr,rotation=0)
    dyc = y/yd*np.hypot(dy/y, dyd/yd)
    yc= y/yd
    plot1d(x,yc,dyc,ax6)
    ax6.set_xlabel("$\\Delta\\phi$ [rad]")
    ax6.set_ylim(*ylimrat)
    
    if not label_other is None:
        plt.text(0, .9, label_other, fontsize='x-large', transform = ax1.transAxes)
    
    if not text is None:
        fig.text(-0.037,0.65,text,bbox=dict(facecolor='white', alpha=1),fontsize='large')
    return fig,[ax1,ax2,ax3,ax4, ax5, ax6]


# 1D plot   Setting df_D to None plots just the nuclear data. 
# Setting the df_D argument plots the ratio of nuclear to deuterium
# Setting the triggers variables to None skips the division by the number of trigger events
# Setting the mixed df variables creates correlation functions.  leaving it None only does the same-event yield
def general_1D_plot(df_A,df_D,df_triggers_A=None,df_triggers_D=None,
                              fig=None,projyrange=(0,2.5),bins1d=10,
              df_mix_A=None,df_mix_D=None,color=None,label=None, plotFourier=False,linestyleFit='-',fillstyle='full',
                   marker='o',markersize=7):
    dphi_range = (0,2*np.pi)
    
    
    
    bins = bins1d
    #now for 1d projections:
    
    df_A = df_A.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
    if df_triggers_A is None:
        denom = len(df_A)*2*np.pi/bins
    else:
        denom = len(df_triggers_A)*2*np.pi/bins
    y, x = np.histogram(offset(df_A.diff_phi_cm,np.pi), bins=bins, range=(0,2*np.pi))
    x = np.add(x[1:],x[:-1])/2
    dy = np.sqrt(y)
    y=np.divide(y,denom)
    dy=np.divide(dy,denom)
    
    if not df_mix_A is None:
        df_mix_A = df_mix_A.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
        y_mix,_ = np.histogram(offset(df_mix_A.diff_phi_cm,np.pi), bins=bins, range=(0,2*np.pi))
        y_mix = np.divide(y_mix,len(df_mix_A)/bins1d)
        y = y/y_mix
        dy = dy/y_mix
    if df_D is None:
        plt.errorbar(x,y,dy,marker='o',color=color,linestyle='',label=label,fillstyle=fillstyle)
        if plotFourier:
            n = plotFourier
            b = [2*sum(y*np.cos(i*x))/len(x) for i in range(0,n+1)]
            b[0]/=2
            s = 0
            
            x = np.linspace(0,2*np.pi,50)
            for i in range(0,n+1):
                s = b[i]*np.cos(i*x)+s
            plt.plot(x,s,linestyle=linestyleFit,color=color)
        return
    
    df_D = df_D.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
    if df_triggers_D is None:
        denom = len(df_D)*2*np.pi/bins1d
    else :
        denom = len(df_triggers_D)*2*np.pi/bins1d
    yd, _ = np.histogram(offset(df_D.diff_phi_cm,np.pi), bins=bins, range=(0,2*np.pi))
    dyd = np.sqrt(yd)
    yd=np.divide(yd,denom)
    dyd=np.divide(dyd,denom)
    
    
    if not df_mix_D is None:
        df_mix_D = df_mix_D.query("diff_rap_cm > %s and diff_rap_cm < %s" %projyrange)
        y_mix,_ = np.histogram(offset(df_mix_D.diff_phi_cm,np.pi), bins=bins, range=(0,2*np.pi))
        y_mix = np.divide(y_mix,len(df_mix_D)/bins1d)
        yd = yd/y_mix
    
    
    dyc = y/yd*np.hypot(dy/y, dyd/yd)
    yc= y/yd
    
    plt.errorbar(x,yc,dyc,marker='o',color=color,linestyle='',label=label,fillstyle=fillstyle)
    
    if plotFourier:
        n = plotFourier
        b = [2*sum(yc*np.cos(i*x))/len(x) for i in range(0,n+1)]
        b[0]/=2
        s = 0
        x = np.linspace(0,2*np.pi,50)
        for i in range(0,n+1):
            s = b[i]*np.cos(i*x)+s
        plt.plot(x,s,linestyle=linestyleFit,color=color)