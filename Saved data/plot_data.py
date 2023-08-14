# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 12:17:59 2021

@author: nandan
"""


import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pylab as pylab

params = {'legend.fontsize': 15,
         'axes.labelsize': 'x-large',
         'axes.labelpad' : 15,
         'axes.titlesize':'x-large',
         'xtick.labelsize':20,
         'ytick.labelsize':20}
pylab.rcParams.update(params)





def threshold_values(model,experiment):
    filename=model+'_'+experiment
    output=np.load(os.path.join(folder_save_output,filename+'.npy'))
    
    Percentage_difference_mean=output[0]
    Polarizarion_ratio_mean=output[1]
    Polarizarion_time_mean=output[3]
    Polarizarion_time_std=output[4]
    
    max_ratio=np.max(Polarizarion_ratio_mean)
    threshold_indx=np.argwhere((Polarizarion_ratio_mean>0.7*max_ratio)==True)[0][0]
    threshold_pd=Percentage_difference_mean[threshold_indx]
    threshold_time=Polarizarion_time_mean[threshold_indx]
    threshold_time_std=Polarizarion_time_std[threshold_indx]
    if model=='Legi':
        # threshold_pd=Percentage_difference_mean[0]
        # threshold_time=Polarizarion_time_mean[np.argwhere(Percentage_difference_mean==threshold_pd)][0][0]
        # threshold_time_std=Polarizarion_time_std[np.argwhere(Percentage_difference_mean==threshold_pd)][0][0]
        
        threshold_pd=Percentage_difference_mean[1]
        threshold_time=Polarizarion_time_mean[np.argwhere(Percentage_difference_mean==threshold_pd)][0][0]
        threshold_time_std=Polarizarion_time_std[np.argwhere(Percentage_difference_mean==threshold_pd)][0][0]
    return threshold_pd,threshold_time,threshold_time_std

def reversal_time(model, experiment):
    filename=model+'_'+experiment
    output=np.load(os.path.join(folder_save_output,filename+'.npy'))
    
    stimulus_ratio_mean=output[0]
    Polarizarion_amplification_mean=output[1]
    RePolarizarion_time_mean=output[3]
    RePolarizarion_time_std=output[4]
    
    pa_1=Polarizarion_amplification_mean[0]
    rpt_1=RePolarizarion_time_mean[0]
    rpt_1_std=RePolarizarion_time_std[0]
    
    pa_2=Polarizarion_amplification_mean[-1]
    rpt_2=RePolarizarion_time_mean[-1]
    rpt_2_std=RePolarizarion_time_std[-1]
    
    return [pa_1,pa_2,rpt_1,rpt_2,rpt_1_std,rpt_2_std]

def resolving_time(model, experiment):
    filename=model+'_'+experiment
    output=np.load(os.path.join(folder_save_output,filename+'.npy'))
    
    stimulus_ratio_mean=output[0]
    Polarizarion_amplification_mean=output[1]
    Resolving_time_mean=output[3]
    Resolving_time_std=output[4]
    
    pa_1=Polarizarion_amplification_mean[0]
    rt_1=Resolving_time_mean[0]
    rt_1_std=Resolving_time_std[0]
    
    pa_2=Polarizarion_amplification_mean[-1]
    rt_2=Resolving_time_mean[-1]
    rt_2_std=Resolving_time_std[-1]
    
    return [pa_1,pa_2,rt_1,rt_2,rt_1_std,rt_2_std]


def bar_plot(data_mean,data_std,ylabel,ylimit,yticks,save=None):
    
    xlabels=['Legi','Turing','WavePinning','SubPB']
    
    plt.figure()
    if data_std is not None:
        plt.bar([1,2,3,4],data_mean,width=0.5,color=Colors,yerr=data_std)
    else:
        plt.bar([1,2,3,4],data_mean,width=0.5,color=Colors,yerr=data_std)
    plt.xticks(ticks=[1,2,3],labels=xlabels,fontsize=20,fontname='Arial')
    plt.xticks(rotation=45,ha='right')
    plt.yticks(ticks=yticks,fontsize=15,fontname='Arial')
    plt.xlim(0,5);plt.ylim(ylimit[0],ylimit[1])
    plt.ylabel(ylabel,fontsize=20,fontname='Arial')
    # plt.title(data_mean)
    
    if save==True:
        plt.savefig(os.path.join(folder_save_fig,'all_model_type_'+experiment+'_bar_plot.svg'))
    plt.show()

def bar_plot_composite(data_x,data_y_mean,data_y_std,width,ylabel,ylimit,yticks,xlimit,save=None):
    
    xlabels=['Legi','Turing','WavePinning','SubPB']
    
    plt.figure()
    if data_y_std is not None:
        plt.bar(data_x,data_y_mean,width=width,color=Colors,yerr=data_y_std)
    else:
        plt.bar(data_x,data_y_mean,width=width,color=Colors,yerr=data_y_std)
    
    plt.xticks(ticks=data_x,labels=np.round(data_x,2),fontsize=20,fontname='Arial')
    plt.xlabel('Polarization amplification',fontsize=20,fontname='Arial')
    # plt.xticks(rotation=45,ha='right')
    plt.yticks(ticks=yticks,fontsize=20,fontname='Arial')
    plt.xlim(xlimit[0],xlimit[1]);plt.ylim(ylimit[0],ylimit[1])
    plt.ylabel(ylabel,fontsize=20,fontname='Arial')
    # plt.title(data_mean)
    
    if save==True:
        plt.savefig(os.path.join(folder_save_fig,'all_model_type_'+experiment+'_bar_plot.svg'))
    
    plt.show()

def bar_plot3D_resolving(data_x,data_y_mean):
    xlabels=['Legi','Turing','WavePinning','SubPB']
    
    x3 = [0.75,1.75,2.75,3.75]
    y3 = data_x
    z3 = np.zeros(4)
    
    dx = np.ones(4)*0.5
    dy = np.ones(4)*1
    dz = data_y_mean
 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(x3, y3, z3, dx, dy, dz,color=Colors) 
    
    ax.w_xaxis.set_ticklabels(xlabels,fontsize=20,fontname='Arial',rotation=45)
    ax.set_ylabel('Polarization amplification',fontsize=15,fontname='Arial')
    ax.set_zlabel('Resolving time(sec)',fontsize=15,fontname='Arial')
    
    # ax.xaxis.pane.fill = False
    # ax.yaxis.pane.fill = False
    # ax.zaxis.pane.fill = False
    
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    ax.grid(color='k',lw=1)
        
    plt.show()
    
def bar_plot3D_reversing(data_x,data_y_mean):
    xlabels=['Legi','Turing','WavePinning','SubPB']
    
    x3 = [0.75,1.75,2.75,3.75]
    y3 = data_x
    z3 = np.zeros(4)
    
    dx = np.ones(4)*0.5
    dy = np.ones(4)*0.0015
    dz = data_y_mean
 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(x3, y3, z3, dx, dy, dz,color=Colors) 
    
    ax.w_xaxis.set_ticklabels(xlabels,fontsize=20,fontname='Arial',rotation=45)
    ax.set_ylabel('Polarization amplification',fontsize=15,fontname='Arial')
    ax.set_zlabel('Reversal time(sec)',fontsize=15,fontname='Arial')
    
    # ax.xaxis.pane.fill = False
    # ax.yaxis.pane.fill = False
    # ax.zaxis.pane.fill = False
    
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    ax.grid(color='k',lw=1)
    
    ax.set_ylim(1,1.03)
    ax.set_zlim(0,500)
        
    plt.show()
    

def colocalization_plot(model, experiment):
    filename=model+'_'+experiment+'_colocalization'
    output=np.load(os.path.join(folder_save_output,filename+'.npy'))
    
    N,M=np.shape(output)
    
    colocalization_index=np.zeros(M)*np.nan
    pd=np.arange(0,0.05+0.001,0.001)*100
    for m in range(M):
        if len(np.where(np.isnan(output[:,m]))[0])==N:
            colocalization_index[m]=np.nan
        else:
            colocalization_index[m]=len(np.where(output[:,m]==1)[0])/N
    
    plt.figure()
    plt.plot(pd,colocalization_index,'ko')
    plt.plot(pd,colocalization_index,'k-')
    plt.xticks(ticks=np.arange(0,5+0.1,0.2),fontsize=15,fontname='Arial')
    plt.xlim(0,5)  
    plt.yticks(ticks=[0,0.2,0.4,0.6,0.8,1],fontsize=15,fontname='Arial')
    plt.xlabel('Percentage difference(%)',fontsize=20,fontname='Arial')
    plt.ylabel('Colocalization-index',fontsize=20,fontname='Arial')
    plt.ylim(0,1.1)
    plt.show()
    
  
###############################################################

# folder_save_output='D:\\polarity paper\\PNAS\\output\\'
folder_save_output=os.getcwd()+'\\final_output\\'
folder_save_fig=folder_save_output+'figs\\'    

#Models=['Legi']
#Experiments=['single_gradient']
#Colors=['black']

Models=['Legi','Otsuji','WavePinning','SubPB']

Experiments=['single_gradient'] 
# Experiments=['gradient_reversal']
# Experiments=['simultaneous_gradients']
# Experiments=['single_gradient_with_offset']

Colors=['magenta','green','blue','red']

save_fig=None

Alphas=[0.7,0.7]

plt.figure()

for j, experiment in enumerate(Experiments):
    for i, model in enumerate(Models):
        filename=model+'_'+experiment
        output=np.load(os.path.join(folder_save_output,filename+'.npy'))
        
        x1_mean=output[0]
        y1_mean=output[1]
        y1_std=output[2]
        y2_mean=output[3]
        y2_std=output[4]
        
        ### plot Fig. 3A
        plt.plot(x1_mean,y1_mean,color=Colors[i],lw=1)
        plt.errorbar(x1_mean,y1_mean,yerr=y1_std,linewidth=1,fmt='s-',color=Colors[i],alpha=Alphas[j])
        plt.xlim(0,5)
        plt.xticks(fontsize=15,fontname='Arial')
        plt.xlabel('stimulus difference(sd, %)',fontsize=20,fontname='Arial')
        plt.ylim(1,10)
        plt.yticks(ticks=[0,5,10],fontsize=15,fontname='Arial')
        plt.ylabel('polarization ratio',fontsize=20,fontname='Arial')
        

### single gradient and localized

if experiment=='single_gradient' or experiment=='localized' or experiment=='single_gradient_with_offset':
  
    if experiment=='single_gradient' :
        
        
        thresh_legi,thresh_time_legi,thresh_time_std_legi=threshold_values('Legi',experiment=Experiments[0])
        thresh_otsuji,thresh_time_otsuji,thresh_time_std_otsuji=threshold_values('Otsuji',experiment=Experiments[0])
        thresh_wavepinning,thresh_time_wavepinning,thresh_time_std_wavepinning=threshold_values('WavePinning',experiment=Experiments[0])
        thresh_subpb,thresh_time_subpb,thresh_time_std_subpb=threshold_values('SubPB',experiment=Experiments[0])
        
        
        # plt.title('Persistent gradient vs localized',fontsize=20,fontname='Arial')
        
        
        if save_fig==True:
            plt.savefig(os.path.join(folder_save_fig,'all_model_type_'+experiment+'_pr.svg'))
        
        data_mean_ta=np.array([thresh_legi,thresh_otsuji,thresh_wavepinning,thresh_subpb])
        data_mean_pt=np.array([thresh_time_legi,thresh_time_otsuji,thresh_time_wavepinning,thresh_time_subpb])
        data_std_pt=np.array([thresh_time_std_legi,thresh_time_std_otsuji,thresh_time_std_wavepinning,thresh_time_std_subpb])
        
        ### plot Fig. 3B
        plt.figure()
        plt.errorbar(data_mean_ta,data_mean_pt,yerr=data_std_pt,fmt='o',ecolor='k',elinewidth=1,capsize=5,ms=10,color='k')
        plt.xticks(ticks=np.arange(0,1.5,0.1),fontsize=15,fontname='Arial')  
        plt.yticks(ticks=np.arange(0,1200,200),fontsize=15,fontname='Arial')
        plt.xlabel(r'threshold for activation($sd_{thresh}$, %)',fontsize=20,fontname='Arial')
        plt.ylabel('polarization time(sec)',fontsize=20,fontname='Arial')
        plt.ylim(0,1100)
        plt.xlim(-0.1,1.5) 
        plt.xticks(ticks=[0,0.5,1],fontsize=15,fontname='Arial')
        if save_fig==True:
            plt.savefig(os.path.join(folder_save_fig,'all_model_type_'+experiment+'_bar_plot.svg'))
        plt.show()
    
else:
    
    plt.xlim(1,2)
    plt.xlabel('Stimulus ratio',fontsize=20,fontname='Arial')
    # plt.ylim(0.99,1.04)
    # if experiment=='gradient_reversal':
    #     plt.yticks(ticks=[0.99,1,1.02,1.04],fontsize=20,fontname='Arial')
    # elif experiment=='simultaneous_gradients':
    #     plt.yticks(ticks=[1,5,10,15,20],fontsize=20,fontname='Arial')   
    
    plt.ylabel('Polarization amplification',fontsize=20,fontname='Arial')
    
    
    # plt.xlim(1,2)
    # plt.xlabel('Stimulus ratio',fontsize=20,fontname='Arial')
    # plt.ylim(0,700)
    # plt.yticks(ticks=[0,100,200,300,400,500,600,700],fontsize=15,fontname='Arial')
    # plt.ylabel('Re-polarization time(sec)',fontsize=20,fontname='Arial')
    
    # plt.show()
    
    # if save_fig==True:
    #     plt.savefig(os.path.join(folder_save_fig,'all_model_type_'+experiment+'.svg'))
    
    ######## Reversing stimulus
    
    if experiment=='gradient_reversal':
    
        [pa_1_legi,pa_2_legi,rpt_1_legi,rpt_2_legi,rpt_1_std_legi,rpt_2_std_legi] = reversal_time('Legi','gradient_reversal')
        [pa_1_otsuji,pa_2_otsuji,rpt_1_otsuji,rpt_2_otsuji,rpt_1_std_otsuji,rpt_2_std_otsuji] = reversal_time('Otsuji','gradient_reversal')
        [pa_1_wavepinning,pa_2_wavepinning,rpt_1_wavepinning,rpt_2_wavepinning,rpt_1_std_wavepinning,rpt_2_std_wavepinning] = reversal_time('WavePinning','gradient_reversal')
        [pa_1_subpb,pa_2_subpb,rpt_1_subpb,rpt_2_subpb,rpt_1_std_subpb,rpt_2_std_subpb] = reversal_time('SubPB_new','gradient_reversal')
        
        data_mean_pa_1=np.array([pa_1_legi,pa_1_otsuji,pa_1_wavepinning,pa_1_subpb])
        data_mean_pa_2=np.array([pa_2_legi,pa_2_otsuji,pa_2_wavepinning,pa_2_subpb])
        data_mean_rpt_1=np.array([rpt_1_legi,rpt_1_otsuji,rpt_1_wavepinning,rpt_1_subpb])
        data_mean_rpt_2=np.array([rpt_2_legi,rpt_2_otsuji,rpt_2_wavepinning,rpt_2_subpb])
        data_mean_rpt_1_std=np.array([rpt_1_std_legi,rpt_1_std_otsuji,rpt_1_std_wavepinning,rpt_1_std_subpb])
        data_mean_rpt_2_std=np.array([rpt_2_std_legi,rpt_2_std_otsuji,rpt_2_std_wavepinning,rpt_2_std_subpb])       
        
        bar_plot_composite(data_mean_pa_2,data_mean_rpt_2,data_y_std=data_mean_rpt_2_std,ylabel='Reversing time(sec)',width=0.05,ylimit=[0,500],yticks=np.arange(0,500,100),xlimit=[data_mean_pa_2[0]-0.1,data_mean_pa_2[-1]+0.1],save=save_fig)
        # bar_plot3D_reversing(data_mean_pa_2,data_mean_rpt_2)
        
        plt.figure()
        plt.plot(data_mean_pa_2,data_mean_rpt_2,color=Colors[i],lw=1)
        plt.errorbar(data_mean_pa_2,data_mean_rpt_2,yerr=data_mean_rpt_2_std,fmt='o',ecolor='k',elinewidth=1,capsize=5,ms=10,color='k')
        plt.xticks(ticks=np.arange(1,1.3,0.1),fontsize=15,fontname='Arial')  
        plt.yticks(ticks=np.arange(0,500,100),fontsize=15,fontname='Arial')
        plt.xlabel('Polarization amplification',fontsize=20,fontname='Arial')
        plt.ylabel('Reversing time(sec)',fontsize=20,fontname='Arial')
        plt.xlim(0.95,1.35)
        plt.ylim(0,500)
        plt.show()
    
    
    ######## Resolving stimulus
    
    if experiment=='simultaneous_gradients':
    
        [pa_1_legi,pa_2_legi,rt_1_legi,rt_2_legi,rt_1_std_legi,rt_2_std_legi] = resolving_time('Legi','simultaneous_gradients')
        [pa_1_otsuji,pa_2_otsuji,rt_1_otsuji,rt_2_otsuji,rt_1_std_otsuji,rt_2_std_otsuji] = resolving_time('Otsuji','simultaneous_gradients')
        [pa_1_wavepinning,pa_2_wavepinning,rt_1_wavepinning,rt_2_wavepinning,rt_1_std_wavepinning,rt_2_std_wavepinning] = resolving_time('WavePinning','simultaneous_gradients')
        [pa_1_subpb,pa_2_subpb,rt_1_subpb,rt_2_subpb,rt_1_std_subpb,rt_2_std_subpb] = resolving_time('SubPB_new','simultaneous_gradients')
         
        data_mean_pa_1=np.array([pa_1_legi,pa_1_otsuji,pa_1_wavepinning,pa_1_subpb])
        data_mean_pa_2=np.array([pa_2_legi,pa_2_otsuji,pa_2_wavepinning,pa_2_subpb])
        data_mean_rpt_1=np.array([rt_1_legi,rt_1_otsuji,rt_1_wavepinning,rt_1_subpb])
        data_mean_rpt_2=np.array([rt_2_legi,rt_2_otsuji,rt_2_wavepinning,rt_2_subpb]) 
        data_mean_rpt_1_std=np.array([rt_1_std_legi,rt_1_std_otsuji,rt_1_std_wavepinning,rt_1_std_subpb])
        data_mean_rpt_2_std=np.array([rt_2_std_legi,rt_2_std_otsuji,rt_2_std_wavepinning,rt_2_std_subpb])
        
        # bar_plot(data_mean_pa_1,data_std=None,ylabel='Polarizarion amplification',ylimit=[0.99,20],yticks=[1,5,10,15,20])
        # bar_plot(data_mean_rpt_1,data_std=data_mean_rpt_1_std,ylabel='Resolving time(sec)',ylimit=[0,1000],yticks=np.arange(0,1000,100))
    
        # bar_plot(data_mean_pa_2,data_std=None,ylabel='Polarizarion amplification',ylimit=[0.99,20],yticks=[1,5,10,15,20])
        # bar_plot(data_mean_rpt_2,data_std=data_mean_rpt_2_std,ylabel='Resolving time(sec)',ylimit=[0,1000],yticks=np.arange(0,1000,100))
        
        bar_plot_composite(data_mean_pa_2,data_mean_rpt_2,data_y_std=data_mean_rpt_2_std,ylabel='Resolving time(sec)',width=0.5,ylimit=[0,1000],yticks=np.arange(0,1000,100),xlimit=[data_mean_pa_2[0]-0.5,np.max(data_mean_pa_2)+0.5],save=save_fig)
        # bar_plot3D(data_mean_pa_2,data_mean_rpt_2)
        
        plt.figure()
        plt.plot(data_mean_pa_2,data_mean_rpt_2,color=Colors[i],lw=1)
        plt.errorbar(data_mean_pa_2,data_mean_rpt_2,yerr=data_mean_rpt_2_std,fmt='o',ecolor='k',elinewidth=1,capsize=5,ms=10,color='k')
        plt.xticks(ticks=np.arange(0,8,1),fontsize=15,fontname='Arial')  
        plt.yticks(ticks=np.arange(0,1000,200),fontsize=15,fontname='Arial')
        plt.xlabel('Polarization amplification',fontsize=20,fontname='Arial')
        plt.ylabel('Reversing time(sec)',fontsize=20,fontname='Arial')
        # plt.xlim(0.95,1.35)
        # plt.ylim(0,500)
        plt.show()
