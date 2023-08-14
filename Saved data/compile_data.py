# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 16:58:57 2022

@author: nandan
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from colormaps import parula_map
import matplotlib.pylab as pylab

params = {'legend.fontsize': 15,
         'axes.labelsize': 20,
         'axes.labelpad' : 15,
         'axes.titlesize':'x-large',
         'xtick.labelsize':20,
         'ytick.labelsize':20}
pylab.rcParams.update(params)

def find_stimulus_minmax(stimulus):
    max_indx=np.argmax(stimulus[:,-1]) 
    min_indx=np.argmin(stimulus[:,-1]) 
    
    indxs_with_max=np.where(stimulus[:,-1]==np.max(stimulus[:,-1]))[0]
    if len(indxs_with_max)==1:
        pass
    else:
        max_indx=indxs_with_max[0]+len(indxs_with_max)//2

    return max_indx, min_indx


def find_response_minmax(response,time):
    max_indx=np.argmax(response[:,time]) 
    min_indx=np.argmin(response[:,time]) 
    
    indxs_with_max=np.where(response[:,time]==np.max(response[:,time]))[0]
    
    if len(indxs_with_max)==1: # check whether response max is a plateau
        pass
    else:
        max_indx=indxs_with_max[0]+len(indxs_with_max)//2

    return max_indx, min_indx   

def max_min(data, not_stimulus=None):
    
    if not_stimulus is not None:
        data=(data-np.min(data))/(np.max(data)-np.min(data))
        
        indxs_with_max=np.where(data>0.95)[0]
        indxs_with_min=np.where(data<0.05)[0]     
    else:
        indxs_with_max=np.where(data==np.max(data))[0]
        indxs_with_min=np.where(data==np.min(data))[0]
    
    if len(indxs_with_max)==1: # check whether response max is a plateau
        max_indx=indxs_with_max[0]
        pass
    else:
        if len(np.where(indxs_with_max==0)[0])==0:
            max_indx=indxs_with_max[0]+len(indxs_with_max)//2 
        else: 
            indxs_with_max_new=indxs_with_max[np.where(indxs_with_max<=10)]

            rest=indxs_with_max[np.where(indxs_with_max>10)]
            rest_N=len(rest)
            indxs_with_max_new=np.pad(indxs_with_max_new,(rest_N,0))
            indxs_with_max_new[:rest_N]=rest
            
            max_indx=indxs_with_max_new[len(indxs_with_max_new)//2]
            
        
    if len(indxs_with_min)==1: # check whether response max is a plateau
        min_indx=indxs_with_min[0]
        pass
    else:
        
        if len(np.where(indxs_with_min==0)[0])==0:
            min_indx=indxs_with_min[0]+len(indxs_with_min)//2
 
        else:
            indxs_with_min_new=indxs_with_min[np.where(indxs_with_min<=10)]
            
            rest=indxs_with_min[np.where(indxs_with_min>10)]
            rest_N=len(rest)
            indxs_with_min_new=np.pad(indxs_with_min_new,(rest_N,0))
            indxs_with_min_new[:rest_N]=rest
            
            min_indx=indxs_with_min_new[len(indxs_with_min_new)//2]

    return max_indx, min_indx   


def time_series(response, stimulus, tp):
    
    N=len(response[:,0])
    
    # sf_indx=np.argmax(stimulus[:,tp]) 
    # sb_indx=np.argmin(stimulus[:,tp])
    # rf_indx=np.argmax(response[:,tp]) 
    # rb_indx=np.argmin(response[:,tp])
    

    sf_indx, sb_indx = max_min(stimulus[:,tp])
    # sb_indx=10
    rf_indx, rb_indx = max_min(response[:,tp], not_stimulus=True)
    # rb_indx=10

    
    if sb_indx==0:
        sf_set=response[:,tp][[sf_indx-1,sf_indx,sf_indx+1]]
        sb_set=response[:,tp][[N-1,sb_indx,sb_indx+1]]
    elif sb_indx==N-1:
        sf_set=response[:,tp][[sf_indx-1,sf_indx,sf_indx+1]]
        sb_set=response[:,tp][[sb_indx-1,sb_indx,0]]
    elif sf_indx==0:
        sf_set=response[:,tp][[N-1,sf_indx,sf_indx+1]]
        sb_set=response[:,tp][[sb_indx-1,sb_indx,sb_indx+1]]
    elif sf_indx==N-1:
        sf_set=response[:,tp][[sf_indx-1,sf_indx,0]]
        sb_set=response[:,tp][[sb_indx-1,sb_indx,sb_indx+1]]
    else:    
        sf_set=response[:,tp][[sf_indx-1,sf_indx,sf_indx+1]]
        sb_set=response[:,tp][[sb_indx-1,sb_indx,sb_indx+1]]
    sf_set = np.sort(sf_set)
    sb_set = np.sort(sb_set)
    
    if (sf_set[0]<=response[:,tp][rf_indx]<=sf_set[-1]) and ((sb_set[0]<=response[:,tp][rb_indx]<=sb_set[-1])):
        max_indx=rf_indx
        min_indx=rb_indx
        
        reponse_front=response[max_indx]
        response_back=response[min_indx]
        alignment_type=True
    else:
        max_indx=np.nan
        min_indx=np.nan
        
        reponse_front=response[rf_indx]
        response_back=response[rb_indx]
        
        alignment_type=False
      
    
    return reponse_front, response_back, alignment_type

############################################################### single gradient quantities
def polarization_time(response_front,response_back,tol):
    
    # normalization of timeseries between min and max
    res_max=np.max(response_front)
    res_min=np.min(response_back)
    response_front=(response_front-res_min)/(res_max-res_min)
    response_back=(response_back-res_min)/(res_max-res_min)
    
    # front_ss=response_front[-1]
    # back_ss=response_back[-1]
    
    front_ss=np.mean(response_front[-10:])
    back_ss=np.mean(response_back[-10:])
    
    response_front_delta=abs(front_ss-response_front)
    response_back_delta=abs(response_back-back_ss)
    
    ss_indx1=np.argwhere((response_front_delta<tol)==True)
    ss_indx2=np.argwhere((response_back_delta<tol)==True)
    front_t= ss_indx1[0][0]
    back_t= ss_indx2[0][0]   
    pt=np.max([front_t,back_t])   
    return pt


def polarization_ratio(time,response_front,response_back,stimulus,pt):
    pol_ind=np.where(time==pt)[0][0]
    pol_max=np.mean(response_front[pol_ind:-1])
    pol_min=np.mean(response_back[pol_ind:-1])
        
    stimulus_max=np.max(stimulus[:,-1])
    stimulus_min=np.min(stimulus[:,-1])
    
    pd=(stimulus_max-stimulus_min)*100
    # fc=int(stimulus_max/(stimulus_max-stimulus_min))
    fc=(stimulus_min/(stimulus_max-stimulus_min))
    return pd,pol_max,pol_min,fc

################################################################ gradient reversal quantities

def reversal_time(response_front,response_back,tol,treversal):
    
    # normalization of timeseries between min and max
    res_max=np.max(response_front)
    res_min=np.min(response_back)
    response_front=(response_front-res_min)/(res_max-res_min)
    response_back=(response_back-res_min)/(res_max-res_min)
    
    # front_ss=response_front[-1]
    # back_ss=response_back[-1]
    
    front_ss=np.mean(response_front[-100:])
    back_ss=np.mean(response_back[-100:])
    
    # front_ss=response_front[-1]
    # back_ss=response_back[-1]
    if front_ss>back_ss:
        rt=np.nan
    else:
        response_front_delta=abs(front_ss-response_front)
        response_back_delta=abs(response_back-back_ss)
        
        ss_indx1=np.argwhere((response_front_delta<tol)==True)
        ss_indx2=np.argwhere((response_back_delta<tol)==True)
        front_t= ss_indx1[0][0]
        back_t= ss_indx2[0][0]   
        rt=np.max([front_t,back_t]) 
    return rt


def polarization_amplification_reversal(response_front,response_back,stimulus,treversal):
    
    sfront=np.max(stimulus[:,treversal])
    sfront_indx=np.argmax(stimulus[:,treversal])
    
    srev=np.max(stimulus[:,-1])
    srev_indx=np.argmax(stimulus[:,-1])
    
    drg=srev/sfront
    
    ufront=response_front[treversal]
    uback=response_back[-1]
    
    pa=uback/ufront
    # pa = (abs(UR-UL))*100
    
    
    
    return drg,pa

############################################################### resolving quantities


def resolving_quantities(response_front,response_back,stimulus, thesh_stimu, tol_time):
    
    Sf=np.max(stimulus[:,-1])
    drg=Sf/thesh_stimu
    
    # estimating the actual polarization apmplification
    front_ss_wt=np.mean(response_front[-100:])
    back_ss_wt=np.mean(response_back[-100:])
    pa_wt=front_ss_wt/back_ss_wt
    
   # normalization of timeseries between min and max
    res_max=np.max(response_front)
    res_min=np.min(response_back)
    response_front=(response_front-res_min)/(res_max-res_min)
    response_back=(response_back-res_min)/(res_max-res_min)
    
    # front_ss=response_front[-1]
    # back_ss=response_back[-1]
    
    front_ss=np.mean(response_front[-1:])
    back_ss=np.mean(response_back[-1:])
    
    # if abs(front_ss-back_ss)*100<tol_res: # if front and back percentage difference of less than tol_res dont count as resolved.
    #     rt=np.nan;pa=np.nan
    if pa_wt<0.75:
        rt=np.nan;pa=np.nan
    else:
        response_front_delta=abs(front_ss-response_front)
        response_back_delta=abs(response_back-back_ss)
        
        ss_indx1=np.argwhere((response_front_delta<tol_time)==True)
        ss_indx2=np.argwhere((response_back_delta<tol_time)==True)
        front_t= ss_indx1[0][0]
        back_t= ss_indx2[0][0]   
        rt=np.max([front_t,back_t]) 
        
        pa=pa_wt # if the to_res criteria is satisfied by the normalized series, use the actual polarization amplification
        
    
    return drg,pa,rt


def polarization_amplification_resolving(response,stimulus,rt):
    
    UL=response[:,-1][0]
    UR=response[:,-1][-1]
    
    if UL>UR:
        pa=np.nan
    else:
        pa=UR/UL
    # pa = (abs(UR-UL))*100
    
    SL=stimulus[:,-1][0]
    SR=stimulus[:,-1][-1]
    drg=SR/SL
    
    return drg,pa

################################################################# plotting
def plot_kymo(u,label,max_indx=None):

    tF,N=np.shape(u.T)    
    
    fig,ax = plt.subplots()
    im = ax.imshow(u.T, extent=[cellmem[0],cellmem[-1],len(u.T),0],cmap=parula_map, aspect = 'auto',vmin=np.min(u),vmax=np.max(u))
    ax.figure.colorbar(im)
    ax.set_ylabel('time(sec)',fontsize=20)
    ax.set_xlabel('plama Membrane bins',fontsize=20)
    # ax.set_xticks(cellmem)
    ax.set_xticklabels(np.round([0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi],2),fontsize=20)
    ax.set_yticks(np.arange(0,len(u.T),50))
    ax.set_yticklabels(np.arange(0,tF,50),fontsize=20)
    if not max_indx is None:
        plt.axvline(x=cellmem[max_indx],color='r')
    #        ax.set_title('@EGFRt'+str(egfrt)+'_spatial profile')
    plt.show()

def plot_gradient_profile(stimulus,max_indx=None):
    
    tF,N=np.shape(stimulus.T) 
      
    plt.figure()
    plt.figure(figsize=(9,3))
    plt.plot(cellmem,data_stimulus[:,100],'--',color='g',lw=3.0)
    plt.plot(cellmem,data_stimulus[:,-1],'-',color='g',lw=3.0)
    plt.ylabel('signal',fontsize=20)
    plt.xlabel('cell membrane contour',fontsize=20)
    plt.xticks(ticks=np.round([0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi],2), labels=np.round([0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi],2),fontsize=20)
    plt.xlim(0,cellmem[-1])
    # plt.ylim(0,5*threshold_stimu+0.05)
    plt.ylim(0,0.048)
    if not max_indx is None:
        plt.axvline(x=cellmem[max_indx],color='r')
    plt.show()

def plot_timeseries(response_front,response_back,polarization_time=None,title=None):
    
    ylim = np.max(response_front)+0.5
    plt.figure()
    plt.figure(figsize=(9,3))
    plt.plot(time,response_front,'k-',lw=3.0)
    plt.plot(time,response_back,'k--',lw=3.0)
    if polarization_time is not None:
        plt.axvline(x=polarization_time,color='r')
    plt.ylabel('response',fontsize=20)
    plt.xlabel('time(sec)',fontsize=20)
    if title is not None:
        plt.title('polarization type '+str(title))
    plt.xlim(0,time[-1])
    # plt.ylim(0.49,0.52)
    plt.show()

def plot_final(xdata,ydata_mean,ydata_std,labels,limits):
    xlabel,ylabel=labels
    xlim,ylim=limits
    plt.figure()
    plt.plot(xdata,ydata_mean,color='k',lw=2)
    plt.fill_between(xdata,ydata_mean-ydata_std,ydata_mean+ydata_std,color='k',alpha=0.3,edgecolor=None)
    plt.xlim(1,xlim)
    plt.xlabel(xlabel,fontsize=20,fontname='Arial')
    plt.ylim(0,ylim)
    plt.ylabel(ylabel,fontsize=20,fontname='Arial')
    plt.show()


###############################


# model = 'Legi';t_rev=800;threshold_stimu=0.005
# model = 'Otsuji';t_rev=800;threshold_stimu=0.001
# model = 'WavePinning';t_rev=300;threshold_stimu=0.003
model = 'SubPB';t_rev=300;threshold_stimu=0.012


experiment='single_gradient'
# experiment='gradient_reversal'
# experiment='simultaneous_gradients'
# experiment='single_gradient_with_offset'


folder_load=os.getcwd()+'\\'+model+'\\'+experiment+'\\'
folder_save_output = os.getcwd()+'\\final_output\\'

plot_single=None

if plot_single==True:
    Stimulus_strengths=[0.012]
    Offset=[2*threshold_stimu]
    # Offset=[0.0]
    # Stimulus_strengths=np.linspace(threshold_stimu,2*threshold_stimu,10)
    Trials=[2]#np.arange(0,10) 
    plot_data=True
    save_data=None
else:
    if experiment =='single_gradient':
        Stimulus_strengths=np.arange(0.0,0.05+0.001,0.001) # single gradient - Otsuji, Wavepinning, SubPB, Legi
        Offset=[0]    
    elif experiment =='gradient_reversal':
        Stimulus_strengths=np.round(np.arange(threshold_stimu,2*threshold_stimu+0.001,0.001),3) # gradient reversal all models
        Offset=[0]    
    elif experiment =='simultaneous_gradients':
        Stimulus_strengths = np.round(np.arange(threshold_stimu,2*threshold_stimu+0.001,0.001),3)# simulataneous gradients Wavepinning, SubPB and Legi
        Offset=[0]      
    elif experiment =='single_gradient_with_offset':
        Stimulus_strengths=[threshold_stimu]
        Offset=np.round(np.arange(0,5*threshold_stimu+0.001,0.001),3)
    
    Trials=np.arange(0,10,1)
    plot_data=None
    save_data=True

x_data=np.zeros((len(Trials),len(Offset)))*np.nan
y1_data=np.zeros((len(Trials),len(Offset)))*np.nan
y2_data=np.zeros((len(Trials),len(Offset)))*np.nan
z_data=np.zeros((len(Trials),len(Offset)))*np.nan

for stimulus_strength in Stimulus_strengths:
    
    stimulus_strength=np.round(stimulus_strength,3)
    print('Stimulus strength '+str(stimulus_strength))
    
    for i, offset in enumerate(Offset):
        
        offset=np.round(offset,3)
        print('offset '+str(offset))
    
        x_temp=np.zeros(len(Trials))*np.nan
        y1_temp = np.zeros(len(Trials))*np.nan
        y2_temp = np.zeros(len(Trials))*np.nan
        pol_type=np.zeros(len(Trials))*np.nan
        
        for j in range(len(Trials)):
            
            trial_num=Trials[j]  
            
            if experiment !='single_gradient_with_offset':            
                if trial_num<5 or model=='Legi':
                    filename=model+'_'+experiment+'_'+str(stimulus_strength)+'_'+str(trial_num)
                else:
                    filename=model+'_'+experiment+'_'+str(stimulus_strength)+'_'+str(offset)+'_'+str(trial_num)
            else:
                filename=model+'_'+experiment+'_'+str(stimulus_strength)+'_'+str(offset)+'_'+str(trial_num)
            
            cellmem=np.load(os.path.join(folder_load,'domain.npy'))
            time=np.load(os.path.join(folder_load,'time.npy'))
            data_response=np.load(os.path.join(folder_load,'response_'+filename+'.npy'))
            data_stimulus=np.load(os.path.join(folder_load,'stimulus_'+filename+'.npy'))
            
            data_response=data_response[:,:len(time)]
            data_stimulus=data_stimulus[:,:len(time)]
            
            if plot_data is not None:
                plot_gradient_profile(data_stimulus)
                plot_kymo(data_response,label='w')
                
     
            if experiment=='single_gradient':
                         
                response_front,response_back,alignment_type=time_series(data_response,data_stimulus,-1)
                
                tol=3e-2
                
                pt=polarization_time(response_front,response_back,tol=tol)
                # pd,pr=polarization_ratio(data_response,data_stimulus)
                pd,pol_max,pol_min,fc=polarization_ratio(time,response_front,response_back,data_stimulus,pt)
                pr=pol_max/pol_min
                
                if model=='WavePinning' or model=='SubPB':
                    if pr<2:
                        pt=np.nan
                
                x_temp[j]=pd
                y1_temp[j]=pr
                y2_temp[j]=pt
                
                if alignment_type==True:
                    title='True alignment'
                elif alignment_type==False:
                    title='False alignment'
                
                if plot_data is not None:
                    plot_timeseries(response_front,response_back,pt, title)
            
            elif experiment=='single_gradient_with_offset':
                         
                response_front,response_back,alignment_type=time_series(data_response,data_stimulus,-1)
                
                tol=3e-2
                
                pt=polarization_time(response_front,response_back,tol=tol)
                # pd,pr=polarization_ratio(data_response,data_stimulus)
                pd,pol_max,pol_min,fc=polarization_ratio(time,response_front,response_back,data_stimulus,pt)
                pr=pol_max/pol_min
                
                if model=='WavePinning' or model=='SubPB':
                    if pr<2:
                        pt=np.nan
                
                x_temp[j]=fc
                y1_temp[j]=pr
                y2_temp[j]=pt
                
                if alignment_type==True:
                    title='True alignment'
                elif alignment_type==False:
                    title='False alignment'
                
                if plot_data is not None:
                    plot_timeseries(response_front,response_back,pt, title)
            
            elif experiment=='gradient_reversal':
                           
                response_front,response_back,pol_type = time_series(data_response,data_stimulus,t_rev)
                
                
                response_front_before=response_front[:t_rev]
                response_back_before=response_back[:t_rev]
                time_before=time[:t_rev]
                data_stimulus_before=data_stimulus[:,:t_rev]
                
                pd,pol_max,pol_min=polarization_ratio(time_before,response_front_before,response_back_before,data_stimulus_before,pt=t_rev-2) # t_rev-2 to avoid empty slicing
                
                drg,pa=polarization_amplification_reversal(response_front,response_back,data_stimulus,treversal=t_rev)
    
                
                if pol_max/pol_min <2 and model!='Legi': # checking for initial polarization
                    title='No initial polarization'
                    rt=np.nan
                else:
                    title=None
                    
                    tol=3e-2
    
                    rt=reversal_time(response_front,response_back,tol,treversal=t_rev)
                
                if pa<=0.75: # filtering for right reversed polarization using polarization amplification.
                    rt=np.nan
                
                x_temp[j]=drg 
                y1_temp[j]=pa
                y2_temp[j]=rt-t_rev
                
                if plot_data is not None:
                    plot_timeseries(response_front,response_back,polarization_time=rt,title=title) 
                
            elif experiment=='simultaneous_gradients':  
                response_front,response_back,align_type = time_series(data_response,data_stimulus,-1)
                
                
                tol_time=3e-2
                if model=='Legi':
                    tol_time=0.5e-3 # for tol_time=0.3e-3 the avarage time was 607 sec. Now it is 506
                
                drg,pa,rst=resolving_quantities(response_front,response_back,data_stimulus,threshold_stimu,tol_time)
    
    
                x_temp[j]=drg 
                y1_temp[j]=pa
                y2_temp[j]=rst
                
                if plot_data is not None:
                    plot_timeseries(response_front,response_back,polarization_time=rst,title=None)
                
        x_data[:,i]=x_temp
        y1_data[:,i]=y1_temp
        y2_data[:,i]=y2_temp
        # Z[:,i]=pol_type
        

x_mean=np.nanmean(x_data,axis=0)    

y1_mean=np.nanmean(y1_data,axis=0)    
y1_std=np.nanstd(y1_data,axis=0)    
y2_mean=np.nanmean(y2_data,axis=0)
y2_std=np.nanstd(y2_data,axis=0)
            

if plot_single==None and plot_data is True:
    
    # plot_final(x_mean,y1_mean,y1_std,labels=['%','PR'],limits=[10,10])
    # plot_final(x_mean,y2_mean,y2_std,labels=['%','time'],limits=[10,1000])
    
    plot_final(x_mean,y1_mean,y1_std,labels=['stimulus ratio','PA'],limits=[5,10])
    plot_final(x_mean,y2_mean,y2_std,labels=['stimulus ratio','time'],limits=[5,1000])

if plot_single==None and save_data==True:
    output=np.zeros((5,len(Offset)))
    output[0]=x_mean
    output[1]=y1_mean
    output[2]=y1_std  
    output[3]=y2_mean
    output[4]=y2_std
    
    np.save(os.path.join(folder_save_output,model+'_'+experiment+'.npy'),output)
            
