# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 08:42:58 2023

@author: Akhilesh Nandan
"""
import numpy as np

class stimulus_single:
    
    '''
    this class generates a single pulse of stimulus
    between two specified time points, stimulus_beg and stimulus_end.
    '''
    
    stimulus_beg=10
    stimulus_end=70
    
    def add_stimulus(self,t,input_params):
        
        
        if self.stimulus_beg<t<self.stimulus_end:
            sL=input_params[1]
            sR=input_params[2]
        else:
            sL=np.min(input_params[1:])
            sR=np.min(input_params[1:])
        return sL,sR
    
class stimulus_ramp:
    
    '''
    this class generates step like ramping pulse of stimuli
    between two specified time points, stimulus_beg and stimulus_end.
    
    the duration of steps, step_size is determined by segmenting the 
    total duration.
    '''
    
    stimulus_beg=10
    stimulus_end=70
    
    def add_stimulus(self,t,input_params):

        step_size=(self.stimulus_end-self.stimulus_beg)/5
        
        if self.stimulus_beg<=t<self.stimulus_beg+step_size:
            sL=input_params[1]/4
            sR=input_params[2]
        elif self.stimulus_beg+step_size<=t<self.stimulus_beg+2*step_size:
            sL=input_params[1]/2
            sR=input_params[2]
        elif self.stimulus_beg+2*step_size<=t<self.stimulus_beg+3*step_size:
            sL=input_params[1]
            sR=input_params[2]
        elif self.stimulus_beg+3*step_size<=t<self.stimulus_beg+4*step_size:
            sL=input_params[1]/2
            sR=input_params[2]
        elif self.stimulus_beg+4*step_size<=t<self.stimulus_beg+5*step_size:
            sL=input_params[1]/4
            sR=input_params[2]
        else:
            sL=np.min(input_params[1:])
            sR=np.min(input_params[1:])
        return sL,sR
    
