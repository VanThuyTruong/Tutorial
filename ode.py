# -*- coding: utf-8 -*-
"""
Van Thuy Truong 

Python 3.8 : Mon Apr 06 17:41:36 2020
"""
from scipy import integrate
import numpy as np
import random
from scipy import integrate, interpolate
from numpy.random import normal

class Ode(object):
    """
    Set up ODE simulation
    """
    def __init__(self, param, samples=None):        
        self.ncell = param.ncell
        self.vec0 = np.array([param.Dose, param.Cp_0, param.Ct_0])
        
        if param.dose_schedule =='Single Dose':
            self.getSolution_single_dose(param)
        else:
            self.getSolution_multiple_dose(param)
            
        self.drugeffect = (param.Imax*self.tumour**param.h)/\
                            ((param.IC_50**param.h)+(self.tumour**param.h))
                            
        if samples is None:
            self.samples = np.array([param.sample])
        else:
            self.samples = samples
            
        self.sample_size = self.samples.shape[0]
        

        
        random.seed(9001)   #this one lead to one random distribution
        if param.vary_pERK0_uniform:   #pERK0 from uniform distribution
            print('if 1')
            self.ncells = np.empty((param.ncell,len(self.times)))
            self.pERKt = np.empty((param.ncell,len(self.times)))
            for i in range(param.ncell):
                param.pERK0 = np.array([random.uniform(param.minbaseline, param.maxbaseline)])              
                self.pERKt[i]=(100-self.drugeffect)*param.pERK0/100
                self.pERKt_function = interpolate.CubicSpline(self.times, self.pERKt[i])
                self.ncells[i] = self.ode_run_simulation(param).T #Total number of tumour
            self.num_cell_average=np.sum(self.ncells,axis=0)/param.ncell
            self.pERKtaverage=np.sum(self.pERKt,axis=0)/param.ncell 
        
        elif param.vary_pERK0_bimodal: #pERK0 from bimodal distribution
            print('if 2')
            self.ncells = np.empty((param.ncell,len(self.times)))
            self.pERKt = np.empty((param.ncell,len(self.times)))
            
            random.seed(9001)   
            self.ncells1 = np.empty((param.TN01, len(self.times)))
            self.pERKt1 = np.empty((param.TN01, len(self.times)))

            for i in range(param.TN01):
                param.TN0=param.TN01
                #mean1, sd1 = np.log(10), np.log(np.e)
                mean1, sd1 = 60, 10
                param.pERK01=random.gauss(mean1, sd1)
                while param.pERK01 > 200 or param.pERK01 < 0:
                    param.pERK01=random.gauss(mean1, sd1)
                self.pERKt1[i]=(100-self.drugeffect)*param.pERK01/100
                self.pERKt_function = interpolate.CubicSpline(self.times, self.pERKt1[i])
                self.ncells1[i] = self.ode_run_simulation(param).T #Total number of tumour

            random.seed(9001)   #this one lead to one random distribution
            self.ncells2 = np.empty((param.TN02, len(self.times)))
            self.pERKt2 = np.empty((param.TN02, len(self.times)))
    
            for i in range(param.TN02):
                param.TN0=param.TN02
                #mean2, sd2 = np.log(190), np.log(np.e)
                mean2, sd2 = 190, 10                
                param.pERK02=random.gauss(mean2, sd2)
                while param.pERK02 > 200 or param.pERK02 < 0:
                    param.pERK02=random.gauss(mean2, sd2)
                #print(param.pERK02)

                self.pERKt2[i]=(100-self.drugeffect)*param.pERK02/100
                self.pERKt_function = interpolate.CubicSpline(self.times, self.pERKt2[i])
                self.ncells2[i] = self.ode_run_simulation(param).T #Total number of tumour


            self.ncells=np.concatenate((self.ncells1,self.ncells2))
            print(self.ncells[:,0], self.ncells.shape)
            self.pERKt=np.concatenate((self.pERKt1,self.pERKt2))
            self.num_cell_average=np.sum(self.ncells,axis=0)/param.ncell
            self.pERKtaverage=np.sum(self.pERKt,axis=0)/param.ncell 

        elif param.vary_lambda: #different lambdas
            print('if vary')
            self.ncells = np.empty((param.no_lamda,len(self.times)))
            self.pERKt = np.empty((param.no_lamda,len(self.times)))
            random.seed(9001)   #this one lead to one random distribution
            for i in range(param.no_lamda):
                mean, sd = np.log(0.1), np.log(np.e)
                param.lamda=random.lognormvariate(mean, sd)
                print(round(param.lamda, 5))
            
                self.pERKt[i]=(100-self.drugeffect)*param.pERK0/100
                self.pERKt_function = interpolate.CubicSpline(self.times, self.pERKt[i])
                self.ncells[i] = self.ode_run_simulation(param).T #Total number of tumour
            self.num_cell_average=np.sum(self.ncells,axis=0)/param.no_lamda
            self.pERKtaverage=np.sum(self.pERKt,axis=0)/param.no_lamda 

            
        else:    #only one value for pERK0
            print('else')
            self.pERKt=(100-self.drugeffect)*param.pERK0/100
            self.pERKt_function = interpolate.CubicSpline(self.times, self.pERKt)
            self.ncells = self.ode_run_simulation(param).T[0] #Total number of tumour


    def ode_run_simulation(self, param): #run ODE simulation
            return integrate.odeint(tumour_growth, param.TN0, self.times,
                                    args=(param, self))

    def getSolution_single_dose(self, param): #single dose treatment
        self.times = np.linspace(0, param.tmax, param.Nt+1)
        vec = integrate.odeint(diff_eqns, self.vec0, self.times, args=(param,))
        self.tumour, self.central =  vec.T[2], vec.T[1]
    
    def getSolution_multiple_dose(self, param): #multiple dose treatment
        if param.tmax < param.num_repeats*param.repetition_fixed:
            raise TypeError("Increase param.tmax")

        output = np.zeros((3, int(1+param.tmax/param.step)))  #empty array for storing output for 3 variable
        self.times = np.zeros(int(1+param.tmax/param.step))
        #from 0 to num repeats*repetition fixed
        init = self.vec0
        for i in range(param.num_repeats):
            nt = int(1+param.repetition_fixed/param.step)
            times = np.linspace(param.repetition_fixed*i,
                                param.repetition_fixed*(i+1), nt+1)
            vec = integrate.odeint(diff_eqns, init, times, args=(param,))
            output[:,i*nt:(i+1)*nt] = vec[:-1].T
            self.times[i*nt:(i+1)*nt] = times[:-1]
            #update next repeat
            init = self.vec0 + vec[-1]
        nt = int(1+(param.tmax-param.num_repeats*\
                param.repetition_fixed)/param.step)
        
        nt2 = int(param.num_repeats*param.repetition_fixed/param.step)
        times = np.linspace(param.num_repeats*param.repetition_fixed,
                                 param.tmax, nt)
        vec = integrate.odeint(diff_eqns, init, times, args=(param,)).T
        output[:, nt2:] = vec
        self.times[nt2:] = times
        self.X0, self.central, self.tumour = output
        
def tumour_growth(ncells, t, param, Ode): #ODE for tumour population growth
    return (param.lamda*\
           ncells-((param.max_pERK-Ode.pERKt_function(t))*param.mu*ncells)/\
            param.max_pERK)


def diff_eqns(Vec, t, param): #ODE for PK
    return np.array([-param.ka*Vec[0],
                     param.ka*Vec[0]+param.Cltp*Vec[2]-param.Clpt*Vec[1]-\
                         param.Cl*Vec[1],
                     param.Clpt*Vec[1]-param.Cltp*Vec[2]])

