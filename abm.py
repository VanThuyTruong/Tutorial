#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Van Thuy Truong

Python 3.8 : Mon Apr 06 17:41:36 2020
"""

from scipy import integrate
import numpy as np

class Abm(object):
    """
    Set up ABM simulation
    """
    def __init__(self, param, samples=None):
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
        self.birth_count = np.zeros(self.sample_size)
        self.death_count = np.zeros(self.sample_size)

    def run_simulation(self, param):
        """
        Set up Gillespie algorithm
        """
        import random, math
        from cells import Cells
        from scipy.interpolate import CubicSpline

        #information to pass forward when plotting
        self.ind_max = 0
        self.death_pools = np.empty((self.sample_size, self.ind_max))#array with shape 1 or number of samples, max timestep
        self.division_pools = np.empty((self.sample_size, self.ind_max))
        self.quiesence_pools = np.empty((self.sample_size, self.ind_max))
        self.inds = np.empty(self.sample_size)
        self.tts = np.empty((self.sample_size, self.ind_max))
        self.tumour_tts = np.empty((self.sample_size, self.ind_max))
        self.central_tts = np.empty((self.sample_size, self.ind_max))
        self.drugeffect_tts= np.empty((self.sample_size, self.ind_max))
        self.pERK_vals = np.empty((self.sample_size, self.ind_max))
        self.ncells = np.empty((self.sample_size, self.ind_max), dtype=int)
        self.image_step = np.empty(self.sample_size, dtype=int)
        total_ncells = 0

        #Loopthrough samples
        for i, sample in enumerate(self.samples):
            param.sample = sample
            cells, ind = Cells(param), 0
            tt, pERK_vals, ncells = np.zeros(1), \
                cells.pERK, np.array(cells.ncell, dtype=int)
            
            death_pool, division_pool, quiesence_pool = \
                np.array([cells.n_death]),  np.array([cells.n_rep]), \
                np.array([cells.n_quiesence])
                
            #ABM Simulation with Gillespie algorithm
            alphan_arr = np.empty(0)
            death_arr = np.empty(0)
            while tt[-1] < param.tmax:
                t_string = 'Time={:.4f} h'.format(tt[-1])
                if param.print_info:
                    print('\t', t_string)
                    
                assert len(np.where(cells.pERK > param.trep)[0]) == \
                    cells.n_rep
                assert len(np.where(cells.pERK < param.tdeath)[0]) == \
                    cells.n_death
                    
                
                if cells.n_rep==cells.n_death==0:
                    assert len(np.where(cells.pERK < param.tdeath)[0]) == 0
                    assert len(np.where(cells.pERK > param.trep)[0]) == 0
                    dt = param.dt #all cells are in the quiescence pool, 
                                  #continue with timesteps                
                else:
                    birthrate = param.λ*0.01*cells.division_sumpERK
                    deathrate = param.μ*0.01*cells.death_sumpERK
                    alphan = birthrate + deathrate
                    alphan_arr = np.append(alphan_arr, alphan)
                    death_arr = np.append(death_arr, deathrate)
        
                    if random.random() < deathrate/alphan:
                        cells.death(param)
                        self.death_count[i]+=1
                    else:
                        cells.division(param, tt[-1])
                        self.birth_count[i]+=1
                    n = random.random()
                    a = -math.log(n)
                    dt = a/alphan

                newtime = tt[-1] + dt
                index = np.searchsorted(self.times, newtime)-1
                pERK_decrease = self.drugeffect[index]

        
                #connect abm with pd
                cells.update_pERK(pERK_decrease, newtime)
                pERK_vals = np.append(pERK_vals, cells.pERK)
                ncells = np.append(ncells, cells.ncell)
                cells.update_pool(param)
                death_pool = np.append(death_pool, cells.n_death)
                division_pool = np.append(division_pool, cells.n_rep)
                quiesence_pool = np.append(quiesence_pool, cells.n_quiesence)
                tt = np.append(tt, newtime)
                ind += 1
                
                if cells.ncell == 0:
                    break

            ncells_sum = np.sum(ncells)  #####
            if ncells_sum > total_ncells:
                empty = np.empty((self.sample_size, ncells_sum-\
                                  total_ncells))
                empty[:] = np.NaN
                self.pERK_vals = np.concatenate((self.pERK_vals, empty),
                                                axis=1)
                total_ncells = ncells_sum

            #Append simulation data from samples
            if ind >= self.ind_max:
                #extend samples data if necessary
                self.extend(ind+1)
            self.image_step[i] = int(np.log10(ind))
            self.pERK_vals[i, :ncells_sum] = pERK_vals
            self.central_tts[i, :ind+1] = CubicSpline(self.times,
                            self.central)(tt)
            self.tumour_tts[i, :ind+1] = CubicSpline(self.times, 
                           self.tumour)(tt)

            self.drugeffect_tts[i, :ind+1] = CubicSpline(self.times,
                            self.drugeffect)(tt)
            self.tts[i, :ind+1] = tt
            self.death_pools[i, :ind+1] = death_pool
            self.division_pools[i, :ind+1] = division_pool
            self.quiesence_pools[i, :ind+1] = quiesence_pool
            self.inds[i] = ind
            self.ncells[i, :ind+1] = ncells
            self.cells = cells

           
    def extend(self, ind):  ###
        empty = np.empty((self.sample_size, ind-self.ind_max))
        empty[:] = np.NaN
        self.ind_max = ind
        for att in ["tumour_tts", "central_tts","drugeffect_tts", "tts", 
                    "death_pools","division_pools", "quiesence_pools", 
                    "ncells","pERK_vals"]: #list of arrays you want to update
            arr = getattr(self, att)
            setattr(self, att, np.concatenate((arr, empty), axis=1))

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


def diff_eqns(Vec, t, param): #ODE for PK
    return np.array([-param.ka*Vec[0],
                     param.ka*Vec[0]+param.Cltp*Vec[2]-param.Clpt*Vec[1]-\
                         param.Cl*Vec[1],
                     param.Clpt*Vec[1]-param.Cltp*Vec[2]])
