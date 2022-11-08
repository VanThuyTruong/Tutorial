#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Van Thuy Truong 

Python 3.8 : Mon Apr 06 17:41:36 2020
"""

import os, random, namelist
import numpy as np

param = namelist.Namelist()
###plots###
param.plot_scale = 1.0 #scale of plots w.r.t screen
param.cmap = 'seismic' #Colour map for plots
param.save_plot = True  #Save plot figures?
param.print_info = False #Print Information relating to system parameters?
param.plot_pool = True #plot death, quiescence, division pool
param.pkpd_plot = True #plot pkpd profile
param.total_number_cells = True  #plot total cell number in one population
param.plot_pERK = True  #plot pERK scatterplot
param.create_video = True #creae video with pERK plots
param.create_video_all = False #create video with 4 plots (pERK, total cell number, pool, pkpd)
param.plot_all = False   #plot 4 plots with pERK, total cell number, pool, pkpd
param.pERK_t_plot = True  #plot pERK of each cell over time

param.samples_plot = False  #multiple samples
param.sample_size = 100 #how many different samples?
#pERK values if multiple samples are chosen
param.multiple_unimodal_samples = False #samples where the pERK inside the population is similar 
param.multiple_unimodal = False #unimodal sample with pERK values of 100%
param.multiple_bimodal = False #bimodal sample with pERK values of 60% and 190%
param.multiple_uniform = True #uniformly distributed pERK

#####drug and dose#####
param.M = 531.3 #molecular weight of the drug (g/mol)
param.Cp_0, param.Ct_0 = 0, 0 #initial plasma/tumour amount

# param.X1, param.X3, param.X10 = 1, 3, 10 #Dose (mg/kg)
# param.D1, param.D3, param.D10 = (param.X1/param.M)*1e3, \
#         (param.X3/param.M)*1e3, (param.X10/param.M)*1e3 #Dose (micromol/kg)
param.X = 3 #1,10 #Dose (mg/kg)
param.Dose = (param.X/param.M)*1e3 #dose of the drug #Dose (micromol/kg)
param.dose_schedule='Single Dose' #'Multiple Dose' or 'Single Dose'

#multiple dose
param.step = 0.01   #number of timesteps of the multiple dose treatment
param.num_repeats = 10 #number of repeats of the multiple dose treatment
param.repetition_fixed = 24 #duration of a cycle for the multiple dose treatment

#pk
param.ka = 1.08 #absorption constant h**-1  
param.Cl = 3.13 #elimination constant(L/h/kg)  
param.Clpt = 1.74 #intercompartmental distribution rate from plasma to tumor compartment (L/h/kg) 
param.Cltp = 4e-2 #intercompartmental distribution rate from tumor to plasma compartment(L/h/kg)
#param.VF = 4.05e1 #apparent distribution volume (L/kg)

#pd
param.Imax  = 97# Percentage of maximum value
param.IC_50 = .78 #Drug concentration in tumour compartment which causes pERK decrease half of Imax (micromol/L)
param.h = 1 #Hill coefficient

######abm######
param.μ, param.λ =0.0828/24, 0.0828/24 #hours, abm death and birth rate
param.ncell = 100  #initial number of tumor cells
param.minbaseline, param.maxbaseline = 0, 2e2 #min, max baseline of pERK value inside cell

#division and death thresholds
param.trep, param.tdeath = 0.5*param.maxbaseline, 0.25*param.maxbaseline

#random distribution of pERK values
param.single_pop_uniform = True #uniformly distributed pERK
param.single_pop_bimodal =False #bimodaly distributed pERK

if param.single_pop_bimodal:
    from numpy.random import normal
    np.random.seed(8900)
    param.sample1 = normal(loc=60, scale=10, size=int(param.ncell*0.5))
    param.sample2 = normal(loc=190, scale=10, size=int(param.ncell*0.5))
    param.tsample = np.hstack((param.sample1, param.sample2))
    param.sample = [param.tsample[i] for i in range(param.ncell)]

if param.single_pop_uniform:
    random.seed(8900)
    param.sample = np.array([random.uniform(param.minbaseline, param.maxbaseline) \
                    for i in range(param.ncell)])


######ODE model#######
param.max_pERK = 200 #max pERK in the ODE model
param.lamda=0.0828/24 #division rate
param.mu=0.16/24 #death rate
param.TN0=100 #total number of cells in a population
param.pERK0=100 #initial pERK value

#ode vary pERK0 (initial pERK value)
param.vary_pERK0_uniform = True #uniform distribution of pERK0
param.vary_lambda = False  #vary lamda?
param.no_lamda=10   #number of lamda variation

#bimodal pERK0
param.vary_pERK0_bimodal = False  #bimodal distribution of pERK0
param.TN01=50 #number of cells in one distribution
param.TN02=50 #number of cells in the other distribution

#ode plots
param.ode_plot_pkpd = True  #plot pkpd against time
param.ode_individual_pERK = True #plot pERK values against time
param.ode_total_cell_number = True #plot number of cells in a population 


param.sample_distr = "Uniform"   #name folder for saving the plots
param.file_dir = os.path.join(param.sample_distr, 'Simulations 2',
                              '3mgkg single dose') #name subfolders

##time of the simulation##
param.tmax, param.Nt = 250, int(1e4) #Final time and number of time steps
param.fps = 5 #40 faster if larger number  frames/second (how many frames to show per second)
param.dt = 1#0.1 #param.tmax/param.Nt #time step for quiescence


param.model = 'ABM'   #run 'ABM' or 'ODE' model?

def setup():
    param.checkall()
    return param
