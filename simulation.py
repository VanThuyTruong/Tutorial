#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Van Thuy Truong 

Python 3.8 : Mon Apr 06 17:41:36 2020
"""

import numpy as np
import os
import sys
np.set_printoptions(threshold=sys.maxsize)

def main(param, samples=None):
    """
    Run ABM with a single sample or multiple samples or ODE model,
    choose which plots are created

    Parameters
    ----------
    param : Namelist object
        System parameters.
    """

    if param.save_plot:
        from PPP.File_Management import dir_assurer
        dir_assurer(param.file_dir)
    
    if param.model == 'ABM':
        from abm import Abm
        from abm_plots import plot_pERK, plot_pool, pkpd_plot, \
                total_number_cells, plot_all, samples_plot

        if samples is None:
            ABM = Abm(param)
            
        else:
            ABM = Abm(param, samples)
            
        ABM.run_simulation(param)
        
        cells = ABM.cells
        if param.pERK_t_plot:
            cells.plot_pERK(param,ABM)
        if param.plot_pERK:
            plot_pERK(param, ABM)
        if param.create_video:
            create_video(param, ABM)
        if param.total_number_cells:
            total_number_cells(param, ABM)
        if param.plot_pool:
            plot_pool(param, ABM)
        if param.pkpd_plot:
            pkpd_plot(param, ABM)
        if param.plot_all:
            plot_all(param, ABM)
        if param.create_video_all:
            create_video_all(param, ABM)
        if param.samples_plot:
            samples_plot(param,ABM)
        
    else: #"Has to be ODE"
        from ode import Ode
        from ode_plots import ode_pkpd_plot, ode_individual_pERK, \
                ode_total_number_cells
        if samples is None:
            ODE = Ode(param)
            
        else:
            ODE = Ode(param, samples)
            
        ODE.ode_run_simulation(param)
        
        if param.ode_plot_pkpd:
            ode_pkpd_plot(param, ODE)
        if param.ode_individual_pERK:
            ode_individual_pERK(param, ODE)
        if param.ode_total_cell_number:
            ode_total_number_cells(param, ODE)

def create_video(param, ABM): 
    """
    Creating a video from all scatterplots containing pERK values of each cell

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """
    from PPP.Animate import animate_images

    for i in range(ABM.sample_size):
        animate_images(os.path.join(param.file_dir, 'Sample{}'.format(i+1),
                                    'Cell Numbers'),
                       int(ABM.inds[i]), param.fps, 
                       os.path.join(param.file_dir, 'Sample{}'.format(i+1))
                       )
def create_video_all(param, ABM): 
    """
    Creating a video from 4 plots with pERK, total cell number, different pools, 
    and pkpd against time

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """
    from PPP.Animate import animate_images

    for i in range(ABM.sample_size):
        animate_images(os.path.join(param.file_dir, 'Sample_all_video{}'.format(i+1),
                                    'All plots'),
                       int(ABM.inds[i]), param.fps, 
                       os.path.join(param.file_dir, 'Sample_all_video{}'.format(i+1))
                       )

if __name__=="__main__":
    import config
    from numpy.random import normal
    param = config.setup()

    
    #population with simular pERK from 10-200
    if param.multiple_unimodal_samples:
        sample_array = np.array([normal(loc=param.max_pERK*i/param.sample_size,
                            scale=5, size=param.ncell) for i in \
                            range(1, param.sample_size+1)])
        
    #only 100 pERK
    if param.multiple_unimodal:
        sample_array = np.full((10,100),100)
    
    # #bimodal sample with pERK 60 and 190 percent   
    if param.multiple_bimodal:
        sample_array=np.zeros((param.sample_size,param.ncell))
        for i in range(param.sample_size):
            sample1 = normal(loc=60, scale=10, size=int(param.ncell*0.5))
            
            sample2 = normal(loc=190, scale=10, size=int(param.ncell*0.5))
            tsample = np.hstack((sample1, sample2))
            sample_array[i] = np.vstack(tsample)[:, 0]
        
    #print(sample.shape, sample)

    #population with uniform pERK
    if param.multiple_uniform:
        sample_array =  np.random.uniform(param.minbaseline, param.maxbaseline, 
                                          size=(param.sample_size, param.ncell))


    
    if param.samples_plot:
        main(param, sample_array) #multiple samples
    else:
        main(param) #single sample