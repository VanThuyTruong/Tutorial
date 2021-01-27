#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Van Thuy Truong : PhD Researcher in Clinical Pharmacology Biologics
                  and Bioanalysis : van.thuytruong@astrazeneca.com
Aaron Klug Building, Granta Park, Cambridge, CB21 6GH

Python 3.8 : Mon Apr 06 17:41:36 2020
"""
from graphic.Plots import plot_setup

def ode_total_number_cells(param,ODE):
    """
    Plots the total number of cells in multiple populations against time.

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ODE : ODE object
        Ordinary differential equation simulation.
    """
    fig, ax = plot_setup('Time (hour)', 'Number of tumour cells', 
                                 scale=param.plot_scale)

    ax.set(xlim=(0,param.tmax))
    ax.grid(linewidth=0.5)
    label1= 'Number of tumour cells per pERK scenario'
    if param.vary_pERK0_uniform or param.vary_pERK0_bimodal:
        for i in range(param.ncell):
            ax.plot(ODE.times,ODE.ncells[i],color='blue',label=label1)
            label1= "_nolegend_"
            ax.set_ylim([0, max(ODE.ncells[i])])
            #ax.set_ylim([0, 200])
            #ax.set_ylim(ymin=0)
        print(ODE.ncells.shape, ODE.ncells[:,0])
        ax.plot(ODE.times, ODE.num_cell_average, color='yellow',
                label='Average number of tumour cells')
    elif param.vary_lambda:
        for i in range(param.no_lamda):
            ax.plot(ODE.times,ODE.ncells[i],color='blue',label=label1)
            label1= "_nolegend_"
            #ax.set_ylim(ymin=0)
        ax.plot(ODE.times, ODE.num_cell_average, color='yellow',
                label='Average number of tumour cells')
    else:
        ax.plot(ODE.times,ODE.ncells,color='blue',label=label1)
        #ax.set_ylim(ymin=0)
    
    ax.legend(loc='upper center', bbox_to_anchor=(0.8, -0.3),
              fontsize=16*param.plot_scale)

    
    if param.save_plot:
        from graphic.Plots import save_plot
        save_plot(fig, ax, file_name='total number of cells',
              folder_name=param.file_dir)
    else:
        import matplotlib.pyplot as pt
        ax.legend(loc='upper center', bbox_to_anchor=(0.7, -0.4),fontsize=30*param.plot_scale)
        pt.show()

def ode_individual_pERK(param, ODE):
    """
    Plots the pERK value inside the populations against time.

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ODE : ODE object
        Ordinary differential equation simulation.
    """
    fig, ax = plot_setup('Time (hour)', 'pERK value inside the population', 
                                 scale=param.plot_scale)
    label1= 'pERK profile of one population'
    if param.vary_pERK0_uniform or param.vary_pERK0_bimodal:
        for i in range(param.ncell):
            ax.plot(ODE.times,ODE.pERKt[i],color='purple', label=label1)
            label1= "_nolegend_"
        ax.plot(ODE.times, ODE.pERKtaverage,color='yellow',
                label='Average pERK value inside the populations')
        ax.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15),fontsize=30*param.plot_scale)
    elif param.vary_lambda:
        for i in range(param.no_lamda):
            ax.plot(ODE.times,ODE.pERKt[i],color='purple', label=label1)
            label1= "_nolegend_"
            
    else:
        ax.plot(ODE.times, ODE.pERKt, color='purple', label=label1)
    ax.legend(loc='upper center', bbox_to_anchor=(0.8, -0.3),
              fontsize=16*param.plot_scale)
    
    if param.save_plot:
        ax.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15),fontsize=30*param.plot_scale)
        from graphic.Plots import save_plot
        save_plot(fig, ax, file_name='Individual pERK',
              folder_name=param.file_dir)
        
    else:
        import matplotlib.pyplot as pt
        ax.legend(loc='upper center', bbox_to_anchor=(0.7, -0.15),fontsize=30*param.plot_scale)
        pt.show()
    
def ode_pkpd_plot(param, ODE): 
    """
    Plots the PKPD against time.

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ODE : ODE object
        Ordinary differential equation simulation.
    """
    fig, ax = plot_setup('Time (hour)', 'Amount (μmol/kg)', 
                                 scale=param.plot_scale)
    ax2 = ax.twinx()
    ax.plot(ODE.times, ODE.central, color='red',
            label='Drug amount in the plasma')
    ax.plot(ODE.times, ODE.tumour, color='green',
            label='Drug amount in the tumour tissue')
    ax2.plot(ODE.times, ODE.drugeffect, color='orange',
             label='percentage PERK decrease')
    ax2.set(ylim=(0, 100))
    ax2.set_ylabel('Percentage Decrease pERK', fontsize=param.plot_scale*20)
    ax2.tick_params(axis="y", labelsize=16*param.plot_scale)
    ax.legend(loc='upper center', fontsize=16*param.plot_scale)
    ax2.legend(loc='upper right', fontsize=16*param.plot_scale)
      
    if param.save_plot:
        from graphic.Plots import save_plot
        save_plot(fig, ax, file_name='PKPD',
              folder_name=param.file_dir)
    else:
        import matplotlib.pyplot as pt
        pt.show()
