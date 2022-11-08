#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Van Thuy Truong 

Python 3.8 : Mon Apr 06 17:41:36 2020
"""
from PPP.Plots import plot_setup
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def plot_pool(param, ABM):
    """
    Plot number of cells in the division, death, and quiescence pool
    
    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """
    for i in range(ABM.sample_size):
        fig, ax= plot_setup('Time (hour)', 'Number of tumour cells', 
                                 scale=param.plot_scale)
        ax.plot(ABM.tts[i], ABM.death_pools[i], color='red',
                label='Death pool')
        ax.plot(ABM.tts[i], ABM.division_pools[i], color='green',
                label='Division pool')
        ax.plot(ABM.tts[i], ABM.quiesence_pools[i], color='black',
                label='Quiesence pool')
        ax.legend(loc='upper center', bbox_to_anchor=(0.88, -0.15),
                      fontsize=16*param.plot_scale)
        if param.save_plot:
            from PPP.Plots import save_plot
            save_plot(fig, ax, file_name='Pools_{}'.format(i+1),
                  folder_name=os.path.join(
                          param.file_dir, 'Sample%s'%(i+1)
                          )
                  )
        else:
            import matplotlib.pyplot as pt
            pt.show()
            

def plot_pERK(param, ABM):
    """
    Scatterplot of pERK values against individual cells

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """    
    for j in range(ABM.sample_size):
        ncell_sum = 0
        ncells = ABM.ncells[j][~np.isnan(ABM.ncells[j])] #removes NaN
        max_ncells = np.nanmax(ncells) #maximum ncell excluding NaN

        for i, ncell in enumerate(ncells):
            ncell= int(ncell)
            ncell_sum += ncell
            t_string = 'Time={:.4f} h'.format(ABM.tts[j, i])
            c_t='$X_t$={:.4f} μmol/kg' .format(ABM.tumour_tts[j, i])
            c_p='$X_p$={:.4f} μmol/kg' .format(ABM.central_tts[j, i])
            #print('\t', t_string)
            fig, ax = plot_setup('Individual cell', 'pERK value', 
                                 scale=param.plot_scale,
                                 title="%s, %s, %s" %(t_string,c_t,c_p))
            ax.scatter(range(ncell),
                       ABM.pERK_vals[j, ncell_sum-ncell:ncell_sum])
            
            ax.axhline(y=param.trep, color='g', linestyle='-',
                        label='Division threshold')
            ax.axhline(y=param.tdeath, color='r', linestyle='-',
                        label='Death threshold')
            ax.legend(fontsize=16*param.plot_scale,loc='upper right')
            ax.set_ylim([param.minbaseline, param.maxbaseline])
            ax.set_xlim([0, max_ncells])
            if param.save_plot:
                from PPP.Plots import save_plot
                save_plot(fig, ax, file_name=str(i),
                          folder_name=os.path.join(
                          param.file_dir, 'Sample%s'%(j+1)
                          )
                          )
            else:
                import matplotlib.pyplot as pt
                pt.show()
                
              
def plot_all(param, ABM):
    """
    Create 4 plots with pERK, total cell number, different pools, and 
    pkpd against time

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """    
    for j in range(ABM.sample_size):
        ncell_sum = 0
        ncells = ABM.ncells[j][~np.isnan(ABM.ncells[j])] #removes NaN
        max_ncells = np.nanmax(ncells) #maximum ncell excluding NaN
        drug=np.concatenate((ABM.tumour_tts[j], ABM.central_tts[j]),axis=0)
        max_drug=np.nanmax(drug)
        cells_pool=np.concatenate((ABM.death_pools[j],  ABM.division_pools[j], ABM.quiesence_pools[j]), axis=0)
        max_pool=np.nanmax(cells_pool)

        for i, ncell in enumerate(ncells):
            ncell= int(ncell)
            ncell_sum += ncell
            t_string = 'Time={:.3f} h'.format(ABM.tts[j, i])
            c_t='$X_t$={:.2f} μmol/kg' .format(ABM.tumour_tts[j, i])
            c_p='$X_p$={:.2f} μmol/kg' .format(ABM.central_tts[j, i])
            pERK_decrease='{:.2f} % pERK decrease' .format(ABM.drugeffect_tts[j, i])
            active_pERK ='{:.2f} % active pERK' .format(100-ABM.drugeffect_tts[j, i])
            ncells='{} cells' .format(int(ABM.ncells[j, i]))
            d_p='Death pool: {} cells' .format(int(ABM.death_pools[j, i]))
            b_p='Division pool: {} cells' .format(int(ABM.division_pools[j, i]))
            q_p='Quiesence pool: {} cells' .format(int(ABM.quiesence_pools[j, i]))
            

            fig2 = plt.figure(figsize=(12.8, 9.6))  #
            spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig2)
            spec2.update(wspace=0.5, hspace=0.5)
            
            
            f2_ax1 = fig2.add_subplot(spec2[1, 0])
            f2_ax1.grid()
            f2_ax1.set_title("%s, %s" %(t_string,ncells))
            f2_ax1.set_xlabel('Individual cell')
            f2_ax1.set_ylabel('pERK value')
            f2_ax1.scatter(range(ncell),
                       ABM.pERK_vals[j, ncell_sum-ncell:ncell_sum])
            
            f2_ax1.axhline(y=param.trep, color='g', linestyle='-',
                        label='Division threshold')
            f2_ax1.axhline(y=param.tdeath, color='r', linestyle='-',
                        label='Death threshold')
            f2_ax1.legend(loc='upper left')
            f2_ax1.set_ylim([param.minbaseline, param.maxbaseline])
            f2_ax1.set_xlim([0, max_ncells])
            
            f2_ax2= fig2.add_subplot(spec2[0, 0])
            f2_ax2.grid()
            ax2 = f2_ax2.twinx()
            f2_ax2.set_title("%s, %s, %s" %(c_p,c_t,pERK_decrease))
            f2_ax2.set(xlabel='Time (h)', ylabel='Amount (μmol/kg)')
            f2_ax2.plot(ABM.tts[j,:i+1], ABM.tumour_tts[j,:i+1], color='green',
                  label='Drug amount in the tumour tissue')
            f2_ax2.plot(ABM.tts[j,:i+1], ABM.central_tts[j,:i+1], color='red',
                        label='Drug amount in the plasma')
            f2_ax2.set_ylim([0, max_drug])
            f2_ax2.set_xlim([0, param.tmax])
            ax2.plot(ABM.tts[j,:i+1], ABM.drugeffect_tts[j,:i+1], color='orange',
                label='Percentage PERK decrease')
            ax2.set(ylim=(0, 100))
            ax2.set_ylabel('Percentage Decrease pERK',
                       fontsize=param.plot_scale*10)
            ax2.tick_params(axis="y", labelsize=10*param.plot_scale)
            lines_1, labels_1 = f2_ax2.get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            lines = lines_1 + lines_2
            labels = labels_1 + labels_2
            f2_ax2.legend(lines, labels,loc='upper left')


            f2_ax3 = fig2.add_subplot(spec2[1, 1])                       
            f2_ax3.grid()
            f2_ax3.set_title("%s, %s, %s" %(b_p,q_p,d_p))
            f2_ax3.set(xlabel='Time (h)', ylabel='Number of tumour cells')
            f2_ax3.plot(ABM.tts[j,:i+1], ABM.death_pools[j,:i+1], color='red',
                        label='Death pool')
            f2_ax3.plot(ABM.tts[j,:i+1], ABM.division_pools[j,:i+1], color='green',
                        label='Division pool')
            f2_ax3.plot(ABM.tts[j,:i+1], ABM.quiesence_pools[j,:i+1], color='black',
                        label='Quiesence pool')
            f2_ax3.set_xlim([0, param.tmax])
            f2_ax3.set_ylim(ymin=0)
            f2_ax3.legend(loc='upper left')
            
           
            f2_ax4 = fig2.add_subplot(spec2[0, 1])                       
            f2_ax4.grid()
            ax2 = f2_ax4.twinx()
            f2_ax4.set_title("%s, %s" %(ncells,active_pERK))
            f2_ax4.set(xlabel='Time (h)', ylabel='Number of tumour cells')
            f2_ax4.plot(ABM.tts[j,:i+1], ABM.ncells[j,:i+1], color='blue',
                         label='Number of tumour cells')
            f2_ax4.set_xlim([0, param.tmax])
            f2_ax4.set_ylim(ymin=0)
            ax2.set_ylabel('Percentage of active pERK')
            ax2.plot(ABM.tts[j,:i+1], 100-ABM.drugeffect_tts[j,:i+1], color='black',
                     label='Percentage of active pERK')
            ax2.set(ylim=(0, 100))
            lines_1, labels_1 = f2_ax4.get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            lines = lines_1 + lines_2
            labels = labels_1 + labels_2
            f2_ax4.legend(lines, labels,loc='upper left',)


            if param.save_plot:
                folder = os.path.join(
                        param.file_dir, 'Sample_all_video%s'%(j+1)
                        )
                file_name = os.path.join(folder, str(i)
                          )
                from PPP.File_Management import dir_assurer
                dir_assurer(folder)
                fig2.savefig(file_name)
                plt.close(fig2)
            else:
                import matplotlib.pyplot as pt
                pt.show()
                

def pkpd_plot(param, ABM):
    """
    Plot of the PKPD against time

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """
    for i in range(ABM.sample_size):
        fig, ax = plot_setup('Time (hour)', 'Amount (μmol/kg)', 
                                 scale=param.plot_scale)
        ax2 = ax.twinx()
        ax.plot(ABM.times, ABM.central, color='red',
                label='Drug amount in the plasma')
        ax.plot(ABM.times, ABM.tumour, color='green',
                label='Drug amount in the tumour tissue')
        ax2.plot(ABM.times, ABM.drugeffect, color='orange',
                 label='percentage PERK decrease')
        ax2.set(ylim=(0, 100))
        ax2.set_ylabel('Percentage Decrease pERK',
                       fontsize=param.plot_scale*20)
        ax2.tick_params(axis="y",labelsize=16*param.plot_scale)
        #fig.subplots_adjust(bottom=0.25) 
        
        lines_1, labels_1 = ax.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        lines = lines_1 + lines_2
        labels = labels_1 + labels_2
        ax.legend(lines, labels,bbox_to_anchor=(1.05, -0.12),
                  fontsize=16*param.plot_scale)
        
        if param.save_plot:
            from PPP.Plots import save_plot
            save_plot(fig, ax, file_name='PKPD_{}'.format(i+1),
                  folder_name=os.path.join(
                          param.file_dir, 'Sample{}'.format(i+1)
                          )
                  )
        else:
            import matplotlib.pyplot as pt
            pt.show()


def total_number_cells(param, ABM):
    """
    Plots the total number of cells against time.

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """
    for i in range(ABM.sample_size):
        fig, ax = plot_setup('Time (hour)', 'Number of tumour cells', 
                                 scale=param.plot_scale)
        ax.plot(ABM.tts[i], ABM.ncells[i], color='blue',
                label='Number of tumour cells in one realisation')
        ax.set_xlim([0, param.tmax])
        #ax.set_ylim([0,200])#ymin=0)
        ax.set_ylim(ymin=0)
        ax.legend(bbox_to_anchor=(1.0, -0.15),fontsize=16*param.plot_scale)


        if param.save_plot:
            from PPP.Plots import save_plot
            save_plot(fig, ax, file_name='Percentage pERK decrease',
                  folder_name=os.path.join(
                          param.file_dir, 'Sample{}'.format(i+1)
                          )
                  )
        else:
            import matplotlib.pyplot as pt
            pt.show()

def samples_plot(param,ABM):
    """
    Plots the total number of cells in multiple populations against time.

    Parameters
    ----------
    param : Namelist object
        System parameters.
    ABM : Abm object
        Agent-based model simulation.
    """
    fig, ax = plot_setup('Time (hour)', 'Number of tumour cells', 
                             scale=param.plot_scale)
    label1= 'Number of tumour cells in one realisation'
    ex_arr = np.empty(0)
    for i in range(ABM.sample_size):
        ax.plot(ABM.tts[i], ABM.ncells[i], color='blue', label=label1)
        label1="_nolegend_"
        ex_arr = np.append(ex_arr,ABM.tts[i][np.isfinite(ABM.tts[i])][-1])
    ax.set_xlim([0, param.tmax])
    #ax.set_ylim([0,200])
    ax.set_ylim(ymin=0)
    #ax.legend(bbox_to_anchor=(1.0, -0.15),fontsize=16*param.plot_scale)
    if param.save_plot:
        from PPP.Plots import save_plot
        save_plot(fig, ax, file_name='Overall cell number_multiple sample',
               folder_name=param.file_dir)
    else:
        import matplotlib.pyplot as pt
        pt.show()
