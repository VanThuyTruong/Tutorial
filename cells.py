#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Van Thuy Truong 

Python 3.8 : Mon Apr 06 17:41:36 2020
"""

from cell import Cell
import numpy as np

class Cells(object):
    """
    Set up a tumour cell population
    """
    def __init__(self, param):
        self.cells = np.array([Cell(param.sample[i]) for i in\
                               range(param.ncell)]) #array of cell objects according to prior samplep
        #print(self.cells)
        self.pERK = np.array([cell.pERK0 for cell in self.cells]) #initialize numpy array with float
        #print(self.pERK)
        self.ncell = param.ncell
        self.active_cells = list(range(self.ncell))
        self.update_pool(param)
        
    def update_pool(self, param): #find tumour cells in the division, quiescence and death pool
        self.rep_finder(param)
        self.death_finder(param)
        self.n_quiesence=self.ncell-self.n_death-self.n_rep

    def rep_finder(self, param): #find cells with pERK values over the divsion threshold
        self.rep_indices = np.where(self.pERK > param.trep)[0] #because we want array of indices not the returned tuple which includes data type
        self.n_rep = len(self.rep_indices)
        self.division_sumpERK = np.sum(self.pERK[self.rep_indices])

    def death_finder(self, param): #find cells with pERK values below the death threshold
        self.death_indices = np.where(self.pERK < param.tdeath)[0]
        self.n_death = len(self.death_indices)
        self.death_sumpERK = np.sum(param.maxbaseline-self.pERK[self.death_indices])

        
    def division(self, param, time=None): #divsion reaction, choose one cell for division
        weights = self.pERK[self.rep_indices]/self.division_sumpERK
        Cind = np.random.choice(self.rep_indices, p=weights)
        self.cells = np.append(self.cells,
                        np.array([Cell(self.pERK[Cind], time)]))
        self.active_cells.append(len(self.cells)-1) #adds a new cell number to our cell list
        self.pERK = np.append(self.pERK,
                        np.array([self.pERK[Cind]]))
        self.ncell = len(self.active_cells)
                                    
    def death(self, param): #death reaction, choose one cell for death
        weights = (param.maxbaseline-self.pERK[self.death_indices])/self.death_sumpERK
        Cind = np.random.choice(self.death_indices, p=weights)
        self.active_cells.pop(Cind) #removes the current cell from active cell list
        self.ncell = len(self.active_cells)


    def update_pERK(self, pERK_decrease, time=None): #apply pERK decrease caused by the drug to remaining tumour cells
        self.pERK = np.empty(self.ncell)
        for i, active_cell_index in enumerate(self.active_cells):
            #i is the 0,1,...,n index, while j is the index of an active cells
            self.cells[active_cell_index].updatepERK(pERK_decrease, time)
            self.pERK[i] = self.cells[active_cell_index].pERKt #call attribute with of cell with .pERKt
    
    def plot_pERK(self,param,ABM):
        """
        Plots the pERK values of all cells in time of a given simulation, showing
        when each cell is created or dies.
        """
        import matplotlib.pyplot as pt
        from PPP.Plots import plot_setup
        from matplotlib.lines import Line2D
        for i in range(ABM.sample_size):
            fig, ax= plot_setup('Time (hour)', 'pERK value', 
                                     scale=param.plot_scale)
            legend_elements = [Line2D([0], [0], color='black', lw=1, 
                                      label='Individual pERK profile of one agent'),
                       Line2D([0], [0], marker='x', color='black', label='Timepoint of the Gillespie algorithm',
                              markerfacecolor='g', markersize=15)]
            for cell in self.cells:
                pt.plot(cell.times, cell.pERKt_values, '-x')
            ax.set_ylim([param.minbaseline, param.maxbaseline])
            #fig.subplots_adjust(bottom=0.3)
            ax.legend(handles=legend_elements,
                      fontsize=16*param.plot_scale)
            pt.show()

            if param.save_plot:
                import os
                from PPP.Plots import save_plot
                save_plot(fig, ax, file_name='pERK_tracking_{}'.format(i+1),
                      folder_name=os.path.join(
                              param.file_dir, 'Sample%s'%(i+1)
                              )
                      )
            else:
                import matplotlib.pyplot as pt
                pt.show()
            
