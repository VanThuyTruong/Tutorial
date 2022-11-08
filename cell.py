#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Van Thuy Truong 

Python 3.8 : Mon Apr 06 17:41:36 2020
"""

class Cell(object):
    """
    Set up one tumour cell
    """
    def __init__(self, pERK0, time=0): #time of birth (time) is 0 by default, unless specified otherwise
        self.pERK0=pERK0  #pERK value of the one tumour cell
        self.pERKt=pERK0
        self.pERKt_values=[self.pERK0]
        self.times=[time]
#        print('pERKt0',self.pERKt)

    def updatepERK(self, pERK_decrease, time=None): #apply pERK decrease caused by the drug to the tumour cell
        self.pERKt = (100-pERK_decrease)*self.pERK0/100
        self.pERKt_values.append(self.pERKt)
        if time is not None:
            self.times.append(time)
