""""just trying some things with the PySONIC module to figure out how it works"""

"-----IMPORTS-----"

# import numpy as np
#import time
#import datetime
import matplotlib.pyplot as plt
#import scipy.io as sio
#import os
#import sys
#import math
#import seaborn as sns
#import logging
#from scipy.interpolate import griddata, interpn
import pandas as pd
from neuron import h, gui
#h("NSTACK_size = 10000")
h.load_file('init.hoc')

# import tempConstants as tc
import tempFunctions as tf
# import prev.functions as fs
# import prev.Interp3Dfield as tt
# import PySONIC as ps
# import MorphoSONIC as ms
from PySONIC.core import PulsedProtocol
from MorphoSONIC.models import RealBase, RealNeuronModel,MRGFiber
from MorphoSONIC.plt import plotFiberXCoords, SectionCompTimeSeries
from MorphoSONIC.sources import IntracellularCurrent, ExtracellularCurrent

"-----CODE-----"
#---create concentric planar disk transducer source
# psource_transd = ms.PlanarDiskTransducerSource(x=(0,0,0),f=0.5e6,rho=tc.RHO,c=tc.C,r=(2e-3)/2,u=1)
# grid_type = 'concentric'
# tf.plt_transdistr(psource_transd,grid_type)

#---read out metadata from lookup files (pickle files)
# pkl_txt =tf.read_pickle(r'C:\\Users\\jgazquez\\PySONIC\\PySONIC\\lookups\\',r'Cm_lkp_32nm.pkl')
# print(pkl_txt['tables']['Cm_rel'].shape)
# pkl_txt =tf.read_pickle(r'C:\\Users\\jgazquez\\PySONIC\\PySONIC\\lookups\\',r'pas_Cm0_1.0uF_cm2_ELeak_-70.0mV_lookups_32nm_20kHz.pkl')
# print("RS\n\n")
# pkl_txt = tf.read_pickle('c:\\users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\test_joa\\','RS_lookups_fs1.00_test.pkl')
# print("realneuron\n\n")
# pkl_txt = tf.read_pickle('c:\\users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\test_joa\\','realneuron_lookups_fs1.00.pkl')

cell_nr = 7
h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
h.cell_chooser(cell_nr)
#get cell name folder based on the cell that has been chosen
cell_name = h.cell_names[cell_nr].s
#print(h.topology()) #print this to decide the code for the cell below
if cell_nr == 2:
    cell = h.bNAC219_L1_NGCDA_e7cec642c3[0] # for cell = 2
elif cell_nr == 3:
    cell = h.bNAC219_L1_NGCDA_46b45974f4[0] # for cell = 3
elif cell_nr == 7:
    cell = h.cADpyr229_L23_PC_8ef1aa6602[0] # for cell = 7

#nodes = tf.coord_dict()
compartments = tf.coord_dict_type()
#print([e[0] for e in nodes.values()])
#print(nodes_type.keys())
nnodes = 21
fiberD = 10e-6  # um

""""comparing the plotFiberXCoords of MRGFiber with RealNeuronModel like is shown in the mrg.ipynb notebook"""
# fig = plotFiberXCoords(MRGFiber(fiberD, nnodes))
# plt.show()
# fig = plotFiberXCoords(RealNeuronModel(compartments))
# plt.show()


def plotSimIClamp(fiber, pp, I=None, ylim=None, title=None):
    if I is not None:
        polarity = 'cathode' if I < 0 else 'anode'
    else:
        polarity = 'anode'    
    source = IntracellularCurrent(fiber.central_ID, I=I, mode=polarity)
    data, meta = fiber.simulate(source, pp)
    fig = SectionCompTimeSeries([(data, meta)], 'Vm', fiber.nodeIDs).render()
    if ylim is not None:
        fig.axes[0].set_ylim(ylim)
    if title is not None:
        fig.axes[0].set_title(title)
    return fig

def plotSimExtSource(fiber, source, pp, ylim=None, title=None):
    data, meta = fiber.simulate(source, pp)
    fig = SectionCompTimeSeries([(data, meta)], 'Vm', fiber.nodeIDs).render()
    if ylim is not None:
        fig.axes[0].set_ylim(ylim)
    if title is not None:
        fig.axes[0].set_title(title)
    return fig

pp = PulsedProtocol(100e-6, 3e-3) 
# fig = plotSimIClamp(MRGFiber(fiberD,nnodes), pp, title='MRGFiber(fiberD,nnodes)')
# plt.show()
fig = plotSimIClamp(RealNeuronModel(compartments), pp, title='RealNeuronModel(compartments)')
plt.show()

