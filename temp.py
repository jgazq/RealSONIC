import numpy as np
import time
import datetime
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import sys
import math
from scipy.interpolate import griddata, interpn
nrn_options = "-NSTACK 10000 -NFRAME 525"
#Aberra recommends: -NSTACK 100000 -NFRAME 20000
os.environ["NEURON_MODULE_OPTIONS"] = nrn_options
from neuron import h
h("NSTACK_size = 10000")
h.load_file('init.hoc')
#h.load_file('thresh4.hoc') 
#h.load_file('get_es2.hoc') #for TMS
#h.load_file('interp_coordinates.hoc') #for TMS (to calc Dx,Dy and Dz)


current_time = datetime.datetime.now()
now = datetime.datetime.strftime(current_time,'%d_%m_%Y_%H_%M_%S')
PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__))


cell_nr = 7
h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
h.cell_chooser(cell_nr)
#print(h.topology()) #print this to decide the code for the cell below
if cell_nr == 2:
    cell = h.bNAC219_L1_NGCDA_e7cec642c3[0] # for cell = 2
elif cell_nr == 3:
    cell = h.bNAC219_L1_NGCDA_46b45974f4[0] # for cell = 3
elif cell_nr == 7:
    cell = h.cADpyr229_L23_PC_8ef1aa6602[0] # for cell = 7


"--------------------------------------------------------------------------------------------------------------------"
#print(h.setParamsAdultRat(1,2,3))
print("myelinate ax: ",h.myelinate_ax) 
#h.nseg = 3
print("nseg: ")
print(h.nseg)

print(h.stim_xtra)

print("topology: ")
print(h.topology())

#h.color_plotmax(1)

print(h.v_init)
h.run() # to run the simulation
#h.graph()
#h.Vector() #to create a vector
#time.sleep(20)

# to calculate the threshold
#all_functions.h.load_file('thresh4.hoc')
#all_functions.h.threshold(0,1,0.1)