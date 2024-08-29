""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt
import re
import copy
import os
import sys
import argparse
#from MorphoSONIC.models import NeuronModel, SpatiallyExtendedNeuronModel, FiberNeuronModel, SingleCableFiber, MRGFiber
# from MorphoSONIC.models import SennFiber, Realnrn
from neuron import h
import tempFunctions as tf
import tempConstants as tc
# import PySONIC as ps


"""to read a pickle file"""

# tf.read_pickle('C:\\Users\\jgazquez\\PySONIC\\PySONIC\\lookups\\test_joa\\','realneuron_lookups_fs1.00_test.pkl')
# tf.read_pickle('C:\\Users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\','realneuron_lookups_fs1.00.pkl')



"""to plot all the effective gating parameter"""
# pkldict = tf.read_pickle('C:\\Users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\','RS_lookups_fs1.00.pkl',True)

# pkldict = tf.read_pickle('C:\\Users\\jgazquez\\PySONIC\\PySONIC\\lookups\\test_joa\\','realneuron_lookups_32nm_500kHz_fs1.00_22_01_2024_13_12_49.pkl',True)
#tf.save_gatingplots(pkldict,"gating")

# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl', prints=True)
# tf.save_gatingplots(pkldict,r'test2',reduced_xrange=False)
"mac"
# pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_lookups_32nm_500kHz_fs0.75_2024_03_08_10_56_59_merged_ext.pkl')
# tf.save_gatingplots(pkldict,"gating_ext",Cm0=0.01,reduced_xrange=1)
"with overtones"
# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged_ext.pkl')
# tf.save_gatingplots_overtones(pkldict,r'test',reduced_xrange=False,heatmap=1)


""""to merge a list of pickle files"""

# tf.merge_LUTlist("/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/",
#                  ['realneuron_lookups_32nm_500kHz_fs0.75_2024_02_27_16_13_33.pkl','realneuron_lookups_32nm_500kHz_fs0.75_2024_02_27_18_33_39.pkl',
#                   'realneuron_lookups_32nm_500kHz_fs0.75_2024_02_27_18_32_51.pkl','realneuron_lookups_32nm_500kHz_fs0.75_2024_02_27_16_14_16.pkl',
#                   'realneuron_lookups_32nm_500kHz_fs0.75_2024_02_28_11_27_30.pkl','realneuron_lookups_32nm_500kHz_fs0.75_2024_02_28_12_10_46.pkl',
#                   'realneuron_lookups_32nm_500kHz_fs0.75_2024_02_28_12_10_50.pkl','realneuron_lookups_32nm_500kHz_fs0.75_2024_02_28_12_11_12.pkl',
#                   'realneuron_lookups_32nm_500kHz_fs0.75_2024_02_28_12_12_34.pkl','realneuron_lookups_32nm_500kHz_fs0.75_2024_02_28_15_58_16.pkl'])



"""to extend a LUT"""
#tf.LUT_extend(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75.pkl')
"mac"
# tf.LUT_extend('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_lookups_32nm_500kHz_fs0.75_2024_03_08_10_56_59_merged.pkl')
"with overtones"
# tf.LUT_extend_1overtone(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged.pkl')



""""down/upsampling of a LUT"""
# tf.downsample_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',[1,1,5,5,1,1])
#tf.upsample_LUT2(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00_ds.pkl',method='nearest')



"""compare original LUT with upsampled one"""
# tf.compare_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',
#                r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\upsampled\(3,5,6,33,2,1)\realneuron_lookups_fs1.00_ds_us_linear.pkl')

"in one step"
#tf.LUT_comparison(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl', factor=[1,1,5,5,1,1], method='linear')



""""to change the units of the LUT"""
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\soma_lookups_1.pkl")
# print(pkldict['tables']['alpham_Ih'][0,0,0,0])
# for e in pkldict['tables']:
#     if '_' in e:
#         #print(e)
#         pkldict['tables'][e] *= 1e3
# print(pkldict['tables']['alpham_Ih'][0,0,0,0])
# tf.load_pickle(pkldict,r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\soma_lookups.pkl")



""""to lookup a value in the tables"""
# print(tf.lookup_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl','alpham_CaHVA'))



""""to test test.hoc"""

# from neuron import h,gui
# h.load_file('test.hoc')



"""to compare the different models"""

# print(f"{'-'*10}Fiber{'-'*10}")
# Fiber = SennFiber(2,2)
# print(f"{'-'*10}MRG{'-'*10}")
# MRG = MRGFiber(fiberD = 2, nnodes = 2)
# print(Fiber.nodes)

# print(f"{'-'*10}Realnrn{'-'*10}")
# realnrn = Realnrn(cell_nr=7,se=0)  
# print(realnrn.seclist())

#print(list(set(dir(Fiber)) & set(dir(MRG))))



""" to test the functionality of the passive mechanism"""

# soma = h.Section(name='soma')
# soma.insert('pas')
# soma.insert('pas_eff')
# print(soma.psection())


"""to check if there are 'ghost sections' present """

#check if duplicate sections:
# soma_dupl = []
# for i,sec in enumerate(h.allsec()):
#     #print('first sec: ',sec) if (i==0) else None
#     if 'soma' in str(sec):
#         soma_dupl.append(sec)
#print('last sec: ',sec)
#print('somas:\t',soma_dupl) # to check if soma is defined multiple times


"""to plot the expsyn and exp2syn mechanisms"""

# from neuron import h, gui
# import matplotlib.pyplot as plt

# # Load NEURON standard run library
# h.load_file('stdrun.hoc')

# # Create a section
# soma = h.Section(name='soma')

# # Create an ExpSyn and an Exp2Syn at different locations in the soma
# expsyn = h.ExpSyn(0.3, sec=soma)
# exp2syn = h.Exp2Syn(0.7, sec=soma)

# # Set parameters for ExpSyn
# expsyn.tau = 5.0  # Time constant in milliseconds
# expsyn.e = 0      # Reversal potential in mV

# # Set parameters for Exp2Syn
# exp2syn.tau1 = 1.0  # Rise time constant in milliseconds
# exp2syn.tau2 = 5.0  # Decay time constant in milliseconds
# exp2syn.e = 0       # Reversal potential in mV
# # exp2syn.gmax = 0.001  # Set the maximum conductance

# # Create NetCons to deliver currents to the synapses
# stim_expsyn = h.NetStim()
# stim_expsyn.number = 1
# stim_expsyn.start = 10
# stim_expsyn.interval = 0
# stim_expsyn.noise = 0

# stim_exp2syn = h.NetStim()
# stim_exp2syn.number = 1
# stim_exp2syn.start = 20
# stim_exp2syn.interval = 0
# stim_exp2syn.noise = 0

# netcon_expsyn = h.NetCon(stim_expsyn, expsyn)
# netcon_expsyn.weight[0] = 0.001  # Set the weight for ExpSyn

# netcon_exp2syn = h.NetCon(stim_exp2syn, exp2syn)
# netcon_exp2syn.weight[0] = 0.001  # Set the weight for Exp2Syn

# # Record time and conductances
# time_vec = h.Vector()
# expsyn_g_vec = h.Vector()
# exp2syn_g_vec = h.Vector()

# time_vec.record(h._ref_t)
# expsyn_g_vec.record(expsyn._ref_g)
# exp2syn_g_vec.record(exp2syn._ref_g)

# # Run the simulation
# h.tstop = 50  # Simulation time in milliseconds
# h.dt = 0.1    # Time step in milliseconds
# h.steps_per_ms = 1/h.dt
# h.run()

# # Plot the results
# plt.figure(figsize=(8, 6))
# plt.plot(time_vec, expsyn_g_vec, label='ExpSyn Conductance', color='red')
# plt.plot(time_vec, exp2syn_g_vec, label='Exp2Syn Conductance', color='blue')
# plt.xlabel('Time (ms)')
# plt.ylabel('Conductance (uS)')
# plt.legend()
# plt.show()



""""to import a cell with a larger stack size"""

# nrn_options = "-NSTACK 100000 -NFRAME 0"
# nrn_options = "-nogui -NSTACK 3000 -NFRAME 525"
# os.environ["NEURON_MODULE_OPTIONS"] = nrn_options
# from neuron import h, gui

#h("NSTACK_size = 10000")
#h("nrngui -NSTACK 100000 -NFRAME 20000 init.hoc ")

# h.load_file('init.hoc')

# cell_nr = 6
# h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference -> this also sets myelinate_axon = 1
# h.cell_chooser(cell_nr)

# Create an instance of the template
#cell = h.cADpyr229_L23_PC_5ecbf9b163(0)



"""to test a parsing function"""

#     parser example
# def main():
#     # Create an ArgumentParser object
#     parser = argparse.ArgumentParser(description="Custom Argument Example")

#     # Add custom arguments for 't', 'f', and 'z' with default values
#     parser.add_argument("-t", type=int, default=1, help="Custom argument 't'")
#     parser.add_argument("-f", type=float, default=0.0, help="Custom argument 'f'")
#     parser.add_argument("-z", type=str, default="default_value", help="Custom argument 'z'")

#     # Parse the command-line arguments
#     args = parser.parse_args()

#     # Access the values of 't', 'f', and 'z' from the parsed arguments
#     t_value = args.t
#     f_value = args.f
#     z_value = args.z

#     # Now you can use these values in your script
#     print(f"t = {t_value}")
#     print(f"f = {f_value}")
#     print(f"z = {z_value}")

# if __name__ == "__main__":
#     main()
# increases the stack in able to load cell_nr = 6
# import subprocess
# import neuron
# subprocess.run('neuron',shell=True)
# subprocess.run('nrngui -NSTACK 100000 -NFRAME 20000 init.hoc',shell=True)