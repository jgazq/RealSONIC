""""just testing some things"""

import numpy as np
import re
import copy
import os
import sys
import argparse
#from MorphoSONIC.models import NeuronModel, SpatiallyExtendedNeuronModel, FiberNeuronModel, SingleCableFiber, MRGFiber
# from MorphoSONIC.models import SennFiber, Realnrn
#from neuron import h
import tempFunctions as tf
import tempConstants as tc
# import PySONIC as ps

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

"""to read a pickle file"""
# tf.read_pickle('C:\\Users\\jgazquez\\PySONIC\\PySONIC\\lookups\\test_joa\\','realneuron_lookups_fs1.00_test.pkl')
# tf.read_pickle('C:\\Users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\','realneuron_lookups_fs1.00.pkl')
np.set_printoptions(suppress=True)
print(len(np.append(0,np.logspace(np.log10(0.1), np.log10(600), 50))))
