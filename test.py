""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt
import os

#this is added to work on macOS
if "DISPLAY" in  os.environ:
    del os.environ['DISPLAY']

from neuron import h
import tempFunctions as tf
import tempConstants as tc

#tf.LUT_extend('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_opt_lookups_32nm_500kHz_[0_321]_fs0.75_2024_03_14_17_20_14.pkl')
pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_opt_lookups_32nm_500kHz_[0_321]_fs0.75_2024_03_14_17_20_14_ext.pkl')
tf.save_gatingplots(pkldict,"gating_opt",Cm0=0.01,reduced_xrange=1)

