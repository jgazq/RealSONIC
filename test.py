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

#pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_lookups_32nm_500kHz_fs0.75_2024_03_08_10_56_59_merged_ext.pkl',prints=1)
# pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_lookups_32nm_500kHz_fs0.75_2024_02_28_12_11_12.pkl',prints=1)
# print(tf.tcomp_run_time("126.99542808644578 dayss", "1 day, 20:18:59.47"))
tf.write_csv()


#tf.save_gatingplots(pkldict,"gating_ext",Cm0=0.02)