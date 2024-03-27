""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt

from neuron import h
import tempFunctions as tf
import tempConstants as tc

# pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_opt_lookups_32nm_500kHz_[0_321]_fs0.75_2024_03_14_17_20_14_ext.pkl')
# #pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_lookups_16nm_32nm_500kHz_1000kHz_fs0.75_merged_2024_03_19_15_59_08_merged_ext.pkl')
# tf.save_gatingplots(pkldict,"test",Cm0=0.01,reduced_xrange=1)

# pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_lookups_16nm_32nm_500kHz_1000kHz_fs0.75_merged_2024_03_19_15_59_08_merged_ext.pkl')
# #pkldict = tf.read_pickle('/Users/joaquin/Documents/python-virtual-environments/PySONIC/PySONIC/lookups/test_joa/realneuron_lookups_16nm_32nm_500kHz_1000kHz_fs0.75_merged_2024_03_19_15_59_08_merged_ext.pkl')
# tf.save_gatingplots_new(pkldict,"gating/sweeps/radius",1e9,'nm',a="all",reduced_xrange=1)

a = [1,2,3,4,8,9]
print(a[0:-1])
