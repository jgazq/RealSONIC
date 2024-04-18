""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt

from neuron import h
import tempFunctions as tf
import tempConstants as tc


"""to plot all gating parameters of a LUT"""
# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75_works2_us_linear.pkl')
# tf.save_gatingplots(pkldict,r'test2\upsampled',A='all',reduced_xrange=False,Cm0=0.01)


""""LUT manipulations"""
#tf.downsample_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',[1,1,10,10,1,1])
#tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75_works.pkl',prints=True)
#tf.LUT_extend(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75.pkl')
#tf.random_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl')
tf.upsample_LUT2(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00_ds.pkl',method='cubic')

"""compare original with upsampled one"""
#tf.compare_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl', r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00_ds_us_cubic.pkl')

#tf.LUT_comparison(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl', factor=[1,1,5,5,1,1], method='linear')