""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt

from neuron import h
import tempFunctions as tf
import tempConstants as tc


"""to plot all gating parameters of a LUT"""
# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75.pkl')
# tf.save_gatingplots(pkldict,r'test',A='all',reduced_xrange=False,Cm0=0.01)




""""LUT manipulations"""
# tf.downsample_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',[1,1,5,5,1,1])
#tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75_works.pkl',prints=True)
#tf.LUT_extend(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75.pkl')
#tf.random_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl')
#tf.upsample_LUT2(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00_ds.pkl',method='nearest')

"""compare original with upsampled one"""
# tf.compare_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',
#                r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\upsampled\(3,5,6,33,2,1)\realneuron_lookups_fs1.00_ds_us_linear.pkl')

#tf.LUT_comparison(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl', factor=[1,1,5,5,1,1], method='linear')

pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged_ext.pkl')
tf.save_gatingplots_overtones(pkldict,r'test',reduced_xrange=False,heatmap=1)

# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl', prints=True)
# tf.save_gatingplots(pkldict,r'test2',reduced_xrange=False)

# tf.LUT_extend_1overtone(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged.pkl')
