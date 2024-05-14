""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt

from neuron import h
import tempFunctions as tf
import tempConstants as tc


"""compare original LUT with upsampled one"""
# tf.compare_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',
#                r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\upsampled\(3,5,6,33,2,1)\realneuron_lookups_fs1.00_ds_us_linear.pkl')

"""plot the variables during a time acoustic stimulation"""
tf.plot_astim(r'output_csv_archive\realistic_cort_realneuron_16nm_fs75%_f_100kHz_A_1.00MPa_CW_tstim_10ms_toffset_3ms.csv')

#tf.LUT_to_LUT2(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged.pkl', remove_zeros=1)

#tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75.pkl",prints=True)