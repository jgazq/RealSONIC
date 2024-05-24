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
#tf.plot_astim(r'output_csv_archive\realistic_cort_realneuron_16nm_fs75%_f_100kHz_A_1.00MPa_CW_tstim_10ms_toffset_3ms.csv')

#tf.LUT_to_LUT2(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged.pkl',1)

#pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged_LUT2.pkl",prints=True)


# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged_unext_2.pkl")
# tf.save_gatingplots(pkldict,r'test')

tf.plot_astim(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 6\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_3ms_tstart_10ms\2024_05_22_15_26_50_basal0.csv")
