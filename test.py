""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt
import os

from neuron import h
import tempFunctions as tf
import tempConstants as tc


"""compare original LUT with upsampled one"""
# tf.compare_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',
#                r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\upsampled\(3,5,6,33,2,1)\realneuron_lookups_fs1.00_ds_us_linear.pkl')

"""plot the variables during a time acoustic stimulation"""
#tf.plot_astim(r'output_csv_archive\realistic_cort_realneuron_16nm_fs75%_f_100kHz_A_1.00MPa_CW_tstim_10ms_toffset_3ms.csv')

"""LUT to LUT2 (where the Cm0-variable is removed an divided into 2 LUT's)"""
#tf.LUT_to_LUT2(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged.pkl',1)

#pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged_LUT2.pkl",prints=True)

"""plotting the gating parameters"""
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged_LUT2.pkl")
# tf.save_gatingplots(pkldict,r'test',A='all',Cm0=None)

"""analysis of the Q/V mismatch resulting in unphysical values for C"""
tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 7\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_06_21_17_21_11_basal0.csv",debug=True)
#tf.plot_astim(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 7\csv\ASTIM_RS_CW_32nm_f_500kHz_A_100.00kPa_tstim_100ms_toffset_50ms_tstart_10ms_sonic\2024_06_20_11_58_35.csv", debug=1)
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_64nm_100kHz_fs0.75.pkl")
# Q = pkldict['refs']['Q']*1e5
# V = pkldict['tables']['V'][2,0,50,:,0]
# testo = np.interp(-0.00955699, Q, V)
# print(f'np.interp(-0.00955699, Q, V) = {testo}')

"""to modify the NMODL files"""
# for root, dirs, files in os.walk(r"C:\Users\jgazquez\RealSONIC\mechanisms_versions\22 apr 2024\eff - Copy"): #go through all files in the mechanisms folders (all depths)
#     for fil in files:
#         file = os.path.join(root,fil)
#         with open(file,'r') as fyle:
#             flist = list(fyle)
#         flist2 = []
#         for e in flist:
#             if 'print' in e:
#                 continue
#             if '_2' in fil and ('alpha' in e or 'beta' in e) and not 'Stoch' in fil:
#                 repl = fil.split('_eff')[0].replace('_','')
#                 print(e.replace(repl,repl+'2'))
#                 flist2.append(e.replace(repl,repl+'2'))
#             else:
#                 flist2.append(e)
#         with open(file.replace('Copy','Copy2'),'w') as fyle2:
#             fyle2.writelines(flist2)