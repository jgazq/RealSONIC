""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt
import os
import csv

from neuron import h
import tempFunctions as tf
import tempConstants as tc

#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 13\csv\ASTIM_soma_CW_64nm_f_100kHz_A_600.00kPa_tstim_98ms_toffset_10ms_tstart_10ms_fs75%_sonic\2024_09_05_11_29_27.csv",debug=True,variables=['Q','V','i_net']) #,separate=1,folder=r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 13\realistic vs soma\\"
#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 13\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_98ms_toffset_10ms_tstart_10ms\2024_09_05_11_31_01_soma0.csv",debug=True) #,separate=1,folder=r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 13\realistic vs soma\\"
#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 12\csv\ASTIM_soma_CW_64nm_f_100kHz_A_600.00kPa_tstim_0s_tstart_10ms_fs75%_sonic\2024_09_05_14_10_20.csv",debug=True,variables=['Q','V','i_net']) #,separate=1,folder=r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 13\realistic vs soma\\" #only SKv3_1
#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_12_14_57_31_node0.csv")#,compare=1) #,separate=1,folder=r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 13\realistic vs soma\\"

#tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 15\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_1ps_toffset_1ps_tstart_10us\2024_09_10_16_54_24_myelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 15\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_1ps_toffset_1ps_tstart_10us\2024_09_10_16_54_26_unmyelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 15\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_1ps_toffset_1ps_tstart_10us\2024_09_10_16_54_28_node0.csv"])
tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_11_15_51_31_myelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_11_15_51_32_node0.csv"],debug=1)
#tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_12_14_57_30_myelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_12_14_57_31_node0.csv"])

#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 9\csv\ASTIM_soma_CW_64nm_f_100kHz_A_600.00kPa_tstim_98ms_toffset_10ms_tstart_10ms_fs75%_sonic\2024_08_07_14_29_35.csv",debug=True)
#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 9\csv\ASTIM_soma_CW_64nm_f_100kHz_A_600.00kPa_tstim_98ms_toffset_10ms_tstart_10ms_fs75%_sonic\2024_08_07_14_53_36.csv",debug=True)

#COMPARE OVERTONES WITH 0 OVERTONES
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups.pkl",prints=True)
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\OneDrive - UGent\PhD\jgazquez - backup\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged_ext.pkl",prints=True)

#tf.LUT_to_LUT2(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged.pkl",1)

#tf.save_gatingplots_group(pkldict,r'test',Cm0=None, reduced_yrange=False)

#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 10\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_08_12_09_56_01_soma0.csv",debug=1)#,variables = ['V','Q','i_net'])


"""compare original LUT with upsampled one"""
# tf.compare_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',
#                r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\upsampled\(3,5,6,33,2,1)\realneuron_lookups_fs1.00_ds_us_linear.pkl')

"""plot the variables during a time acoustic stimulation"""
#tf.plot_astim(r'output_csv_archive\realistic_cort_realneuron_16nm_fs75%_f_100kHz_A_1.00MPa_CW_tstim_10ms_toffset_3ms.csv')

"""LUT to LUT2 (where the Cm0-variable is removed an divided into 2 LUT's)"""
#tf.LUT_to_LUT2(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged.pkl',1)

#pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_merged_LUT2.pkl",prints=True)

"""plotting the gating parameters"""
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups.pkl")
# tf.save_gatingplots(pkldict,r'test',A='all',Cm0=None)

"""analysis of the Q/V mismatch resulting in unphysical values for C"""
#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 7\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_06_25_17_57_25_soma0.csv",debug=True,variables=['V','Q'],separate=1)
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