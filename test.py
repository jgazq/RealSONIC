""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import pprint
from scipy.optimize import brentq, least_squares, fmin, minimize
from scipy.linalg import block_diag

from neuron import h, numpy_element_ref
import tempFunctions as tf
import tempConstants as tc

# tf.LUT_Jacob(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_1ov.pkl")
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_1ov_Jac.pkl",prints=1)
# print(pkldict['Jacobians'].keys())
# for e in pkldict['Jacobians'].values():
#     print(e.shape)

a = [x for i in range(1, 10) for x in ((5,5),(2,2))]
print(a)


#tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_1ov.pkl",prints=1)

#tf.minimize_ov(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups.pkl", r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_1ov.pkl")


# tf.LUT_to_LUT2_1overtone(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged.pkl",1)

#tf.save_gatingplots_overtones(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged_ext.pkl",r'test')

"""activation zones: thresholds"""
# tf.analyze_over_sections(r"C:\Users\jgazquez\RealSONIC\pickledump_RealSIM\pickledump\stimulation_zones_anlysis\dump_75.0__32.0nm_500.0kHz_112.99551984243782kPa_100.0ms_10.0ms_1000.0Hz_0.53DC.pkl")
# tf.analyze_over_sections(r"C:\Users\jgazquez\RealSONIC\pickledump_RealSIM\pickledump\stimulation_zones_anlysis\dump_75.0__32.0nm_500.0kHz_119.55039638327578kPa_100.0ms_10.0ms_1000.0Hz_0.53DC.pkl")

"""to combine two seperated calculated LUT"""
# pkl_Cm01 = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\Cm0_separate\realneuron_lookups_32nm_500kHz_[0_-1]_fs0.75_0.01_2024_12_03_09_55_02.pkl")
# pkl_Cm02 = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\Cm0_separate\realneuron_lookups_32nm_500kHz_[0_-1]_fs0.75_0.02_2024_12_02_18_26_19.pkl")
# del pkl_Cm01['tables']['tcomp']
# del pkl_Cm02['tables']['tcomp']

# for e,f in pkl_Cm01['tables'].items():
#     nonzeros = np.nonzero(f[0,0,0,:,0])[0]
#     pkl_Cm01['tables'][e] = f[:,:,:,:,0,:][:,:,:,nonzeros]

# for e,f in pkl_Cm02['tables'].items():
#     nonzeros2 = np.nonzero(f[0,0,0,:,0])[0]
#     pkl_Cm01['tables'][e+'2'] = f[:,:,:,:,0,:][:,:,:,nonzeros2]

# pkl_Cm01['tables']['Q_ext'] = pkl_Cm02['refs']['Q'][nonzeros2]
# pkl_Cm01['refs']['Q'] = pkl_Cm01['refs']['Q'][nonzeros]
# print(pkl_Cm01['tables']['V'].shape, pkl_Cm02['tables']['V'].shape)
# tf.load_pickle(pkl_Cm01,r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\Cm0_separate\realneuron_lookups_32nm_500kHz_[0_-1]_fs0.75_LUT2.pkl")
""""""

#tf.plot_titration_curves_cells(r"C:\Users\jgazquez\RealSONIC\titrate.pkl")


#tf.plot_astim2(r"C:\Users\jgazquez\RealSONIC\pickledump\75.0%_16.0nm_100.0kHz_20.0ms_20.0ms_100.0Hz_0.5DC\csv\dump_75.0%_16.0nm_100.0kHz_20.0ms_20.0ms_100.0Hz_0.5DC_600.0kPa.csv")

#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 10\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_08_12_09_56_01_soma0.csv",debug=1)#,variables = ['V','Q','i_net'])

# tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 20\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_25_17_22_43_unmyelin0.csv") #
# pkl = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_2024_09_25_17_48_48.pkl")
# print(pkl['refs']['Q'])
# print(pkl['tables']['V'])

"""empty titrate.pkl, read from txt and plot the titration curve"""
# tf.load_pickle({},r"C:\Users\jgazquez\RealSONIC\titrate.pkl")

# tf.txt_to_titration(r"C:\Users\jgazquez\RealSONIC\titrate", r"C:\Users\jgazquez\RealSONIC\titrate.pkl", save=1)

# tit = tf.read_pickle(r"C:\Users\jgazquez\RealSONIC\titrate.pkl")
#pprint.pprint(tit)
# tf.plot_titration_curves(r"C:\Users\jgazquez\RealSONIC\titrate.pkl")
""""""

# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\RS_lookups.pkl",prints=True)
# LUV = tf.lookup_LUT(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\RS_lookups_fs1.00.pkl",lookups=[32*1e-9, 500*1e3, 600*1e3, 50*1e-5,1.],Lemaire=1)
# print(LUV)
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\Downloads\realneuron_lookups_2024_12_03_17_44_00.pkl")
# pprint.pprint(pkldict)


# tf.LUT_extend(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\Cm0 separate\realneuron_lookups_32nm_500kHz_[0_-1]_fs0.75_0.01_2024_12_03_09_55_02.pkl")
# tf.merge_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\Cm0 separate',
#              r"\realneuron_lookups_32nm_500kHz_[0_-1]_fs0.75_0.01_2024_12_03_09_55_02_ext.pkl",
#              r'\realneuron_lookups_32nm_500kHz_[0_-1]_fs0.75_0.02_2024_12_02_18_26_19.pkl')
# tf.LUT_to_LUT2(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\Cm0_together\realneuron_lookups_2024_12_05_03_57_14.pkl",1)

"""titrate input data for HPC worker submission"""
# with open ("titrate_data.txt",'w') as tf:
#     i = 0
#     for a in [32]:
#         for f in [500]:
#             for PRF in [1000]:
#                 for DC in range(15,101):
#                     tf.write(f'7,{a},{f},{PRF},{DC}\n')
#                     i+=1
# print(i)

"""titrate data for worker wsub"""
# with open ("titrate_data.txt",'w') as tf:
#     i = 0
#     for cell in range(1,26):
#         for a in [32,64]:
#             for f in [500, 1000]:
#                 for PRF in [10, 100, 1000]:
#                     for DC in [50,100]:
#                         tf.write(f'{cell},{a},{f},{PRF},{DC}\n')
#                         i+=1
# print(i)


"""to plot two C curves by looking the values up in 2 LUT's"""
# Q_range = np.arange(-107,50,1)
# lookup = np.array([tf.lookup_LUT(r"C:\Users\jgazquez\RealSIM\PySONIC\lookups\downloaded\realneuron_lookups_32nm_500kHz_[0_321]_fs0.75_2024_03_14_17_20_14.pkl",'V',[32*1e-9, 500*1e3, 0*1e3, e*1e-5, 0.02]) for e in Q_range])
# Q_range = Q_range[abs(lookup) > 1e-5]
# lookup = lookup[abs(lookup) > 1e-5]
# plt.plot(Q_range,Q_range/ lookup,label='first')
# plt.xticks(fontsize=36)
# plt.yticks(fontsize=36)
# plt.show()

# lookup = np.array([tf.lookup_LUT(r"C:\Users\jgazquez\RealSIM\PySONIC\lookups\downloaded\realneuron_lookups_32nm_20kHz_[0_-1]_fs0.75_2024_05_29_16_44_56.pkl",'V',[32*1e-9, 500*1e3, 0*1e3, e*1e-5, 0.01]) for e in Q_range])
# Q_range = Q_range[abs(lookup) > 1e-5]
# lookup = lookup[abs(lookup) > 1e-5]
# plt.plot(Q_range,Q_range/ lookup,label='second')
# plt.xticks(fontsize=36)
# plt.yticks(fontsize=36)
# plt.legend()
# plt.show()


"""to investigate if first simulation influences the second onwards"""
# tf.analyze_over_sections(r"C:\Users\jgazquez\RealSONIC\pickledump\75.0%_64.0nm_100.0kHz_20.0ms_20.0ms_100.0Hz_1.0DC\pkl\dump_75.0%_64.0nm_100.0kHz_20.0ms_20.0ms_100.0Hz_1.0DC_600.0kPa.pkl", 
#                          r"C:\Users\jgazquez\RealSONIC\pickledump\75.0%_64.0nm_100.0kHz_20.0ms_20.0ms_100.0Hz_1.0DC\pkl\dump_75.0%_64.0nm_100.0kHz_20.0ms_20.0ms_100.0Hz_1.0DC_10.0kPa.pkl")
#tf.analyze_over_sections(r"C:\Users\jgazquez\RealSONIC\pickledump\pkldump\dump_75.0%_64.0nm_100.0kHz_20.0kPa_20.0ms_10.0ms_200.0Hz_1.0DC.pkl")
#tf.analyze_over_sections(r"C:\Users\jgazquez\RealSONIC\pickledump\pkldump\dump_75.0%_64.0nm_2000.0kHz_20.0kPa_100.0ms_10.0ms_100.0Hz_1.0DC.pkl")
#tf.analyze_over_sections(r"C:\Users\jgazquez\RealSONIC\pickledump\analysis\analysis1\dump_75.0%_32.0nm_500.0kHz_103.984375kPa_100.0ms_10.0ms_1000.0Hz_0.53DC.pkl")

"""plots the different probes parameters for various titration/amplitude values"""
# directory = r"C:\Users\jgazquez\RealSONIC\pickledump\75.0%_16.0nm_100.0kHz_20.0ms_100.0ms_100.0Hz_0.5DC\csv"
# dirlist = os.listdir(directory)
# dirlist_amp = [float(e.split('_')[-1].split('kPa')[0]) for e in dirlist]
# dirlist_sorted = [name for _, name in sorted(zip(dirlist_amp, dirlist))]

# for filename in dirlist_sorted:
#     f = os.path.join(directory, filename)
#     # checking if it is a file
#     if os.path.isfile(f):
#         print(f)
#         tf.plot_astim2(f)


""""Lemaire2019 imitation"""
# result_dict = {'refs': {'fs': [0.75], 'radius': [16*1e-9, 32*1e-9, 64*1e-9], 'freq': [100*1e3,1000*1e3], 'DC': [0.5, 1.],
#                          'PRF': [100., 200.]}, 'table': np.array([[[[[599.4573974609375,9.9951171875],[460.0042724609375,9.9951171875]],[[599.8675537109375,9.9951171875],[460.0042724609375,9.9951171875]]],[[[20.0048828125,20.0048828125],[20.0048828125,20.0048828125]],[[20.0048828125,20.0048828125],[20.0048828125,20.0048828125]]], [[[15.0048828125,17.5048828125],[15.0048828125, 15.0048828125]],[[19.9951171875, 20.0048828125],[20.0048828125, 19.9951171875]]]]])}
# result_matr = result_dict['table']
# for i,e in enumerate(result_dict['refs']['radius']):
#     print(result_matr[0,i,0,:,0])
#     plt.plot(result_dict['refs']['DC'],result_matr[0,i,0,:,0],label=e*1e9)
#     plt.yscale('log')
# plt.xlabel('DC')
# plt.ylabel('Amplitude [kPa]')
# plt.legend()
# plt.show()

# for i,e in enumerate(result_dict['refs']['freq']):
#     print(result_matr[0,1,i,:,0])
#     plt.plot(result_dict['refs']['DC'],result_matr[0,1,i,:,0],label=e*1e-3)
#     plt.yscale('log')
# plt.xlabel('DC')
# plt.ylabel('Amplitude [kPa]')
# plt.legend()
# plt.show()

# for i,e in enumerate(result_dict['refs']['PRF']):
#     print(result_matr[0,1,0,:,i])
#     plt.plot(result_dict['refs']['DC'],result_matr[0,1,0,:,i],label=e)
#     plt.yscale('log')
# plt.xlabel('DC')
# plt.ylabel('Amplitude [kPa]')
# plt.legend()
# plt.show()
""" """

"""connected plots of FEARS poster"""
# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_01_soma0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_03_unmyelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_06_apical0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_08_basal0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_10_axon0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_12_myelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_14_node0.csv"])

# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_06_apical0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_08_basal0.csv",])

# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_03_unmyelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_12_myelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_14_node0.csv"])

# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_01_soma0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_12_10_axon0.csv",])

"""disconnected plots of FEARS poster"""
# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_03_soma0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_05_unmyelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_07_apical0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_09_basal0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_11_axon0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_13_myelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_15_node0.csv"])

# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_07_apical0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_09_basal0.csv",])

# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_05_unmyelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_13_myelin0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_15_node0.csv"])

# tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_03_soma0.csv",
#                         r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 23\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_11_08_16_17_11_axon0.csv",])


"""plotting the simulation of multiple sections on 1 plot"""
#tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 21\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_26_15_53_05_myelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 21\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_26_15_53_07_node0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 21\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_26_15_53_10_soma0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 21\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_26_15_53_12_unmyelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 21\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_26_15_53_15_apical0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 21\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_26_15_53_18_basal0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 21\csv\realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_26_15_53_20_axon0.csv"])
#tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_11_15_51_31_myelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_11_15_51_32_node0.csv"],debug=1)
#tf.plot_astim_sections([r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_12_14_57_30_myelin0.csv",r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 16\csv\mini_realistic_cort_realneuron_64nm_fs75%_f_100kHz_A_600.00kPa_CW_tstim_100ms_toffset_10ms_tstart_10ms\2024_09_12_14_57_31_node0.csv"])

#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 9\csv\ASTIM_soma_CW_64nm_f_100kHz_A_600.00kPa_tstim_98ms_toffset_10ms_tstart_10ms_fs75%_sonic\2024_08_07_14_29_35.csv",debug=True)
#tf.plot_astim2(r"C:\Users\jgazquez\OneDrive - UGent\PhD\Figures\self_made\run_realistic_astim output\try 9\csv\ASTIM_soma_CW_64nm_f_100kHz_A_600.00kPa_tstim_98ms_toffset_10ms_tstart_10ms_fs75%_sonic\2024_08_07_14_53_36.csv",debug=True)

"""COMPARE OVERTONES WITH 0 OVERTONES"""
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups.pkl",prints=True)
# pkldict = tf.read_pickle(r"C:\Users\jgazquez\OneDrive - UGent\PhD\jgazquez - backup\PySONIC\PySONIC\lookups\test_joa\1overtone\realneuron_lookups_32nm_500kHz_fs0.75_1overtones_2024_04_24_11_09_41_merged_ext.pkl",prints=True)


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