""""just testing some things"""

import numpy as np
import matplotlib.pyplot as plt

from neuron import h
import tempFunctions as tf
import tempConstants as tc

"""padding a LUT so it has the same dimensions in order to merge"""
# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\realneuron_lookups_16nm_3000kHz_[0_-1]_fs0.75_2024_03_27_17_40_44.pkl')
# for e in pkldict['refs']:
#     print(e)
# pkldict['refs']['Q'] = np.pad(pkldict['refs']['Q'],(0,1))
# for e in pkldict['tables']:
#     pkldict['tables'][e] = np.pad(pkldict['tables'][e],((0,0),(0,0),(0,0),(0,1),(0,0),(0,0)))
#     print(pkldict['tables'][e].shape)   
# tf.load_pickle(pkldict,r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\realneuron_lookups_16nm_3000kHz_[0_-1]_fs0.75_2024_03_27_17_40_44_padded.pkl')


"""to plot all gating parameters of a LUT"""
# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\0overtones\realneuron_lookups_merged_ext.pkl')
# tf.save_gatingplots(pkldict,r'test\radius',a='all',reduced_xrange=False)

""""look up the value in the table"""
# pkldict = tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_32nm_500kHz_fs0.75.pkl')
# for e,f in pkldict['refs'].items():
#     print(e,f[265:275]) if e == 'Q' else None
# a=32*1e-9
# f=500*1e3
# A=101.64826607291788*1e3
# Cm0=0.01
# Q = 50*1e-5
# var_dict = {'a': a, 'f': f, 'A': A, 'Q': Q, 'Cm0': Cm0}
# ind_list = [int(np.where(abs((pkldict['refs'][k]-v)/(v+1e-10))<0.01)[0][0]) for (k,v) in var_dict.items()]
# ind_list = f'[{ind_list[0]}, {ind_list[1]}, {ind_list[2]}, {ind_list[3]}, {ind_list[4]}, {0}]'
# print(f"pkldict['tables']['alpham_CaHVA']{ind_list}")
# value = eval(f"pkldict['tables']['alpham_CaHVA']{ind_list}")
# print(value)

tf.downsample_LUT(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_fs1.00.pkl',[3,5,10,10,1,1])
#tf.read_pickle(r'C:\Users\jgazquez\PySONIC\PySONIC\lookups\test_joa\realneuron_lookups_16nm_100kHz_fs0.75_works.pkl',prints=True)