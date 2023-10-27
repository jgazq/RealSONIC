""""just trying some things with the PySONIC module to figure out how it works"""

"-----IMPORTS-----"

import numpy as np
#import time
#import datetime
import matplotlib.pyplot as plt
#import scipy.io as sio
#import os
#import sys
#import math
#import seaborn as sns
#import logging
#from scipy.interpolate import griddata, interpn
import pandas as pd

import tempConstants as tc
import tempFunctions as tf
import prev.functions as fs
import prev.Interp3Dfield as tt
import PySONIC as ps
import MorphoSONIC as ms

"-----CODE-----"
#---create concentric planar disk transducer source
# psource_transd = ms.PlanarDiskTransducerSource(x=(0,0,0),f=0.5e6,rho=tc.RHO,c=tc.C,r=(2e-3)/2,u=1)
# grid_type = 'concentric'
# tf.plt_transdistr(psource_transd,grid_type)

#---read out metadata from lookup files (pickle files)
# pkl_txt =tf.read_pickle(r'C:\\Users\\jgazquez\\PySONIC\\PySONIC\\lookups\\',r'Cm_lkp_32nm.pkl')
# print(pkl_txt['tables']['Cm_rel'].shape)
# pkl_txt =tf.read_pickle(r'C:\\Users\\jgazquez\\PySONIC\\PySONIC\\lookups\\',r'pas_Cm0_1.0uF_cm2_ELeak_-70.0mV_lookups_32nm_20kHz.pkl')
print("RS\n\n")
pkl_txt = tf.read_pickle('c:\\users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\test_joa\\','RS_lookups_fs1.00_test.pkl')
print("realneuron\n\n")
pkl_txt = tf.read_pickle('c:\\users\\jgazquez\\RealSONIC\\PySONIC\\lookups\\test_joa\\','realneuron_lookups_fs1.00.pkl')

