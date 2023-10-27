""""this file contains a process in which all the necessary information is read out of the mod files of the Aberra/BBP cells
    that is used to create the RealisticNeuron model -> translation of NEURON to Python"""

"-----IMPORTS-----"

import numpy as np
#import time
#import datetime
import matplotlib.pyplot as plt
#import scipy.io as sio
import os
#import sys
#import math
#import seaborn as sns
#import logging
import re
#from scipy.interpolate import griddata, interpn
#nrn_options = "-NSTACK 10000 -NFRAME 525"
#Aberra recommends: -NSTACK 100000 -NFRAME 20000
#os.environ["NEURON_MODULE_OPTIONS"] = nrn_options
from neuron import h, gui
#h("NSTACK_size = 10000")
h.load_file('init.hoc')
#h.load_file('thresh4.hoc') 
#h.load_file('get_es2.hoc') #for TMS
#h.load_file('interp_coordinates.hoc') #for TMS (to calc Dx,Dy and Dz)

import tempConstants as tc
import tempFunctions as tf
import prev.functions as fs
import prev.Interp3Dfield as tt
import PySONIC as ps
import MorphoSONIC as ms

'''NEURON translator of Lemaire (2021)'''
# from PySONIC.core import PointNeuron
# from PySONIC.neurons import Cortical, OtsukaSTN, getPointNeuron, getNeuronsDict
# from MorphoSONIC.core import NmodlTranslator

# # create a NmodlTranslator from MorphoSonic (python->NMODL)
# pndict = getNeuronsDict()
# #print(pndict)
# pn = getPointNeuron('LTS')
# MODtrans = NmodlTranslator(pn)

"-----TIME AND CELL INIT-----"

#current_time = datetime.datetime.now()
#now = datetime.datetime.strftime(current_time,'%d_%m_%Y_%H_%M_%S')
#PROJECT_ROOT = os.path.dirname(os.path.realpath(__file__))


cell_nr = 7
h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
h.cell_chooser(cell_nr)
#get cell name folder based on the cell that has been chosen
# print(h.cell_names)
# for e in h.cell_names:
#     a = str(e)
#     print(a.len())
#print(h.topology()) #print this to decide the code for the cell below
if cell_nr == 2:
    cell = h.bNAC219_L1_NGCDA_e7cec642c3[0] # for cell = 2
elif cell_nr == 3:
    cell = h.bNAC219_L1_NGCDA_46b45974f4[0] # for cell = 3
elif cell_nr == 7:
    cell = h.cADpyr229_L23_PC_8ef1aa6602[0] # for cell = 7

"-----TESTS-----"


#tf.gbar_mechs()
#tf.cell_mechs()
#tf.sec_mechs()
#fs.plot_ap(cell)
# tf.dist_sec_soma()


"-----CODE-----"

"put code out of comment to give the user the option of the cell"
#cell chooser
cell_folder = "L23_PC_cADpyr229_2" #"L4_LBC_cACint209_2" #input("Give folder name of BBP cell:")
mech_folder = "cells/"+cell_folder+"/mechanisms/"

#part of the neuron -> type of section
location = 'somatic' #input("Give the type of compartment: (somatic, axonal, apical, basal)")
''' somatic = soma, axonal = axon, apical = apic, basal = dend'''
if location not in ['somatic', 'axonal', 'apical', 'basal']:
    raise ValueError("Wrong type of location, the options are given between brackets")
mechanics = tf.sec_mechs()[location]

#read all maximal channel conductances from biophysics.hoc
dist_2_soma = 0.001 
g_dict = tf.read_gbars("cells/"+cell_folder+"/",dist_2_soma) #S/m2


# put all .mod files text into 1 list -> 2D list of text, list of mod file names
mod_files, mod_names = tf.read_mod(mech_folder,mechanics)

# read out all variables and their one line (insufficient) formulas and returns where it has found these formulas
l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names)

#calculates all values of alpha, beta, tau and inf for all mechanisms given the membrane potential
# -> gating states kinetics
"put code out of comment to choose the membrane potential"
Vm = 2 #float(input("Give membrane potential Vm (mV):"))
#Vm *= 1e-3 # to convert it to V

mod_number, variables = 0, 0
gating_states_kinetics = {}
for e,f in zip(mod_files,mod_names):
    #print("mod file:",f)
    if mod_number in hits:
        variables = tf.gating_from_PROCEDURES(e,f,Vm=Vm)
    #print(variables,"\n")
    mod_number += 1
    if variables:
        gating_states_kinetics[f.replace('.mod','')] = variables

#states dictionary containing all activation and inactivation gate parameters 
# -> states names & descriptions
states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless

# calculates the steady state values of the gating variables from the gating states kinetics
# -> steady states
steadyStates = tf.steadystates_from_gating(states,gating_states_kinetics) #dimensionless

#calculates the derivative state values of the gating variables from the 4 gating parameters 
# -> states derivatives
"put code out of comment let the user initialize the opening probabilities"
x_dict = {key : 0.5 for key in states.keys()} #{key : float(input(f'Give opening probability x for {key}:')) for key in states.keys()}
derStates = tf.derstates_from_gating(states,gating_states_kinetics,x_dict) #s-1

mod_number, variables = 0, 0
currents_with_param, currents = {}, {}
for e,f in zip(mod_files,mod_names):
    #print("mod file:",f)
    if mod_number in hits:
        variables, reduced = tf.currents_from_BREAKPOINT(e,f,Vm,x_dict,g_dict,location)
    #print(variables,"\n")
    mod_number += 1
    if variables:
        currents_with_param[f.replace('.mod','')] = variables
        currents[f.replace('.mod','')] = reduced


# important prints to check with PySONIC
# print('-----------------------------------------')
# print(states)
# print('-----------------------------------------')
# print(gating_states_kinetics)
# print('-----------------------------------------')
# print(steadyStates)
# print('-----------------------------------------')
# print(derStates)
# print('-----------------------------------------')
# print(g_dict)
# print('-----------------------------------------')
# print(currents_with_param)
# print(currents)

#---------------------------------------------------------------------------------------------------------------#

'''DEPRECATED'''
''' wrong because only takes into consideration one line equations'''
'''formula dictionary containing all extracted formulas -> converts output of filter_mod into a dictionary with
key = variable and value = (one-line) formula'''
# derStates = tf.formulas_from_lists(l_alphas, l_betas, l_taus, l_infs)
# alphas, betas, taus, infs = derStates
# for e in derStates:
#     for f in e:
#         print(f,e[f])

'''this one is already wrong because it takes only oe line equations (not taking into account ifs) and creates
    a lambda function instead of just giving the value (for a given Vm)'''
#steadyStates = tf.steadystates_from_gating_old(alphas,betas,taus,infs,states)
# print(steadyStates.keys())
# print(steadyStates["m_Ca"](5))