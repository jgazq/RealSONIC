"""" this file writes all the necessary code to PySONIC/PySONIC/neurons/real_neuron.py and RealSONIC/PySONIC/neurons/real_neuron.py to create the RealisticNeuron class
    in order to create the LUT for all the different mechanisms"""


import numpy as np
from neuron import h
import shutil
import sys
import os
#sys.path.append("C:\\Users\\jgazquez\\RealSONIC")
import tempFunctions as tf

""""--------------- INPUTS ---------------"""
cell_nr = 7
cell_folder = "L23_PC_cADpyr229_2"
sec_type = 'somatic'
dist_2_soma = 20 #um

""""--------------------------------------"""

mech_folder = "cells/"+cell_folder+"/mechanisms/"
mod_files, mod_names = tf.read_mod(mech_folder)

for i,e in enumerate(mod_names):
    mod_names[i] = tf.rm_us(e)

l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names)
states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless
g_dict = tf.read_gbars("cells/"+cell_folder+"/",dist_2_soma) #S/m2

path = os. getcwd() + "\\PySONIC\\neurons\\real_neuron.py"
with open(path,'w') as filenaam:

# write header of the file and imports
    path2 = os.getcwd().replace('\\','\\\\')
    filenaam.write(f"""# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-06-11 15:58:38
# @Last Modified by:   Joaquin Gazquez
# @Last Modified time: 2023-24-10
                   
import numpy as np
from neuron import h
import sys
sys.path.append("{path2}")
import tempFunctions as tf

from ..core import PointNeuron, addSonicFeatures\n\n""")
    
# write the first part of the class

    filenaam.write(f"""@addSonicFeatures
class RealisticNeuron(PointNeuron):
    ''' Realistic neuron class '''

    # Neuron name
    name = 'realneuron'

    # ------------------------------ Biophysical parameters ------------------------------

    # Resting parameters
    Cm0 = 1e-2   # Membrane capacitance (F/m2)
    Vm0 = -71.9  # Membrane potential (mV)\n\n""")
    filenaam.write(f"""         
    # Reversal potentials (mV)
    #TODO\n""")
    filenaam.write(f"""         
    # Maximal channel conductances (S/m2)\n""")
    for e,f in g_dict[sec_type].items():
        filenaam.write(f"    {e} = {f}\n")
    filenaam.write(f"""         
    # Additional parameters
    VT = -56.2  # Spike threshold adjustment parameter (mV)
    dist_2_soma = {dist_2_soma} # Distance from the considered segment to the soma (um?)\n\n""")
    filenaam.write(f"    mod_files, mod_names = tf.read_mod(\"{mech_folder}\")\n")
    filenaam.write(f"    g_dict = tf.read_gbars(\"cells/\"+\"{cell_folder}\"+\"/\",dist_2_soma)\n\n")

# write the states names & descriptions

    filenaam.write("""    # ------------------------------ States names & descriptions ------------------------------
    states = {\n""")
    for e,f in zip(states.keys(), states.values()):
        filenaam.write(f"\t\t'{e}' : '{f}',\n")
    filenaam.write("\t}\n\n")

# write the gating states kinetics

    filenaam.write("    # ------------------------------ Gating states kinetics ------------------------------\n\n")
    mod_number, variables = 0, 0
    for e,f in zip(mod_files,mod_names):
        if mod_number in hits:
            name_mod = f.replace('.mod','')
            #filenaam.write(f"#{name_mod}")
            for x in ['m','h']:
                gating_param =  x+'_'+name_mod
                if gating_param in states:
                    filenaam.write(f"""    @classmethod
    def alpha{gating_param}(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[{mod_number}], mod_name='{f}')
        return variables['{x}'+'alpha']\n\n""")

                    filenaam.write(f"""    @classmethod
    def beta{gating_param}(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[{mod_number}], mod_name='{f}')
        return variables['{x}'+'beta']\n\n""")
            filenaam.write("\n\n\n")
        mod_number += 1 
    
# write the gating states derivatives

    filenaam.write("""    # ------------------------------ States derivatives ------------------------------

    @classmethod
    def derStates(cls):
        return {\n""")
    for e in states:
        filenaam.write(f"\t\t\t'{e}' : lambda Vm, x: cls.alpha{e}(Vm) * (1 - x['{e}']) - cls.beta{e}(Vm) * x['{e}'],\n")
    filenaam.write("\t\t}\n\n")

# write the steady states

    filenaam.write("""    # ------------------------------ Steady states ------------------------------

    @classmethod
    def steadyStates(cls):
        return {\n""")
    for e in states:
        filenaam.write(f"\t\t\t'{e}' : lambda Vm: cls.alpha{e}(Vm) / (cls.alpha{e}(Vm) + cls.beta{e}(Vm)),\n")
    filenaam.write("\t\t}\n\n")

# write the membrane currents

    filenaam.write("""    # ------------------------------ Membrane currents ------------------------------\n""")

    currents_dict = ''
    mod_number, variables = 0, 0
    for e,f in zip(mod_files,mod_names):
        if mod_number in hits:
            name_mod = f.replace('.mod','')
            gating_var = ''
            current_var = ''
            x_actual = []
            for x in ['m','h']:
                gating_param =  x+'_'+name_mod
                if gating_param in states:
                    gating_var += f'{gating_param},'
                    current_var += f'x[\'{gating_param}\'], '
                    x_actual.append(x)
            if gating_var != '':
                #code between quotation marks is original code but has been replaced with single line code
                """return cls.g_{name_mod}*Vm \n\n"""
                filenaam.write(f"""    @classmethod
    def i_{name_mod}(cls,{gating_var}Vm):
        ''' i{name_mod} current '''
        x_dict = {{'e': e for e in [{gating_var[:-1]}]}}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[{mod_number}], mod_name='{f}', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "{sec_type}")
        currents = [e for e in variables.keys() if (e.startswith(\'i\') or e.startswith(\'I\'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0\n\n""")
                currents_dict += f'\n\t\t\t\'i_{name_mod}\': lambda Vm, x: cls.i_{name_mod}({current_var}Vm),'
            
        mod_number += 1 
    filenaam.write("""    @classmethod
    def currents(cls):
        return {""")
    filenaam.write(currents_dict)
    filenaam.write('\n        }')

# copy the file from RealSONIC/PySONIC to PySONIC/PySONIC so both contain the correct real_neuron.py file
shutil.copy(path, os. getcwd().replace('RealSONIC','PySONIC') + "\\PySONIC\\neurons\\real_neuron.py")
