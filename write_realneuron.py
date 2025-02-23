"""" this file writes all the necessary code to PySONIC/PySONIC/neurons/real_neuron.py and RealSONIC/PySONIC/neurons/real_neuron.py to create the RealisticNeuron class
    in order to create the LUT for all the different mechanisms"""

import os
import numpy as np
from neuron import h
import shutil
import sys
import datetime
#sys.path.append("C:\\Users\\jgazquez\\RealSONIC")
import tempFunctions as tf
h.load_file("init.hoc")

""""--------------- INPUTS ---------------"""
cell_nr = 7
sec_type = 'somatic'
dist_2_soma = 20 #um

""""--------------------------------------"""

cell_folder = h.cell_names[cell_nr-1].s #OR h.cell_names.o(cell_nr-1).s        (cell_folder = "L23_PC_cADpyr229_2")
mech_folder = "mechanisms/"#"cells/"+cell_folder+"/mechanisms/"
mod_files, mod_names = tf.read_mod(mech_folder)

for i,e in enumerate(mod_names):
    mod_names[i] = tf.rm_us(e)

l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names)
states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless
g_dict = tf.read_gbars("cells/"+cell_folder+"/",dist_2_soma) #S/m2
bio_dict = tf.read_biophysics("cells/"+cell_folder+"/",dist_2_soma)

current_time = datetime.datetime.now()
now = datetime.datetime.strftime(current_time,'%Y-%m-%d %H:%M:%S')

this_path = os.getcwd()
print(this_path)
path = os.path.abspath(os.path.join(this_path, '..', "PySONIC/PySONIC/neurons/real_neuron.py"))
print(path)
#path = this_path + "/PySONIC/neurons/real_neuron.py"
with open(path,'w') as filenaam:

# write header of the file and imports
    #path2 = os.getcwd().replace('\\','\\\\')
    filenaam.write(f"""# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-06-11 15:58:38
# @Last Modified by:   Joaquin Gazquez
# @Last Modified time: {now}
                   
import numpy as np
from neuron import h
import sys
sys.path.append(r"{this_path}")
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
    Vm0 = -75  # Membrane potential (mV)\n\n""")
    # filenaam.write(f"""         
    # # Reversal potentials (mV)
    # #TODO\n""")
    # filenaam.write(f"""         
    # # Maximal channel conductances (S/m2)\n""")
    # # for e,f in g_dict[sec_type].items(): #this is skipped as the g_bar differs per section type and distance to the soma
    # #     filenaam.write(f"    {e} = {f}\n")
    # filenaam.write(f"""         
    # # Additional parameters
    # VT = -56.2  # Spike threshold adjustment parameter (mV)
    # #dist_2_soma = {dist_2_soma} # Distance from the considered segment to the soma (um?)\n\n""") #skipped
    #filenaam.write(f"    mod_files, mod_names = tf.read_mod(\"{mech_folder}\")\n")
    #filenaam.write(f"    #g_dict = tf.read_gbars(\"cells/\"+\"{cell_folder}\"+\"/\",dist_2_soma)\n\n") #also skipped (see comment above)

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
    def alpha{gating_param}(cls,Vm):\n""")
                    for e in tf.gating_from_PROCEDURES_list(list_mod=mod_files[mod_number], mod_name=f):
                        filenaam.write(f"""{e}\n""")
                    filenaam.write(f"""
        return {x}alpha\n\n""")

                    filenaam.write(f"""    @classmethod
    def beta{gating_param}(cls,Vm):\n""")
                    for e in tf.gating_from_PROCEDURES_list(list_mod=mod_files[mod_number], mod_name=f):
                        filenaam.write(f"""{e}\n""")
                    filenaam.write(f"""
        return {x}beta\n\n""")
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
    conductances_dict = ''
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
                curr_lines = tf.currents_from_BREAKPOINT_list(list_mod=mod_files[mod_number], mod_name=f, gating_var=x_actual) #lines of code that are written to the electric current method
                cond_lines = tf.conductances_from_BREAKPOINT_list(list_mod=mod_files[mod_number], mod_name=f, gating_var=x_actual) #lines of code that are written to the electric current method
                gbar = ''
                for i,line in enumerate(curr_lines):
                    if 'g' in line and 'bar' in line:
                        gbar = curr_lines.pop(i).strip() #remove the gbar line as this will be added to the input arguments of the method instead of putting it in the file
                        break

                filenaam.write(f"""    @classmethod
    def i_{name_mod}(cls,{gating_var}Vm,{gbar}):
        ''' i{name_mod} current '''\n""")
                for e in curr_lines[:-1]:
                    filenaam.write(f"""{e}\n""")
                filenaam.write(f'\n{curr_lines[-1]}\n\n')
                currents_dict += f'\n\t\t\t\'i_{name_mod}\': lambda Vm, x, g_bar: cls.i_{name_mod}({current_var}Vm) if g_bar is None else cls.i_{name_mod}({current_var}Vm, g_bar),'

                filenaam.write(f"""    @classmethod
    def g_{name_mod}(cls,{gating_var}Vm,{gbar}):
        ''' g{name_mod} conductance '''\n""")
                for e in cond_lines[:-1]:
                    filenaam.write(f"""{e}\n""")
                filenaam.write(f'\n{cond_lines[-1]}\n\n')
                conductances_dict += f'\n\t\t\t\'g_{name_mod}\': lambda Vm, x, g_bar: cls.g_{name_mod}({current_var}Vm) if g_bar is None else cls.g_{name_mod}({current_var}Vm, g_bar),'
            
        mod_number += 1 
    filenaam.write(f"""    @classmethod
    def i_pas(cls, Vm):
        ''' ipas current '''\n""")
    for e in bio_dict['all']:
        if str(e) == 'g_pas':
            filenaam.write(f"""        {e} = {bio_dict['all'][e]} * 1e4 #S/cm2 -> S/m2 \n""")
        else:
            filenaam.write(f"""        {e} = {bio_dict['all'][e]}\n""")
    filenaam.write(f"""        ipas = g_pas*(Vm-e_pas)\n""")
    filenaam.write(f"""\n        return ipas\n\n""")
    filenaam.write(f"""    @classmethod
    def g_pas(cls, Vm):
        ''' gpas conductance '''\n""")
    for e in bio_dict['all']:
        if str(e) == 'g_pas':
            filenaam.write(f"""        {e} = {bio_dict['all'][e]} * 1e4 #S/cm2 -> S/m2 \n""")
        else:
            filenaam.write(f"""        {e} = {bio_dict['all'][e]}\n""")
    filenaam.write(f"""\n        return g_pas\n\n""")
    currents_dict += f'\n\t\t\t\'i_pas\': lambda Vm, x, g_bar: cls.i_pas(Vm)' #added dummy argument variables for uniformity
    conductances_dict += f'\n\t\t\t\'g_pas\': lambda Vm, x, g_bar: cls.g_pas(Vm)' #added dummy argument variables for uniformity

    filenaam.write("""    @classmethod
    def currents(cls):
        return {""")
    filenaam.write(currents_dict)
    filenaam.write('\n        }\n\n')

    filenaam.write("""    @classmethod
    def conductances(cls):
        return {""")
    filenaam.write(conductances_dict)
    filenaam.write('\n        }')    

# copy the file from RealSONIC/PySONIC to PySONIC/PySONIC so both contain the correct real_neuron.py file
#shutil.copy(path, os. getcwd().replace('RealSONIC','PySONIC') + "\\PySONIC\\neurons\\real_neuron.py")