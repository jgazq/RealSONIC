# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-06-11 15:58:38
# @Last Modified by:   Joaquin Gazquez
# @Last Modified time: 2023-24-10
                   
import numpy as np
from neuron import h
import sys
sys.path.append("C:\\Users\\jgazquez\\RealSONIC")
import tempFunctions as tf

from ..core import PointNeuron, addSonicFeatures

@addSonicFeatures
class CaNeuron(PointNeuron):
    '''Only calcium neuron class '''

    # Neuron name
    name = 'ca_only'

    # ------------------------------ Biophysical parameters ------------------------------

    # Resting parameters
    Cm0 = 1e-2   # Membrane capacitance (F/m2)
    Vm0 = -71.9  # Membrane potential (mV)

         
    # Reversal potentials (mV)
    #TODO
         
    # Maximal channel conductances (S/m2)
    gskv3_1bar = 1025.17
    gsk_e2bar = 994.3299999999999
    gca_hvabar = 3.7399999999999998
    gnats2_tbar = 9267.05
    gihbar = 0.8
    gca_lvastbar = 7.78
         
    # Additional parameters
    VT = -56.2  # Spike threshold adjustment parameter (mV)
    dist_2_soma = 20 # Distance from the considered segment to the soma (um?)

    mod_files, mod_names = tf.read_mod("cells/L23_PC_cADpyr229_2/mechanisms/")
    g_dict = tf.read_gbars("cells/"+"L23_PC_cADpyr229_2"+"/",dist_2_soma)

    # ------------------------------ States names & descriptions ------------------------------
    states = {
		'm_CaHVA' : 'CaHVA activation gate',
		'h_CaHVA' : 'CaHVA inactivation gate',
		'm_CaLVAst' : 'CaLVAst activation gate',
		'h_CaLVAst' : 'CaLVAst inactivation gate',
	}

    # ------------------------------ Gating states kinetics ------------------------------

    @classmethod
    def alpham_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[1], mod_name='CaHVA')
        return variables['m'+'alpha']

    @classmethod
    def betam_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[1], mod_name='CaHVA')
        return variables['m'+'beta']

    @classmethod
    def alphah_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[1], mod_name='CaHVA')
        return variables['h'+'alpha']

    @classmethod
    def betah_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[1], mod_name='CaHVA')
        return variables['h'+'beta']




    @classmethod
    def alpham_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaLVAst')
        return variables['m'+'alpha']

    @classmethod
    def betam_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaLVAst')
        return variables['m'+'beta']

    @classmethod
    def alphah_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaLVAst')
        return variables['h'+'alpha']

    @classmethod
    def betah_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaLVAst')
        return variables['h'+'beta']


    # ------------------------------ States derivatives ------------------------------

    @classmethod
    def derStates(cls):
        return {
			'm_CaHVA' : lambda Vm, x: cls.alpham_CaHVA(Vm) * (1 - x['m_CaHVA']) - cls.betam_CaHVA(Vm) * x['m_CaHVA'],
			'h_CaHVA' : lambda Vm, x: cls.alphah_CaHVA(Vm) * (1 - x['h_CaHVA']) - cls.betah_CaHVA(Vm) * x['h_CaHVA'],
			'm_CaLVAst' : lambda Vm, x: cls.alpham_CaLVAst(Vm) * (1 - x['m_CaLVAst']) - cls.betam_CaLVAst(Vm) * x['m_CaLVAst'],
			'h_CaLVAst' : lambda Vm, x: cls.alphah_CaLVAst(Vm) * (1 - x['h_CaLVAst']) - cls.betah_CaLVAst(Vm) * x['h_CaLVAst'],
		}

    # ------------------------------ Steady states ------------------------------

    @classmethod
    def steadyStates(cls):
        return {
			'm_CaHVA' : lambda Vm: cls.alpham_CaHVA(Vm) / (cls.alpham_CaHVA(Vm) + cls.betam_CaHVA(Vm)),
			'h_CaHVA' : lambda Vm: cls.alphah_CaHVA(Vm) / (cls.alphah_CaHVA(Vm) + cls.betah_CaHVA(Vm)),
			'm_CaLVAst' : lambda Vm: cls.alpham_CaLVAst(Vm) / (cls.alpham_CaLVAst(Vm) + cls.betam_CaLVAst(Vm)),
			'h_CaLVAst' : lambda Vm: cls.alphah_CaLVAst(Vm) / (cls.alphah_CaLVAst(Vm) + cls.betah_CaLVAst(Vm)),
		}

    # ------------------------------ Membrane currents ------------------------------
    @classmethod
    def i_CaHVA(cls,m_CaHVA,h_CaHVA,Vm):
        ''' iCaHVA current '''
        x_dict = {'e': e for e in [m_CaHVA,h_CaHVA]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[1], mod_name='CaHVA', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_CaLVAst(cls,m_CaLVAst,h_CaLVAst,Vm):
        ''' iCaLVAst current '''
        x_dict = {'e': e for e in [m_CaLVAst,h_CaLVAst]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[2], mod_name='CaLVAst', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def currents(cls):
        return {
			'i_CaHVA': lambda Vm, x: cls.i_CaHVA(x['m_CaHVA'], x['h_CaHVA'], Vm),
			'i_CaLVAst': lambda Vm, x: cls.i_CaLVAst(x['m_CaLVAst'], x['h_CaLVAst'], Vm),
        }