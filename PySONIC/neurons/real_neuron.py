# -*- coding: utf-8 -*-
# @Author: Theo Lemaire
# @Email: theo.lemaire@epfl.ch
# @Date:   2019-06-11 15:58:38
# @Last Modified by:   Joaquin Gazquez
# @Last Modified time: 2023-12-04 11:33:54
                   
import numpy as np
from neuron import h
import sys
sys.path.append(r"C:\Users\jgazquez\RealSONIC")
import tempFunctions as tf

from ..core import PointNeuron, addSonicFeatures

@addSonicFeatures
class RealisticNeuron(PointNeuron):
    ''' Realistic neuron class '''

    # Neuron name
    name = 'realneuron'

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

    mod_files, mod_names = tf.read_mod("mechanisms/")
    g_dict = tf.read_gbars("cells/"+"L23_PC_cADpyr229_2"+"/",dist_2_soma)

    # ------------------------------ States names & descriptions ------------------------------
    states = {
		'm_Ca' : 'Ca activation gate',
		'h_Ca' : 'Ca inactivation gate',
		'm_CaHVA' : 'CaHVA activation gate',
		'h_CaHVA' : 'CaHVA inactivation gate',
		'm_Ih' : 'Ih activation gate',
		'm_Im' : 'Im activation gate',
		'm_NapEt2' : 'NapEt2 activation gate',
		'h_NapEt2' : 'NapEt2 inactivation gate',
		'm_NaTat' : 'NaTat activation gate',
		'h_NaTat' : 'NaTat inactivation gate',
		'm_NaTs2t' : 'NaTs2t activation gate',
		'h_NaTs2t' : 'NaTs2t inactivation gate',
		'm_CaLVAst' : 'CaLVAst activation gate',
		'h_CaLVAst' : 'CaLVAst inactivation gate',
		'm_KPst' : 'KPst activation gate',
		'h_KPst' : 'KPst inactivation gate',
		'm_KTst' : 'KTst activation gate',
		'h_KTst' : 'KTst inactivation gate',
		'm_SKv31' : 'SKv31 activation gate',
		'm_KdShu2007' : 'KdShu2007 activation gate',
		'h_KdShu2007' : 'KdShu2007 inactivation gate',
	}

    # ------------------------------ Gating states kinetics ------------------------------

    @classmethod
    def alpham_Ca(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[0], mod_name='Ca')
        return variables['m'+'alpha']

    @classmethod
    def betam_Ca(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[0], mod_name='Ca')
        return variables['m'+'beta']

    @classmethod
    def alphah_Ca(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[0], mod_name='Ca')
        return variables['h'+'alpha']

    @classmethod
    def betah_Ca(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[0], mod_name='Ca')
        return variables['h'+'beta']




    @classmethod
    def alpham_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaHVA')
        return variables['m'+'alpha']

    @classmethod
    def betam_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaHVA')
        return variables['m'+'beta']

    @classmethod
    def alphah_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaHVA')
        return variables['h'+'alpha']

    @classmethod
    def betah_CaHVA(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[2], mod_name='CaHVA')
        return variables['h'+'beta']




    @classmethod
    def alpham_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[3], mod_name='CaLVAst')
        return variables['m'+'alpha']

    @classmethod
    def betam_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[3], mod_name='CaLVAst')
        return variables['m'+'beta']

    @classmethod
    def alphah_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[3], mod_name='CaLVAst')
        return variables['h'+'alpha']

    @classmethod
    def betah_CaLVAst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[3], mod_name='CaLVAst')
        return variables['h'+'beta']




    @classmethod
    def alpham_Ih(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[4], mod_name='Ih')
        return variables['m'+'alpha']

    @classmethod
    def betam_Ih(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[4], mod_name='Ih')
        return variables['m'+'beta']




    @classmethod
    def alpham_Im(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[5], mod_name='Im')
        return variables['m'+'alpha']

    @classmethod
    def betam_Im(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[5], mod_name='Im')
        return variables['m'+'beta']




    @classmethod
    def alpham_KdShu2007(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KdShu2007')
        return variables['m'+'alpha']

    @classmethod
    def betam_KdShu2007(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KdShu2007')
        return variables['m'+'beta']

    @classmethod
    def alphah_KdShu2007(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KdShu2007')
        return variables['h'+'alpha']

    @classmethod
    def betah_KdShu2007(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KdShu2007')
        return variables['h'+'beta']




    @classmethod
    def alpham_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='KPst')
        return variables['m'+'alpha']

    @classmethod
    def betam_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='KPst')
        return variables['m'+'beta']

    @classmethod
    def alphah_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='KPst')
        return variables['h'+'alpha']

    @classmethod
    def betah_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='KPst')
        return variables['h'+'beta']




    @classmethod
    def alpham_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='KTst')
        return variables['m'+'alpha']

    @classmethod
    def betam_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='KTst')
        return variables['m'+'beta']

    @classmethod
    def alphah_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='KTst')
        return variables['h'+'alpha']

    @classmethod
    def betah_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='KTst')
        return variables['h'+'beta']




    @classmethod
    def alpham_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NapEt2')
        return variables['m'+'alpha']

    @classmethod
    def betam_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NapEt2')
        return variables['m'+'beta']

    @classmethod
    def alphah_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NapEt2')
        return variables['h'+'alpha']

    @classmethod
    def betah_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NapEt2')
        return variables['h'+'beta']




    @classmethod
    def alpham_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[10], mod_name='NaTat')
        return variables['m'+'alpha']

    @classmethod
    def betam_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[10], mod_name='NaTat')
        return variables['m'+'beta']

    @classmethod
    def alphah_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[10], mod_name='NaTat')
        return variables['h'+'alpha']

    @classmethod
    def betah_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[10], mod_name='NaTat')
        return variables['h'+'beta']




    @classmethod
    def alpham_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[11], mod_name='NaTs2t')
        return variables['m'+'alpha']

    @classmethod
    def betam_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[11], mod_name='NaTs2t')
        return variables['m'+'beta']

    @classmethod
    def alphah_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[11], mod_name='NaTs2t')
        return variables['h'+'alpha']

    @classmethod
    def betah_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[11], mod_name='NaTs2t')
        return variables['h'+'beta']







    @classmethod
    def alpham_SKv31(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[14], mod_name='SKv31')
        return variables['m'+'alpha']

    @classmethod
    def betam_SKv31(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[14], mod_name='SKv31')
        return variables['m'+'beta']




    # ------------------------------ States derivatives ------------------------------

    @classmethod
    def derStates(cls):
        return {
			'm_Ca' : lambda Vm, x: cls.alpham_Ca(Vm) * (1 - x['m_Ca']) - cls.betam_Ca(Vm) * x['m_Ca'],
			'h_Ca' : lambda Vm, x: cls.alphah_Ca(Vm) * (1 - x['h_Ca']) - cls.betah_Ca(Vm) * x['h_Ca'],
			'm_CaHVA' : lambda Vm, x: cls.alpham_CaHVA(Vm) * (1 - x['m_CaHVA']) - cls.betam_CaHVA(Vm) * x['m_CaHVA'],
			'h_CaHVA' : lambda Vm, x: cls.alphah_CaHVA(Vm) * (1 - x['h_CaHVA']) - cls.betah_CaHVA(Vm) * x['h_CaHVA'],
			'm_Ih' : lambda Vm, x: cls.alpham_Ih(Vm) * (1 - x['m_Ih']) - cls.betam_Ih(Vm) * x['m_Ih'],
			'm_Im' : lambda Vm, x: cls.alpham_Im(Vm) * (1 - x['m_Im']) - cls.betam_Im(Vm) * x['m_Im'],
			'm_NapEt2' : lambda Vm, x: cls.alpham_NapEt2(Vm) * (1 - x['m_NapEt2']) - cls.betam_NapEt2(Vm) * x['m_NapEt2'],
			'h_NapEt2' : lambda Vm, x: cls.alphah_NapEt2(Vm) * (1 - x['h_NapEt2']) - cls.betah_NapEt2(Vm) * x['h_NapEt2'],
			'm_NaTat' : lambda Vm, x: cls.alpham_NaTat(Vm) * (1 - x['m_NaTat']) - cls.betam_NaTat(Vm) * x['m_NaTat'],
			'h_NaTat' : lambda Vm, x: cls.alphah_NaTat(Vm) * (1 - x['h_NaTat']) - cls.betah_NaTat(Vm) * x['h_NaTat'],
			'm_NaTs2t' : lambda Vm, x: cls.alpham_NaTs2t(Vm) * (1 - x['m_NaTs2t']) - cls.betam_NaTs2t(Vm) * x['m_NaTs2t'],
			'h_NaTs2t' : lambda Vm, x: cls.alphah_NaTs2t(Vm) * (1 - x['h_NaTs2t']) - cls.betah_NaTs2t(Vm) * x['h_NaTs2t'],
			'm_CaLVAst' : lambda Vm, x: cls.alpham_CaLVAst(Vm) * (1 - x['m_CaLVAst']) - cls.betam_CaLVAst(Vm) * x['m_CaLVAst'],
			'h_CaLVAst' : lambda Vm, x: cls.alphah_CaLVAst(Vm) * (1 - x['h_CaLVAst']) - cls.betah_CaLVAst(Vm) * x['h_CaLVAst'],
			'm_KPst' : lambda Vm, x: cls.alpham_KPst(Vm) * (1 - x['m_KPst']) - cls.betam_KPst(Vm) * x['m_KPst'],
			'h_KPst' : lambda Vm, x: cls.alphah_KPst(Vm) * (1 - x['h_KPst']) - cls.betah_KPst(Vm) * x['h_KPst'],
			'm_KTst' : lambda Vm, x: cls.alpham_KTst(Vm) * (1 - x['m_KTst']) - cls.betam_KTst(Vm) * x['m_KTst'],
			'h_KTst' : lambda Vm, x: cls.alphah_KTst(Vm) * (1 - x['h_KTst']) - cls.betah_KTst(Vm) * x['h_KTst'],
			'm_SKv31' : lambda Vm, x: cls.alpham_SKv31(Vm) * (1 - x['m_SKv31']) - cls.betam_SKv31(Vm) * x['m_SKv31'],
			'm_KdShu2007' : lambda Vm, x: cls.alpham_KdShu2007(Vm) * (1 - x['m_KdShu2007']) - cls.betam_KdShu2007(Vm) * x['m_KdShu2007'],
			'h_KdShu2007' : lambda Vm, x: cls.alphah_KdShu2007(Vm) * (1 - x['h_KdShu2007']) - cls.betah_KdShu2007(Vm) * x['h_KdShu2007'],
		}

    # ------------------------------ Steady states ------------------------------

    @classmethod
    def steadyStates(cls):
        return {
			'm_Ca' : lambda Vm: cls.alpham_Ca(Vm) / (cls.alpham_Ca(Vm) + cls.betam_Ca(Vm)),
			'h_Ca' : lambda Vm: cls.alphah_Ca(Vm) / (cls.alphah_Ca(Vm) + cls.betah_Ca(Vm)),
			'm_CaHVA' : lambda Vm: cls.alpham_CaHVA(Vm) / (cls.alpham_CaHVA(Vm) + cls.betam_CaHVA(Vm)),
			'h_CaHVA' : lambda Vm: cls.alphah_CaHVA(Vm) / (cls.alphah_CaHVA(Vm) + cls.betah_CaHVA(Vm)),
			'm_Ih' : lambda Vm: cls.alpham_Ih(Vm) / (cls.alpham_Ih(Vm) + cls.betam_Ih(Vm)),
			'm_Im' : lambda Vm: cls.alpham_Im(Vm) / (cls.alpham_Im(Vm) + cls.betam_Im(Vm)),
			'm_NapEt2' : lambda Vm: cls.alpham_NapEt2(Vm) / (cls.alpham_NapEt2(Vm) + cls.betam_NapEt2(Vm)),
			'h_NapEt2' : lambda Vm: cls.alphah_NapEt2(Vm) / (cls.alphah_NapEt2(Vm) + cls.betah_NapEt2(Vm)),
			'm_NaTat' : lambda Vm: cls.alpham_NaTat(Vm) / (cls.alpham_NaTat(Vm) + cls.betam_NaTat(Vm)),
			'h_NaTat' : lambda Vm: cls.alphah_NaTat(Vm) / (cls.alphah_NaTat(Vm) + cls.betah_NaTat(Vm)),
			'm_NaTs2t' : lambda Vm: cls.alpham_NaTs2t(Vm) / (cls.alpham_NaTs2t(Vm) + cls.betam_NaTs2t(Vm)),
			'h_NaTs2t' : lambda Vm: cls.alphah_NaTs2t(Vm) / (cls.alphah_NaTs2t(Vm) + cls.betah_NaTs2t(Vm)),
			'm_CaLVAst' : lambda Vm: cls.alpham_CaLVAst(Vm) / (cls.alpham_CaLVAst(Vm) + cls.betam_CaLVAst(Vm)),
			'h_CaLVAst' : lambda Vm: cls.alphah_CaLVAst(Vm) / (cls.alphah_CaLVAst(Vm) + cls.betah_CaLVAst(Vm)),
			'm_KPst' : lambda Vm: cls.alpham_KPst(Vm) / (cls.alpham_KPst(Vm) + cls.betam_KPst(Vm)),
			'h_KPst' : lambda Vm: cls.alphah_KPst(Vm) / (cls.alphah_KPst(Vm) + cls.betah_KPst(Vm)),
			'm_KTst' : lambda Vm: cls.alpham_KTst(Vm) / (cls.alpham_KTst(Vm) + cls.betam_KTst(Vm)),
			'h_KTst' : lambda Vm: cls.alphah_KTst(Vm) / (cls.alphah_KTst(Vm) + cls.betah_KTst(Vm)),
			'm_SKv31' : lambda Vm: cls.alpham_SKv31(Vm) / (cls.alpham_SKv31(Vm) + cls.betam_SKv31(Vm)),
			'm_KdShu2007' : lambda Vm: cls.alpham_KdShu2007(Vm) / (cls.alpham_KdShu2007(Vm) + cls.betam_KdShu2007(Vm)),
			'h_KdShu2007' : lambda Vm: cls.alphah_KdShu2007(Vm) / (cls.alphah_KdShu2007(Vm) + cls.betah_KdShu2007(Vm)),
		}

    # ------------------------------ Membrane currents ------------------------------
    @classmethod
    def i_Ca(cls,m_Ca,h_Ca,Vm):
        ''' iCa current '''
        x_dict = {'e': e for e in [m_Ca,h_Ca]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[0], mod_name='Ca', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_CaHVA(cls,m_CaHVA,h_CaHVA,Vm):
        ''' iCaHVA current '''
        x_dict = {'e': e for e in [m_CaHVA,h_CaHVA]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[2], mod_name='CaHVA', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
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
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[3], mod_name='CaLVAst', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_Ih(cls,m_Ih,Vm):
        ''' iIh current '''
        x_dict = {'e': e for e in [m_Ih]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[4], mod_name='Ih', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_Im(cls,m_Im,Vm):
        ''' iIm current '''
        x_dict = {'e': e for e in [m_Im]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[5], mod_name='Im', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_KdShu2007(cls,m_KdShu2007,h_KdShu2007,Vm):
        ''' iKdShu2007 current '''
        x_dict = {'e': e for e in [m_KdShu2007,h_KdShu2007]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[6], mod_name='KdShu2007', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_KPst(cls,m_KPst,h_KPst,Vm):
        ''' iKPst current '''
        x_dict = {'e': e for e in [m_KPst,h_KPst]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[7], mod_name='KPst', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_KTst(cls,m_KTst,h_KTst,Vm):
        ''' iKTst current '''
        x_dict = {'e': e for e in [m_KTst,h_KTst]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[8], mod_name='KTst', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_NapEt2(cls,m_NapEt2,h_NapEt2,Vm):
        ''' iNapEt2 current '''
        x_dict = {'e': e for e in [m_NapEt2,h_NapEt2]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[9], mod_name='NapEt2', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_NaTat(cls,m_NaTat,h_NaTat,Vm):
        ''' iNaTat current '''
        x_dict = {'e': e for e in [m_NaTat,h_NaTat]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[10], mod_name='NaTat', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_NaTs2t(cls,m_NaTs2t,h_NaTs2t,Vm):
        ''' iNaTs2t current '''
        x_dict = {'e': e for e in [m_NaTs2t,h_NaTs2t]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[11], mod_name='NaTs2t', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def i_SKv31(cls,m_SKv31,Vm):
        ''' iSKv31 current '''
        x_dict = {'e': e for e in [m_SKv31]}
        variables = tf.currents_from_BREAKPOINT(list_mod=cls.mod_files[14], mod_name='SKv31', Vm=Vm, x_dict = x_dict, g_dict = cls.g_dict, location = "somatic")
        currents = [e for e in variables.keys() if (e.startswith('i') or e.startswith('I'))]
        print(currents)
        if currents:
            return variables[currents[0]]
        else:
            return 0

    @classmethod
    def currents(cls):
        return {
			'i_Ca': lambda Vm, x: cls.i_Ca(x['m_Ca'], x['h_Ca'], Vm),
			'i_CaHVA': lambda Vm, x: cls.i_CaHVA(x['m_CaHVA'], x['h_CaHVA'], Vm),
			'i_CaLVAst': lambda Vm, x: cls.i_CaLVAst(x['m_CaLVAst'], x['h_CaLVAst'], Vm),
			'i_Ih': lambda Vm, x: cls.i_Ih(x['m_Ih'], Vm),
			'i_Im': lambda Vm, x: cls.i_Im(x['m_Im'], Vm),
			'i_KdShu2007': lambda Vm, x: cls.i_KdShu2007(x['m_KdShu2007'], x['h_KdShu2007'], Vm),
			'i_KPst': lambda Vm, x: cls.i_KPst(x['m_KPst'], x['h_KPst'], Vm),
			'i_KTst': lambda Vm, x: cls.i_KTst(x['m_KTst'], x['h_KTst'], Vm),
			'i_NapEt2': lambda Vm, x: cls.i_NapEt2(x['m_NapEt2'], x['h_NapEt2'], Vm),
			'i_NaTat': lambda Vm, x: cls.i_NaTat(x['m_NaTat'], x['h_NaTat'], Vm),
			'i_NaTs2t': lambda Vm, x: cls.i_NaTs2t(x['m_NaTs2t'], x['h_NaTs2t'], Vm),
			'i_SKv31': lambda Vm, x: cls.i_SKv31(x['m_SKv31'], Vm),
        }