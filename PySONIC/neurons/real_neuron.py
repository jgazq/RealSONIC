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

    mod_files, mod_names = tf.read_mod("cells/L23_PC_cADpyr229_2/mechanisms/")
    g_dict = tf.read_gbars("cells/"+"L23_PC_cADpyr229_2"+"/",dist_2_soma)

    # ------------------------------ States names & descriptions ------------------------------
    states = {
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




    @classmethod
    def alpham_Ih(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[3], mod_name='Ih')
        return variables['m'+'alpha']

    @classmethod
    def betam_Ih(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[3], mod_name='Ih')
        return variables['m'+'beta']




    @classmethod
    def alpham_Im(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[4], mod_name='Im')
        return variables['m'+'alpha']

    @classmethod
    def betam_Im(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[4], mod_name='Im')
        return variables['m'+'beta']




    @classmethod
    def alpham_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[5], mod_name='KPst')
        return variables['m'+'alpha']

    @classmethod
    def betam_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[5], mod_name='KPst')
        return variables['m'+'beta']

    @classmethod
    def alphah_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[5], mod_name='KPst')
        return variables['h'+'alpha']

    @classmethod
    def betah_KPst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[5], mod_name='KPst')
        return variables['h'+'beta']




    @classmethod
    def alpham_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KTst')
        return variables['m'+'alpha']

    @classmethod
    def betam_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KTst')
        return variables['m'+'beta']

    @classmethod
    def alphah_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KTst')
        return variables['h'+'alpha']

    @classmethod
    def betah_KTst(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[6], mod_name='KTst')
        return variables['h'+'beta']




    @classmethod
    def alpham_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='NapEt2')
        return variables['m'+'alpha']

    @classmethod
    def betam_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='NapEt2')
        return variables['m'+'beta']

    @classmethod
    def alphah_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='NapEt2')
        return variables['h'+'alpha']

    @classmethod
    def betah_NapEt2(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[7], mod_name='NapEt2')
        return variables['h'+'beta']




    @classmethod
    def alpham_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='NaTat')
        return variables['m'+'alpha']

    @classmethod
    def betam_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='NaTat')
        return variables['m'+'beta']

    @classmethod
    def alphah_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='NaTat')
        return variables['h'+'alpha']

    @classmethod
    def betah_NaTat(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[8], mod_name='NaTat')
        return variables['h'+'beta']




    @classmethod
    def alpham_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NaTs2t')
        return variables['m'+'alpha']

    @classmethod
    def betam_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NaTs2t')
        return variables['m'+'beta']

    @classmethod
    def alphah_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NaTs2t')
        return variables['h'+'alpha']

    @classmethod
    def betah_NaTs2t(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[9], mod_name='NaTs2t')
        return variables['h'+'beta']







    @classmethod
    def alpham_SKv31(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[12], mod_name='SKv31')
        return variables['m'+'alpha']

    @classmethod
    def betam_SKv31(cls,Vm):
        variables = tf.gating_from_PROCEDURES(Vm=Vm, list_mod=cls.mod_files[12], mod_name='SKv31')
        return variables['m'+'beta']




    # ------------------------------ States derivatives ------------------------------

    @classmethod
    def derStates(cls):
        return {
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
		}

    # ------------------------------ Steady states ------------------------------

    @classmethod
    def steadyStates(cls):
        return {
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
		}

    # ------------------------------ Membrane currents ------------------------------
    @classmethod
    def i_CaHVA(cls,m_CaHVA,h_CaHVA,Vm):
        ''' iCaHVA current '''
        return cls.g_CaHVA*Vm 

    @classmethod
    def i_CaLVAst(cls,m_CaLVAst,h_CaLVAst,Vm):
        ''' iCaLVAst current '''
        return cls.g_CaLVAst*Vm 

    @classmethod
    def i_Ih(cls,m_Ih,Vm):
        ''' iIh current '''
        return cls.g_Ih*Vm 

    @classmethod
    def i_Im(cls,m_Im,Vm):
        ''' iIm current '''
        return cls.g_Im*Vm 

    @classmethod
    def i_KPst(cls,m_KPst,h_KPst,Vm):
        ''' iKPst current '''
        return cls.g_KPst*Vm 

    @classmethod
    def i_KTst(cls,m_KTst,h_KTst,Vm):
        ''' iKTst current '''
        return cls.g_KTst*Vm 

    @classmethod
    def i_NapEt2(cls,m_NapEt2,h_NapEt2,Vm):
        ''' iNapEt2 current '''
        return cls.g_NapEt2*Vm 

    @classmethod
    def i_NaTat(cls,m_NaTat,h_NaTat,Vm):
        ''' iNaTat current '''
        return cls.g_NaTat*Vm 

    @classmethod
    def i_NaTs2t(cls,m_NaTs2t,h_NaTs2t,Vm):
        ''' iNaTs2t current '''
        return cls.g_NaTs2t*Vm 

    @classmethod
    def i_SKv31(cls,m_SKv31,Vm):
        ''' iSKv31 current '''
        return cls.g_SKv31*Vm 

    @classmethod
    def currents(cls):
        return {
			'i_CaHVA': lambda Vm, x: cls.i_CaHVA(x['m_CaHVA'], x['h_CaHVA'], Vm),
			'i_CaLVAst': lambda Vm, x: cls.i_CaLVAst(x['m_CaLVAst'], x['h_CaLVAst'], Vm),
			'i_Ih': lambda Vm, x: cls.i_Ih(x['m_Ih'], Vm),
			'i_Im': lambda Vm, x: cls.i_Im(x['m_Im'], Vm),
			'i_KPst': lambda Vm, x: cls.i_KPst(x['m_KPst'], x['h_KPst'], Vm),
			'i_KTst': lambda Vm, x: cls.i_KTst(x['m_KTst'], x['h_KTst'], Vm),
			'i_NapEt2': lambda Vm, x: cls.i_NapEt2(x['m_NapEt2'], x['h_NapEt2'], Vm),
			'i_NaTat': lambda Vm, x: cls.i_NaTat(x['m_NaTat'], x['h_NaTat'], Vm),
			'i_NaTs2t': lambda Vm, x: cls.i_NaTs2t(x['m_NaTs2t'], x['h_NaTs2t'], Vm),
			'i_SKv31': lambda Vm, x: cls.i_SKv31(x['m_SKv31'], Vm),
        }














# # -*- coding: utf-8 -*-
# # @Author: Theo Lemaire
# # @Email: theo.lemaire@epfl.ch
# # @Date:   2019-06-11 15:58:38
# # @Last Modified by:   Theo Lemaire
# # @Last Modified time: 2020-03-31 18:14:08

# import numpy as np
# from neuron import h
# import sys
# sys.path.append("C:\\Users\\jgazquez\\RealSONIC")
# import tempFunctions as tf

# from ..core import PointNeuron, addSonicFeatures


# # cell_nr = 7
# # h.setParamsAdultHuman() #this needs to go before the cell chooser, otherwise it won't make a difference
# # h.cell_chooser(cell_nr)

# cell_folder = "L23_PC_cADpyr229_2" #"L4_LBC_cACint209_2" #input("Give folder name of BBP cell:")
# mech_folder = "cells/"+cell_folder+"/mechanisms/"
# mod_files, mod_names = tf.read_mod(mech_folder)
# l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names)
# states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless


# @addSonicFeatures
# class RealisticNeuron(PointNeuron):
#     ''' Template neuron class '''

#     # Neuron name
#     name = 'realneuron'

#     # ------------------------------ Biophysical parameters ------------------------------

#     # Resting parameters
#     Cm0 = 1e-2   # Membrane capacitance (F/m2)
#     Vm0 = -71.9  # Membrane potential (mV)

#     # Reversal potentials (mV)
#     # ENa = 50.0     # Sodium
#     # EK = -85.0 #-90.0     # Potassium
#     # ELeak = -70.3  # Non-specific leakage

#     # Maximal channel conductances (S/m2)
#     # gNabar = 560.0  # Sodium
#     # gKdbar = 60.0   # Delayed-rectifier Potassium
#     # gLeak = 0.205   # Non-specific leakage

#     # Additional parameters
#     VT = -56.2  # Spike threshold adjustment parameter (mV)

#     # Aberra mechanisms
#     # gIhbar = 8e-05*1e4 #0.00001*1e4          #Ih (S/cm2 -> S/m2)
#     # ehcn =  -45.0                            #Ih (mV), hyperpolarization-activated 
#     #                                          #cyclic nucleotide-gated channels
#     # gImbar = 0.00074*1e4 #0.00001*1e4        #Im (S/cm2 -> S/m2)
#     # gKPstbar = 0.959296*1e4 #0.00001*1e4     #K_Pst (S/cm2 -> S/m2)
#     # gKTstbar = 0.001035*1e4 #0.00001*1e4     #K_Tst (S/cm2 -> S/m2)
#     # gNapEt2bar = 0.009803*1e4 #0.00001*1e4   #Nap_Et2 (S/cm2 -> S/m2)
#     # gNaTatbar = 3.429725*1e4 #0.00001*1e4    #NaTa_t (S/cm2 -> S/m2)
#     # gNaTs2tbar = 0.926705*1e4 #0.00001*1e4   #NaTs2_t (S/cm2 -> S/m2)
#     "SK are calcium activated so are not taken into consideration"
#     #v =                         #SK_E -> maybe inside this __new__ method?
#     #gSK_E2bar = .000001*1e4     #SK_E (S/cm2 -> S/m2)
#     #zTau = 1                    #SK_E (ms)
#     #ek =                        #SK_E (mV)
#     #cai =                       #SK_E (mM)
#     #gSKv3_1bar = 0.00001*1e4    #SKv3_1 (S/cm2 -> S/m2)

#     VIh = 154.9
#     T_C = 34 # Temperature (°C)
#     gating_kinetics = {}
#     dummy = 1

#     # ------------------------------ States names & descriptions ------------------------------
#     cell_folder = "L23_PC_cADpyr229_2" #"L4_LBC_cACint209_2" #input("Give folder name of BBP cell:")
#     mech_folder = "cells/"+cell_folder+"/mechanisms/"
#     mod_files, mod_names = tf.read_mod(mech_folder)
#     l_alphas, l_betas, l_taus, l_infs, hits = tf.filter_mod(mod_files,mod_names)
#     states = tf.states_from_lists(l_alphas, l_betas, l_taus, l_infs) #dimensionless
#     states = {'m_Ih' : states['m_Ih']}
#     #print(states)
#     states = {
#         'm_Ih' : 'Ih activation gate', #hyperpolarization-activated currents
#         'm_Im' : 'Im activation gate', #muscarinic agonist-activated currents

#         'm_KPst': 'iK_Pst activation gate',
#         'h_KPst': 'iK_Pst inactivation gate',   
#         'm_KTst': 'iK_Tst activation gate',
#         'h_KTst': 'iK_Tst inactivation gate', 

#         'm_NaTat': 'iNaTa_t activation gate',
#         'h_NaTat': 'iNaTa_t inactivation gate', 
#         'm_NaTs2t': 'iNaTs2_t activation gate',
#         'h_NaTs2t': 'iNaTs2_t inactivation gate',
#         'm_NapEt2': 'iNap_Et2 activation gate',
#         'h_NapEt2': 'iNap_Et2 inactivation gate'
#     }

#     # ------------------------------ Gating states kinetics ------------------------------

#     # "Ih and Im"    
#     # @classmethod
#     # def alpham_Ih(cls, Vm):
#     #     return 0.001 * 6.43 * cls.vtrap((Vm + cls.VIh), 11.9) * 1e3  # s-1 -> ms-1??

#     # @classmethod
#     # def betam_Ih(cls, Vm):
#     #     return 0.001 * 193 * np.exp(Vm/33.1) * 1e3  # s-1
    

    
#     # @classmethod
#     # def alpham_Im(cls, Vm):
#     #     return 3.3e-3 * np.exp(2.5*0.04*(Vm - -35)) * 1e3  # s-1

#     # @classmethod
#     # def betam_Im(cls, Vm):
#     #     return 3.3e-3 * np.exp(-2.5*0.04*(Vm - -35)) * 1e3  # s-1
    


#     # "K_Pst and K_Tst"    
#     # @classmethod
#     # def taum_KPst(cls, Vm):
#     #     v = Vm + 10
#     #     qt = 2.3**((cls.T_C-21)/10)
#     #     return (1.25+175.03*np.exp(-v * -0.026))/qt * 1e3 if v<-50 else ((1.25+13*np.exp(-v*0.026)))/qt * 1e3  # s-1

#     # @classmethod
#     # def m_KPstinf(cls, Vm):
#     #     v = Vm + 10
#     #     return (1/(1 + np.exp(-(v+1)/12))) * 1e3  # s-1

#     # @classmethod
#     # def tauh_KPst(cls, Vm):
#     #     v = Vm + 10
#     #     qt = 2.3**((cls.T_C-21)/10)
#     #     return (360+(1010+24*(v+55))*np.exp(-((v+75)/48)**2))/qt * 1e3  # s-1

#     # @classmethod
#     # def h_KPstinf(cls, Vm):
#     #     v = Vm + 10
#     #     return 1/(1 + np.exp(-(v+54)/-11)) * 1e3  # s-1
    


#     # @classmethod
#     # def taum_KTst(cls, Vm):
#     #     v = Vm + 10
#     #     qt = 2.3**((cls.T_C-21)/10)
#     #     return (0.34+0.92*np.exp(-((v+71)/59)**2))/qt * 1e3  # s-1

#     # @classmethod
#     # def m_KTstinf(cls, Vm):
#     #     v = Vm + 10
#     #     return 1/(1 + np.exp(-(v+0)/19)) * 1e3  # s-1

#     # @classmethod
#     # def tauh_KTst(cls, Vm):
#     #     v = Vm + 10
#     #     qt = 2.3**((cls.T_C-21)/10)
#     #     return (8+49*np.exp(-((v+73)/23)**2))/qt * 1e3  # s-1

#     # @classmethod
#     # def h_KTstinf(cls, Vm):
#     #     v = Vm + 10
#     #     return 1/(1 + np.exp(-(v+66)/-10)) * 1e3  # s-1
    

    
#     # "NaTa_t, NaTs2_t and Nap_Et2"    
#     # @classmethod
#     # def alpham_NaTat(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-38) else Vm
#     #     return (0.182 * (v- -38))/(1-(np.exp(-(v- -38)/6))) * 1e3  # s-1

#     # @classmethod
#     # def betam_NaTat(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-38) else Vm
#     #     return (0.124 * (-v -38))/(1-(np.exp(-(-v -38)/6))) * 1e3  # s-1

#     # @classmethod
#     # def alphah_NaTat(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-66) else Vm
#     #     return (-0.015 * (v- -66))/(1-(np.exp((v- -66)/6))) * 1e3  # s-1

#     # @classmethod
#     # def betah_NaTat(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-66) else Vm
#     #     return (-0.015 * (-v -66))/(1-(np.exp((-v -66)/6))) * 1e3  # s-1



#     # @classmethod
#     # def alpham_NaTs2t(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-32) else Vm
#     #     return (0.182 * (v- -32))/(1-(np.exp(-(v- -32)/6))) * 1e3  # s-1

#     # @classmethod
#     # def betam_NaTs2t(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-32) else Vm
#     #     return (0.124 * (-v -32))/(1-(np.exp(-(-v -32)/6))) * 1e3  # s-1

#     # @classmethod
#     # def alphah_NaTs2t(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-60) else Vm
#     #     return (-0.015 * (v- -60))/(1-(np.exp((v- -60)/6))) * 1e3  # s-1

#     # @classmethod
#     # def betah_NaTs2t(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-60) else Vm
#     #     return (-0.015 * (-v -60))/(1-(np.exp((-v -60)/6))) * 1e3  # s-1
    


#     # @classmethod
#     # def alpham_NapEt2(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-38) else Vm
#     #     return (0.182 * (v- -38))/(1-(np.exp(-(v- -38)/6))) * 1e3  # s-1

#     # @classmethod
#     # def betam_NapEt2(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-38) else Vm
#     #     return (0.124 * (-v -38))/(1-(np.exp(-(-v -38)/6))) * 1e3  # s-1

#     # @classmethod
#     # def alphah_NapEt2(cls, Vm):
#     #     v = Vm + 0.0001 if (Vm ==-17 or Vm == -64.4) else Vm
#     #     return -2.88e-6 * (v + 17) / (1 - np.exp((v + 17)/4.63)) * 1e3  # s-1

#     @classmethod
#     def check(cls, Vm):
#         return 0
#     exec(
#     """@classmethod
# def functie(cls, Vm):
#     return 0"""
#     )
    
#     exec('test = functie')
    

#     '''how do I add the @classmethod? -> fixed'''
#     mod_number, variables = 0, 0
#     for e,f in zip(mod_files,mod_names):
#         if mod_number in hits:
#             name_mod = f.replace('.mod','')
#             for x in ['m','h']:
#                 gating_param =  x+'_'+name_mod
#                 if gating_param in states:
#                     function_template_alpha = f"""@classmethod
# def alpha{gating_param}_fct(cls,Vm):
#     variables = tf.gating_from_PROCEDURES(Vm=Vm,list_mod={[str(i) for i in e]},mod_name='{f}')
#     return variables['{x}'+'alpha']"""
#                     function_template_beta = f"""@classmethod
# def beta{gating_param}_fct(cls,Vm):
#     variables = tf.gating_from_PROCEDURES(Vm=Vm,list_mod={[str(i) for i in e]},mod_name='{f}')
#     return variables['{x}'+'beta']"""
#                     exec(function_template_alpha)#,{**globals(),**{'e':e, 'f':f}})
#                     exec(f'alpha{gating_param} = alpha{gating_param}_fct')
#                     gating_kinetics[f'alpha{gating_param}'] = eval(f'alpha{gating_param}')
#                     #print(f'alpha{gating_param} = alpha{gating_param}_fct')
#                     exec(function_template_beta)#,{**globals(),**{'e':e, 'f':f}})
#                     exec(f'beta{gating_param} = beta{gating_param}_fct')
#                     gating_kinetics[f'beta{gating_param}'] = eval(f'beta{gating_param}')
#         #print(variables,"\n")
#         mod_number += 1
#         #print(f"alpham_Ca_HVA: {alpham_Ca_HVA}, check: {check}")

#     # ------------------------------ States derivatives ------------------------------

#     @classmethod
#     def derStates(cls):
#         d_derStates = {}
#         for e in cls.states:
#             #print(e)
#             #print(eval(f"cls.alpha{e}"))
#             #print(f"lambda Vm, x: cls.alpha{e}(Vm) * (1 - x['{e}']) - cls.beta{e}(Vm) * x['{e}']")
#             #exec(f"d_derStates[e] = lambda Vm, x: cls.alpha{e}(Vm) * (1 - x['{e}']) - cls.beta{e}(Vm) * x['{e}']")
#             d_derStates[e] = lambda Vm, x: eval(f"cls.alpha{e}")(Vm) * (1 - x[e]) - eval(f"cls.beta{e}")(Vm) * x[e]
#         #print('derstates 1: ',d_derStates)
        
#         return d_derStates

#         d_derStates = {
#         #     'm_Ca_HVA' : lambda Vm, x: cls.alpham_Ca_HVA(Vm) * (1 - x['m_Ca_HVA']) - cls.betam_Ca_HVA(Vm) * x['m_Ca_HVA'],
#         #     'h_Ca_HVA' : lambda Vm, x: cls.alphah_Ca_HVA(Vm) * (1 - x['h_Ca_HVA']) - cls.betah_Ca_HVA(Vm) * x['h_Ca_HVA'],
#         #     'm_Ca_LVAst' : lambda Vm, x: cls.alpham_Ca_LVAst(Vm) * (1 - x['m_Ca_LVAst']) - cls.betam_Ca_LVAst(Vm) * x['m_Ca_LVAst'],
#         #     'h_Ca_LVAst' : lambda Vm, x: cls.alphah_Ca_LVAst(Vm) * (1 - x['h_Ca_LVAst']) - cls.betah_Ca_LVAst(Vm) * x['h_Ca_LVAst'],

#             'm_Ih': lambda Vm, x:  eval(f"cls.alpham_Ih")(Vm) * (1 - x['m_Ih']) - eval(f"cls.betam_Ih")(Vm) * x['m_Ih'],
#             # 'm_Im': lambda Vm, x: cls.alpham_Im(Vm) * (1 - x['m_Im']) - cls.betam_Im(Vm) * x['m_Im'],

#         #     'm_NaTat': lambda Vm, x: cls.alpham_NaTat(Vm) * (1 - x['m_NaTat']) - cls.betam_NaTat(Vm) * x['m_NaTat'],
#         #     'h_NaTat': lambda Vm, x: cls.alphah_NaTat(Vm) * (1 - x['h_NaTat']) - cls.betah_NaTat(Vm) * x['h_NaTat'],
#         #     'm_NaTs2t': lambda Vm, x: cls.alpham_NaTs2t(Vm) * (1 - x['m_NaTs2t']) - cls.betam_NaTs2t(Vm) * x['m_NaTs2t'],
#         #     'h_NaTs2t': lambda Vm, x: cls.alphah_NaTs2t(Vm) * (1 - x['h_NaTs2t']) - cls.betah_NaTs2t(Vm) * x['h_NaTs2t'],
#         #     'm_NapEt2': lambda Vm, x: cls.alpham_Nap_Et2(Vm) * (1 - x['m_NapEt2']) - cls.betam_Nap_Et2(Vm) * x['m_NapEt2'],
#         #     'h_NapEt2': lambda Vm, x: cls.alphah_NapEt2(Vm) * (1 - x['h_NapEt2']) - cls.betah_NapEt2(Vm) * x['h_NapEt2'],
#         #     'm_SKv3_1': lambda Vm, x: cls.alphah_m_SKv3_1(Vm) * (1 - x['h_m_SKv3_1']) - cls.betah_m_SKv3_1(Vm) * x['h_m_SKv3_1']

#         }
#         # print('derstates : ',d_derStates)
    
#         return d_derStates

#     # ------------------------------ Steady states ------------------------------

#     @classmethod
#     def steadyStates(cls):
#         d_steadyStates = {}
#         #print(getattr(cls,'alpham_Ih'),getattr(cls,'alpham_Ih'),alpham_Ih)
#         for e in cls.states:
#             #print(e)
#             #print(eval(f"cls.alpha{e}"))
#             #print(f"lambda Vm, x: cls.alpha{e}(Vm) * (1 - x['{e}']) - cls.beta{e}(Vm) * x['{e}']")
#             #exec(f"d_derStates[e] = lambda Vm, x: cls.alpha{e}(Vm) * (1 - x['{e}']) - cls.beta{e}(Vm) * x['{e}']")
#             d_steadyStates[e] = lambda Vm : eval(f"cls.alpha{e}")(Vm) / (eval(f"cls.alpha{e}")(Vm) + eval(f"cls.beta{e}")(Vm))
#         #print('steadystates: ',d_steadyStates)
#         #return d_steadyStates
#         #print(globals().keys(),'\n\n\n',locals().keys())
#         return {
#             'm_Ih': lambda Vm: cls.dummy * gating_kinetics[f'alpham_Ih'](Vm) / (gating_kinetics[f'alpham_Ih'](Vm) + gating_kinetics[f'alpham_Ih'](Vm)),
#         #     'm_Im': lambda Vm: cls.alpham_Im(Vm) / (cls.alpham_Im(Vm) + cls.betam_Im(Vm)),

#         #     'm_KPst': lambda Vm: cls.m_KPstinf(Vm),
#         #     'h_KPst': lambda Vm: cls.h_KPstinf(Vm),
#         #     'm_KTst': lambda Vm: cls.m_KTstinf(Vm),
#         #     'h_KTst': lambda Vm: cls.h_KTstinf(Vm),  

#         #     'm_NaTat': lambda Vm: cls.alpham_NaTat(Vm) / (cls.alpham_NaTat(Vm) + cls.betam_NaTat(Vm)),
#         #     'h_NaTat': lambda Vm: cls.alphah_NaTat(Vm) / (cls.alphah_NaTat(Vm) + cls.betah_NaTat(Vm)),
#         #     'm_NaTs2t': lambda Vm: cls.alpham_NaTs2t(Vm) / (cls.alpham_NaTs2t(Vm) + cls.betam_NaTs2t(Vm)),
#         #     'h_NaTs2t': lambda Vm: cls.alphah_NaTs2t(Vm) / (cls.alphah_NaTs2t(Vm) + cls.betah_NaTs2t(Vm)),
#         #     'm_NapEt2': lambda Vm: cls.alpham_NapEt2(Vm) / (cls.alpham_NapEt2(Vm) + cls.betam_NapEt2(Vm)),
#         #     'h_NapEt2': lambda Vm: cls.alphah_NapEt2(Vm) / (cls.alphah_NapEt2(Vm) + cls.betah_NapEt2(Vm)),
#         }

#     # ------------------------------ Membrane currents ------------------------------
    
#     # @classmethod
#     # def ihcn(cls, m_Ih, Vm):
#     #     ''' Ih '''
#     #     gIh = cls.gIhbar*m_Ih
#     #     return gIh*(Vm - cls.ehcn)  # mA/m2 #ehcn != ENa: hyperpolarization-activated
#     #                                 # cyclic nucleotide-gated channels
    
#     # @classmethod
#     # def ik(cls, m_Im, Vm):
#     #     ''' Im '''
#     #     gIm = cls.gImbar*m_Im
#     #     return gIm*(Vm - cls.EK)  # mA/m2

#     # @classmethod
#     # def iKPst(cls, m_KPst, h_KPst, Vm):
#     #     ''' I_K_Pst '''
#     #     gKPst = cls.gKPstbar*m_KPst**2*h_KPst
#     #     return gKPst*(Vm - cls.EK)  # mA/m2
    
#     # @classmethod
#     # def iKTst(cls, m_KTst, h_KTst, Vm):
#     #     ''' I_K_Tst '''
#     #     gKTst = cls.gKTstbar*m_KTst**4*h_KTst
#     #     return gKTst*(Vm - cls.EK)  # mA/m2
    
#     # @classmethod
#     # def iNaTat(cls, m_NaTat, h_NaTat, Vm):
#     #     ''' I_NaTa_t '''
#     #     gNaTat = cls.gNaTatbar*m_NaTat**3*h_NaTat
#     #     return gNaTat*(Vm - cls.ENa)  # mA/m2
    
#     # @classmethod
#     # def iNaTs2t(cls, m_NaTs2t, h_NaTs2t, Vm):
#     #     ''' I_NaTs2_t '''
#     #     gNaTs2t = cls.gNaTs2tbar*m_NaTs2t**3*h_NaTs2t
#     #     return gNaTs2t*(Vm - cls.ENa)  # mA/m2
    
#     # @classmethod
#     # def iNapEt2(cls, m_NapEt2, h_NapEt2, Vm):
#     #     ''' I_Nap_Et2 '''
#     #     gNapEt2 = cls.gNapEt2bar*m_NapEt2**3*h_NapEt2
#     #     return gNapEt2*(Vm - cls.ENa)  # mA/m2

#     @classmethod
#     def currents(cls):
#         return {
#             # 'ihcn': lambda Vm, x: cls.ihcn(x['m_Ih'], Vm),
#             # 'ik': lambda Vm, x: cls.ik(x['m_Im'], Vm),

#             # 'iKPst': lambda Vm, x: cls.iKPst(x['m_KPst'], x['h_KPst'], Vm),
#             # 'iKTst': lambda Vm, x: cls.iKTst(x['m_KTst'], x['h_KTst'], Vm),

#             # 'iNaTat': lambda Vm, x: cls.iNaTat(x['m_NaTat'], x['h_NaTat'], Vm),
#             # 'iNaTs2t': lambda Vm, x: cls.iNaTs2t(x['m_NaTs2t'], x['h_NaTs2t'], Vm),
#             # 'iNapEt2': lambda Vm, x: cls.iNapEt2(x['m_NapEt2'], x['h_NapEt2'], Vm),
#         }
    

# # klasse = RealisticNeuron()
# # print(f"alpham_Ca_HVA: {klasse.alpham_Ca_HVA}, check: {klasse.check}")
# # print(f"alpham_Ca_HVA: {klasse.alpham_Ca_HVA}, check: {klasse.check}")
# # print(f"functie: {klasse.test}")