""""this file just contains garbage from other files that could be useful in the future (probably not but you never know)"""

"how to skip a value/iteration in a for/while loop"

                if_statement = eval(re.search(if_pattern,e).group().replace('exp','np.exp').replace('v','Vm').replace('^','**').lower(),{**globals(),**variables_dict})
                while(not re.search("}",list_mod[i+1])):
                    calc_eq(e,var_pattern,math_pattern,equation_pattern,variables_dict) if if_statement else None
                    [next(i) for _ in range(1)]
                    [next(e) for _ in range(1)]
                if re.search('else',e):
                    while(not re.search("}",list_mod[i+1])):
                        calc_eq(e,var_pattern,math_pattern,equation_pattern,variables_dict) if not if_statement else None
                        [next(i) for _ in range(1)]
                        [next(e) for _ in range(1)]    

----------------------------------------------------------------------------------------   
"first half attempt to dynamically create functions" 
    # @classmethod
    # def betah_NapEt2(cls, Vm):
    #     v = Vm + 0.0001 if (Vm ==-17 or Vm == -64.4) else Vm
    #     return 6.94e-6 * (v + 64.4) / (1 - np.exp(-(v + 64.4)/2.63)) * 1e3  # s-1
    

    # for e in states:
    #     function_template = f""" @classmethod; def beta{e}: TODO"""
    #     def create_funcs():
    #         @classmethod
    #         def betah_testerdetest(cls, Vm):
    #             v = Vm + 0.0001 if (Vm ==-17 or Vm == -64.4) else Vm
    #             return 6.94e-6 * (v + 64.4) / (1 - np.exp(-(v + 64.4)/2.63)) * 1e3  # s-1    
    #         return betah_testerdetest
    # betah_testerdetest = create_func()

----------------------------------------------------------------------------------------
"how to add dynamic method to a class (with the help of Jorn)"
    def add_dynamic_method(self, method_string,dynamic_method_name):

        method_code = compile(method_string, '<string>', 'exec')
        dynamic_method = types.FunctionType(method_code.co_consts[0], globals(), "dynamic_class_method")
        setattr(RealisticNeuron, dynamic_method_name, dynamic_method)

# alpha_method = types.MethodType(alpham_NapEt2,None,RealisticNeuron)
# RealisticNeuron.alpham_NapEt2 = alpha_method

# gating_param = 'm_Ca_HVA'
# function_template_alpha = f"""@classmethod
# def alpha{gating_param}(cls,Vm):
#     variables = tf.gating_from_PROCEDURES(Vm=Vm,list_mod=['test'],mod_name='Ca_HVA')
#     return variables['m'+'alpha']"""
# klasse.add_dynamic_method(function_template_alpha,f"alpha{gating_param}")
-------------------------------------------------------------------------------------------------
"trick to create different lambdas by using the same variables -> if variables are changed inside lambda statement then lambda is also changed"
# dic_funcs = {}

# dates = ['2018', '2019', '2020']

# idx = ['0', '1', '2']

# def printt(text):
#     print(text)


# for d in dates:
#     for i in idx:

#         #print(d, i)

#         dic_funcs[d + '_' + i] = lambda a=i, b=d:printt(b + '_' + a) #dic_funcs[d + '_' + i] = lambda:printt(d + '_' + i)


# for key, item in dic_funcs.items():
#     print(key)
#     item()      
#     print('-----')  


------------------------------------------------------------------------------
"real_neuron.py backup"
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


------------------------------------------------------------------------------------------------------------
def mod_to_eff(mech_folder, restrictions = None):
    """read all .mod files in a directory and create the effective ones in parallel"""

    mod_files = []
    mod_names = []
    for root, dirs, files in os.walk(mech_folder):
        for file in files:
            if restrictions:
                if file.endswith(".mod") and not file.endswith("_eff.mod") and file.replace(".mod","") in restrictions:
                    file_dupl = file.replace(".mod","_eff.mod")
                    with open (os.path.join(root,file)) as f, open(os.path.join(root,file_dupl),'w') as dupl:
                        for line in f:
                            if '{' in line:
                                dupl.write(line)
            elif not file.endswith("_eff.mod"):
                if file.endswith(".mod"):
                    file_dupl = file.replace(".mod","_eff.mod")
                    print(os.path.join(root,file),os.path.join(root,"eff",file_dupl))
                    with open (os.path.join(root,file)) as f, open(os.path.join(root,"eff",file_dupl),'w') as dupl:
                        for line in f:
                            if '{' in line:
                                dupl.write(line)
    return mod_files, mod_names 