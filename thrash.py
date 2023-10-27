""""this file just contains garbage from other files that could be useful in the future (probably not but you never know)"""

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