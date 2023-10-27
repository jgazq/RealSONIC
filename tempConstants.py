import numpy as np
"-----CONSTANTS-----"
RHO = 1e3 #medium density (kg/m3)
C = 1490.0 #medium speed of sound (m)

"-----GATING PARAMETERS-----"
#RESTING PARAMETERS
Cm0 = 1e-2   # Membrane capacitance (F/m2)
Vm0 = -71.9  # Membrane potential (mV)

#REVERSAL PARAMETERS (Mv)
ENa = 50.0     # Sodium
EK = -85.0 #-90.0     # Potassium
ELeak = -70.3  # Non-specific leakage
#ABERRA/BBP
ehcn =  -45.0 #Ih (mV), hyperpolarization-activated 
              #cyclic nucleotide-gated channels

#MAXIMAL CHANNEL CONDUCTANCES (S/m2)
gNabar = 560.0  # Sodium
gKdbar = 60.0   # Delayed-rectifier Potassium
gLeak = 0.205   # Non-specific leakage
#ABERRA/BBP
#BASAL
# gIhbar = 8e-05*1e4 #0.00001*1e4          #Ih (S/cm2 -> S/m2)
# #APICAL
# gNaTs2_tbar = 0.012009*1e4 #0.00001*1e4   #Nap_Et2 (S/cm2 -> S/m2)
# gSKv3_1bar = 0.000513*1e4 #0.00001*1e4   #Nap_Et2 (S/cm2 -> S/m2)
# dist_2_soma = 0.001
# gIhbar_Ih = (-0.869600 + 2.087000*np.exp((dist_2_soma-0.000000)*0.003100))*0.000080
# gImbar = 0.00074*1e4 #0.00001*1e4        #Im (S/cm2 -> S/m2)
# gNaTa_tbar = 3.429725*1e4 #0.00001*1e4    #NaTa_t (S/cm2 -> S/m2)
# gK_Tstbar = 0.001035*1e4 #0.00001*1e4     #K_Tst (S/cm2 -> S/m2)
# gNap_Et2bar = 0.009803*1e4 #0.00001*1e4   #Nap_Et2 (S/cm2 -> S/m2)
# gSK_E2bar = 0.008085*1e4 #0.00001*1e4   #Nap_Et2 (S/cm2 -> S/m2)

# gK_Pstbar = 0.959296*1e4 #0.00001*1e4     #K_Pst (S/cm2 -> S/m2)
# gNaTs2_tbar = 0.926705*1e4 #0.00001*1e4   #NaTs2_t (S/cm2 -> S/m2)

#ADDITIONAL
VT = -56.2  # Spike threshold adjustment parameter (mV)
VIh = 154.9 # Mv
T_C = 34 # Temperature (°C)



class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'