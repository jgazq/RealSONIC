import numpy as np


"-----PHYSICAL CONSTANTS-----"
RHO = 1e3 #medium density (kg/m3)
C = 1490.0 #medium speed of sound (m)

"-----GATING PARAMETERS-----"
"""RESTING PARAMETERS"""
Cm0 = 1e-2   # Membrane capacitance (F/m2)
Vm0 = -71.9  # Membrane potential (mV)

"""REVERSAL PARAMETERS (Mv)"""
ENa = 50.0     # Sodium
EK = -85.0 #-90.0     # Potassium
ELeak = -70.3  # Non-specific leakage
"""ABERRA/BBP"""
ehcn =  -45.0 #Ih (mV), hyperpolarization-activated 
              #cyclic nucleotide-gated channels

"""MAXIMAL CHANNEL CONDUCTANCES (S/m2)"""
gNabar = 560.0  # Sodium
gKdbar = 60.0   # Delayed-rectifier Potassium
gLeak = 0.205   # Non-specific leakage
"""ABERRA/BBP"""
"""BASAL"""
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

"""ADDITIONAL"""
VT = -56.2  # Spike threshold adjustment parameter (mV)
VIh = 154.9 # Mv
T_C = 34 # Temperature (°C)

"-----TEXT FORMATTING-----"
"""COLORING LOGS"""
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

"""-----REGULAR EXPRESSION PATTERNS-----"""
"""WRITE_REALNEURON"""
alpha_x_pattern = ".*[Aa]lpha.*=.*v.*|[Aa]lpha.*.*=.*v.*" # '.*' can be anything, '|': OR
beta_x_pattern = ".*[Bb]eta.*=.*v.*|[Bb]eta.*.*=.*v.*" # '.*' can be anything, '|': OR
tau_x_pattern = ".*[Tt]au.*=.*v.*|[Tt]au.*.*=.*v.*" # '.*' can be anything, '|': OR
x_inf_pattern = ".*[Ii]nf.*=.*v.*|[Ii]nf.*.*=.*v.*" # '.*' can be anything, '|': OR

alpha_pattern = "[a-zA-Z]_[Aa]lpha|[Aa]lpha_[a-zA-Z]|[a-zA-Z][Aa]lpha|[Aa]lpha[a-zA-Z]"
beta_pattern = "[a-zA-Z]_[Bb]eta|[Bb]eta_[a-zA-Z]|[a-zA-Z][Bb]eta|[Bb]eta[a-zA-Z]"
tau_pattern = "[a-zA-Z]_[Tt]au|[Tt]au_[a-zA-Z]|[a-zA-Z][Tt]au|[Tt]au[a-zA-Z]"
inf_pattern = "[a-zA-Z]_[Ii]nf|[Ii]nf_[a-zA-Z]|[a-zA-Z][Ii]nf|[Ii]nf[a-zA-Z]"

var_pattern = "[a-zA-Z_][a-zA-Z0-9_]*"
math_pattern = "[0-9\.\+\-\*/\(\)a-zA-Z]*[ 0-9\.\+\-\*/\(\)a-zA-Z^_]*[0-9\.\+\-\*/\(\)a-zA-Z]*" #removal of \t, \n and spaces around the formula
equation_pattern = r"[0-9.\+\-\*/\(\)a-zA-Z^=_ ]+" #"[0-9\.\+\-\*/\(\)a-zA-Z]*[ 0-9\.\+\-\*/\(\)a-zA-Z^=_]*[0-9\.\+\-\*/\(\)a-zA-Z]*" #removal of \t, \n and spaces around the formula
if_pattern = "\(.*\)"

""""WRITE_NMODL"""
neuronvoltage_pattern = "v.*\(m[Vv]\)"
vRHS_pattern = "\=.*v"
vsecluded_pattern = "\Wv\W" #v pattern that is not part of a word
block_init_pattern = "^[A-Z][A-Z]*.*\{"
block_pattern = "^[A-Z][A-Z]*"
state_pattern = "[a-zA-Z]"
stateder_pattern = "[a-zA-Z]\'"
onelineBLOCK_pattern = "^[a-zA-Z]+.*\{.*\}"