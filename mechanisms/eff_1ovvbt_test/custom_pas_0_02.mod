TITLE Custom passive current

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}
NEURON {
    SUFFIX pas_eff0_02
    NONSPECIFIC_CURRENT i : passive leakage current
    RANGE g, e
    RANGE Adrive, Vm, y, Fdrive, A_t : section specific
    RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
    RANGE q1, f1
    POINTER V 
    RANGE V_val : contains the specific value of LUT V for a specific segment 
}

PARAMETER {
    stimon       : Stimulation state
    Fdrive (kHz) : Stimulation frequency
    Adrive (kPa) : Stimulation amplitude
    detailed     : Simulation type
    g      (S/cm2)
    e      (mV)
}

ASSIGNED {
    v   (nC/cm2)
    Vm  (mV)
    i   (mA/cm2)
    A_t  (kPa)
    y
    q1  (nC/cm2)
    f1  (rad)
    V_val (mV) 
}



BREAKPOINT {
Vm = v/0.02 :Cm0 = 0.02 uF/cm2
y = v
    i = g * (Vm - e)
}