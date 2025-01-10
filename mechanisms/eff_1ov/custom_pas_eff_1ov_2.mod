TITLE Custom passive current

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}
NEURON {
    SUFFIX pas_eff2
    NONSPECIFIC_CURRENT i : passive leakage current
    RANGE g, e
    RANGE Adrive, Vm, y, Fdrive, A_t : section specific
    RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
    RANGE q1, f1
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
}

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE A_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE phi_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (rad)

INCLUDE "update.inc"

BREAKPOINT {
    update()
    i = g * (Vm - e)
}