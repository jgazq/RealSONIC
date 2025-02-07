TITLE K-D
: K-D current for prefrontal cortical neuron ------Yuguo Yu  2007

NEURON {
     THREADSAFE
	SUFFIX KdShu20072
	USEION k WRITE ik
	RANGE  gkbar, ik, ek
	GLOBAL minf, mtau, hinf, htau
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
}

PARAMETER {
	stimon       : Stimulation state
	Fdrive (kHz) : Stimulation frequency
	Adrive (kPa) : Stimulation amplitude
	detailed     : Simulation type
	gkbar = 0.1   	(mho/cm2)	
								
	celsius
	ek = -100	(mV)            : must be explicitly def. in hoc
	v (nC/cm2)
	Vm (mV)
	vhalfm=-43  (mV)
	km=8
	vhalfh=-67  (mV) 
      kh=7.3
	q10=2.3
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	PI = (pi) (1) :in order to use the constant pi
} 

ASSIGNED {
	ik 		(mA/cm2)
	minf 		mtau (ms)	 	
	hinf 		htau (ms)	 	
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)
FUNCTION_TABLE alpham_KdShu20072(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betam_KdShu20072(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_KdShu20072(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betah_KdShu20072(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE A_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)
FUNCTION_TABLE B_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)
 

STATE {
	 m h
}
BREAKPOINT {
	update()
        SOLVE states METHOD cnexp
       ik = gkbar * m*h*(Vm-ek)
} 

INITIAL {
	update()
	m= alpham_KdShu20072(A_t, y, a1, b1) / (alpham_KdShu20072(A_t, y, a1, b1) + betam_KdShu20072(A_t, y, a1, b1))
	h= alphah_KdShu20072(A_t, y, a1, b1) / (alphah_KdShu20072(A_t, y, a1, b1) + betah_KdShu20072(A_t, y, a1, b1))
}

DERIVATIVE states {   
        m' = alpham_KdShu20072(A_t, y, a1, b1) * (1 - m) - betam_KdShu20072(A_t, y, a1, b1) * m
        h' = alphah_KdShu20072(A_t, y, a1, b1) * (1 - h) - betah_KdShu20072(A_t, y, a1, b1) * h
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}