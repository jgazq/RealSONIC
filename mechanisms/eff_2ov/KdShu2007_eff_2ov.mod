TITLE K-D
: K-D current for prefrontal cortical neuron ------Yuguo Yu  2007

NEURON {
     THREADSAFE
	SUFFIX KdShu2007
	USEION k WRITE ik
	RANGE  gkbar, ik, ek
	GLOBAL minf, mtau, hinf, htau
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1, a2, b2 : section (even segment) specific
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
	a2  (nC/cm2)
	b2  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (mV)
FUNCTION_TABLE alpham_KdShu2007(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE betam_KdShu2007(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_KdShu2007(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE betah_KdShu2007(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
 

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
	m= alpham_KdShu2007(A_t, y, a1, b1, a2, b2) / (alpham_KdShu2007(A_t, y, a1, b1, a2, b2) + betam_KdShu2007(A_t, y, a1, b1, a2, b2))
	h= alphah_KdShu2007(A_t, y, a1, b1, a2, b2) / (alphah_KdShu2007(A_t, y, a1, b1, a2, b2) + betah_KdShu2007(A_t, y, a1, b1, a2, b2))
}

DERIVATIVE states {   
        m' = alpham_KdShu2007(A_t, y, a1, b1, a2, b2) * (1 - m) - betam_KdShu2007(A_t, y, a1, b1, a2, b2) * m
        h' = alphah_KdShu2007(A_t, y, a1, b1, a2, b2) * (1 - h) - betah_KdShu2007(A_t, y, a1, b1, a2, b2) * h
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}