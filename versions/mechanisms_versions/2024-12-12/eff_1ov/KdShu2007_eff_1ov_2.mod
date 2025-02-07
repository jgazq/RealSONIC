TITLE K-D
: K-D current for prefrontal cortical neuron ------Yuguo Yu  2007

NEURON {
     THREADSAFE
	SUFFIX KdShu20072
	USEION k WRITE ik
	RANGE  gkbar, ik, ek
	GLOBAL minf, mtau, hinf, htau
	RANGE Adrive, Vm, y, Fdrive, A_t, q1, f1 : section (even segment) specific
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
	q1  (nC/cm2)
	f1  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE A_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE phi_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (rad)
FUNCTION_TABLE alpham_KdShu20072(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betam_KdShu20072(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE alphah_KdShu20072(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betah_KdShu20072(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
 

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
	m= alpham_KdShu20072(A_t, y, q1, f1) / (alpham_KdShu20072(A_t, y, q1, f1) + betam_KdShu20072(A_t, y, q1, f1))
	h= alphah_KdShu20072(A_t, y, q1, f1) / (alphah_KdShu20072(A_t, y, q1, f1) + betah_KdShu20072(A_t, y, q1, f1))
}

DERIVATIVE states {   
        m' = alpham_KdShu20072(A_t, y, q1, f1) * (1 - m) - betam_KdShu20072(A_t, y, q1, f1) * m
        h' = alphah_KdShu20072(A_t, y, q1, f1) * (1 - h) - betah_KdShu20072(A_t, y, q1, f1) * h
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}