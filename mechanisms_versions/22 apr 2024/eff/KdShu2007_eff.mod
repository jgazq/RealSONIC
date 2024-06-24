TITLE K-D
: K-D current for prefrontal cortical neuron ------Yuguo Yu  2007

NEURON {
     THREADSAFE
	SUFFIX KdShu2007
	USEION k WRITE ik
	RANGE  gkbar, ik, ek
	GLOBAL minf, mtau, hinf, htau
	RANGE Adrive, Vm, y, Fdrive, A_t : section specific
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
} 

ASSIGNED {
	ik 		(mA/cm2)
	minf 		mtau (ms)	 	
	hinf 		htau (ms)	 	
	A_t  (kPa)
	y
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_KdShu2007(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_KdShu2007(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_KdShu2007(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betah_KdShu2007(A(kPa), Q(nC/cm2)) (/ms)
 

STATE {
	 m h
}
BREAKPOINT {
printf("KdShu2007.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
        SOLVE states METHOD cnexp
       ik = gkbar * m*h*(Vm-ek)
} 

INITIAL {
printf("KdShu2007.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
printf("KdShu2007.mod: \n")
printf("V = %g, alpha = %g, beta = %g\n",V(A_t,y), alpham_KdShu2007(A_t, y), betam_KdShu2007(A_t, y))
	m= alpham_KdShu2007(A_t, y) / (alpham_KdShu2007(A_t, y) + betam_KdShu2007(A_t, y))
printf("KdShu2007.mod: \n")
printf("V = %g, alpha = %g, beta = %g\n",V(A_t,y), alphah_KdShu2007(A_t, y), betah_KdShu2007(A_t, y))
	h= alphah_KdShu2007(A_t, y) / (alphah_KdShu2007(A_t, y) + betah_KdShu2007(A_t, y))
}

DERIVATIVE states {   
        m' = alpham_KdShu2007(A_t, y) * (1 - m) - betam_KdShu2007(A_t, y) * m
        h' = alphah_KdShu2007(A_t, y) * (1 - h) - betah_KdShu2007(A_t, y) * h
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}