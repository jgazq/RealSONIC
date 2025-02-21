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

	POINTER V_table, alpham_KdShu20072_table, betam_KdShu20072_table, alphah_KdShu20072_table, betah_KdShu20072_table, A_1_table, B_1_table
	POINTER A_arr, Q_arr, A1_arr, B1_arr
	RANGE A_s, Q_s, A1_s, B1_s
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

	V_table  alpham_KdShu20072_table  betam_KdShu20072_table  alphah_KdShu20072_table  betah_KdShu20072_table  A_1_table  B_1_table  
	A_arr  Q_arr  A1_arr    B1_arr
	A_s  Q_s  A1_s  B1_s
}

INCLUDE "update.inc"
INCLUDE "interp.inc"

FUNCTION fV() { 
VERBATIM
	double V_value;
	V_value = interp4D(_p_V_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(V_value);
ENDVERBATIM
	fV = V_value
}

FUNCTION falpham_KdShu20072() { 
VERBATIM
	double alpham_KdShu20072_value;
	alpham_KdShu20072_value = interp4D(_p_alpham_KdShu20072_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_KdShu20072_value);
ENDVERBATIM
	falpham_KdShu20072 = alpham_KdShu20072_value
}

FUNCTION fbetam_KdShu20072() { 
VERBATIM
	double betam_KdShu20072_value;
	betam_KdShu20072_value = interp4D(_p_betam_KdShu20072_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_KdShu20072_value);
ENDVERBATIM
	fbetam_KdShu20072 = betam_KdShu20072_value
}

FUNCTION falphah_KdShu20072() { 
VERBATIM
	double alphah_KdShu20072_value;
	alphah_KdShu20072_value = interp4D(_p_alphah_KdShu20072_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_KdShu20072_value);
ENDVERBATIM
	falphah_KdShu20072 = alphah_KdShu20072_value
}

FUNCTION fbetah_KdShu20072() { 
VERBATIM
	double betah_KdShu20072_value;
	betah_KdShu20072_value = interp4D(_p_betah_KdShu20072_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_KdShu20072_value);
ENDVERBATIM
	fbetah_KdShu20072 = betah_KdShu20072_value
}

FUNCTION fA_1() { 
VERBATIM
	double A_1_value;
	A_1_value = interp4D(_p_A_1_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(A_1_value);
ENDVERBATIM
	fA_1 = A_1_value
}

FUNCTION fB_1() { 
VERBATIM
	double B_1_value;
	B_1_value = interp4D(_p_B_1_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(B_1_value);
ENDVERBATIM
	fB_1 = B_1_value
}

 

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
	m= falpham_KdShu20072() / (falpham_KdShu20072() + fbetam_KdShu20072())
	h= falphah_KdShu20072() / (falphah_KdShu20072() + fbetah_KdShu20072())
}

DERIVATIVE states {   
        m' = falpham_KdShu20072() * (1 - m) - fbetam_KdShu20072() * m
        h' = falphah_KdShu20072() * (1 - h) - fbetah_KdShu20072() * h
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}