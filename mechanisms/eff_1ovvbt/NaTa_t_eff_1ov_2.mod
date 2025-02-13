:Reference :Colbert and Pan 2002
NEURON	{
	SUFFIX NaTa_t2
	USEION na READ ena WRITE ina
	RANGE gNaTa_tbar, gNaTa_t, ina
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_NaTat2_table, betam_NaTat2_table, alphah_NaTat2_table, betah_NaTat2_table, A_1_table, B_1_table
	RANGE V_val, alpham_NaTat2_val, betam_NaTat2_val, alphah_NaTat2_val, betah_NaTat2_val, A_1_val, B_1_val
	POINTER A_arr, Q_arr, A1_arr, B1_arr
	RANGE A_s, Q_s, A1_s, B1_s
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
	PI = (pi) (1) :in order to use the constant pi
}

PARAMETER	{
	stimon       : Stimulation state
	Fdrive (kHz) : Stimulation frequency
	Adrive (kPa) : Stimulation amplitude
	detailed     : Simulation type
	gNaTa_tbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTa_t	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	hInf
	hTau
	hAlpha
	hBeta
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)

	V_table  alpham_NaTat2_table  betam_NaTat2_table  alphah_NaTat2_table  betah_NaTat2_table  A_1_table  B_1_table  
	V_val (mV)  alpham_NaTat2_val (/ms)  betam_NaTat2_val (/ms)  alphah_NaTat2_val (/ms)  betah_NaTat2_val (/ms)  A_1_val (nC/cm2)  B_1_val (nC/cm2)  
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

FUNCTION falpham_NaTat2() { 
VERBATIM
	double alpham_NaTat2_value;
	alpham_NaTat2_value = interp4D(_p_alpham_NaTat2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_NaTat2_value);
ENDVERBATIM
	falpham_NaTat2 = alpham_NaTat2_value
}

FUNCTION fbetam_NaTat2() { 
VERBATIM
	double betam_NaTat2_value;
	betam_NaTat2_value = interp4D(_p_betam_NaTat2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_NaTat2_value);
ENDVERBATIM
	fbetam_NaTat2 = betam_NaTat2_value
}

FUNCTION falphah_NaTat2() { 
VERBATIM
	double alphah_NaTat2_value;
	alphah_NaTat2_value = interp4D(_p_alphah_NaTat2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_NaTat2_value);
ENDVERBATIM
	falphah_NaTat2 = alphah_NaTat2_value
}

FUNCTION fbetah_NaTat2() { 
VERBATIM
	double betah_NaTat2_value;
	betah_NaTat2_value = interp4D(_p_betah_NaTat2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_NaTat2_value);
ENDVERBATIM
	fbetah_NaTat2 = betah_NaTat2_value
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


STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gNaTa_t = gNaTa_tbar*m*m*m*h
	ina = gNaTa_t*(Vm-ena)
}

DERIVATIVE states	{
	m' = falpham_NaTat2() * (1 - m) - fbetam_NaTat2() * m
	h' = falphah_NaTat2() * (1 - h) - fbetah_NaTat2() * h
}

INITIAL{
	update()
	m = falpham_NaTat2() / (falpham_NaTat2() + fbetam_NaTat2())
	h = falphah_NaTat2() / (falphah_NaTat2() + fbetah_NaTat2())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}