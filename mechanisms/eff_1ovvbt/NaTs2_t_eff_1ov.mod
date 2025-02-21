:Reference :Colbert and Pan 2002
:comment: took the NaTa and shifted both activation/inactivation by 6 mv
NEURON	{
	SUFFIX NaTs2_t
	USEION na READ ena WRITE ina
	RANGE gNaTs2_tbar, gNaTs2_t, ina
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_NaTs2t_table, betam_NaTs2t_table, alphah_NaTs2t_table, betah_NaTs2t_table, A_1_table, B_1_table
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
	gNaTs2_tbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTs2_t	(S/cm2)
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

	V_table  alpham_NaTs2t_table  betam_NaTs2t_table  alphah_NaTs2t_table  betah_NaTs2t_table  A_1_table  B_1_table  
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

FUNCTION falpham_NaTs2t() { 
VERBATIM
	double alpham_NaTs2t_value;
	alpham_NaTs2t_value = interp4D(_p_alpham_NaTs2t_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_NaTs2t_value);
ENDVERBATIM
	falpham_NaTs2t = alpham_NaTs2t_value
}

FUNCTION fbetam_NaTs2t() { 
VERBATIM
	double betam_NaTs2t_value;
	betam_NaTs2t_value = interp4D(_p_betam_NaTs2t_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_NaTs2t_value);
ENDVERBATIM
	fbetam_NaTs2t = betam_NaTs2t_value
}

FUNCTION falphah_NaTs2t() { 
VERBATIM
	double alphah_NaTs2t_value;
	alphah_NaTs2t_value = interp4D(_p_alphah_NaTs2t_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_NaTs2t_value);
ENDVERBATIM
	falphah_NaTs2t = alphah_NaTs2t_value
}

FUNCTION fbetah_NaTs2t() { 
VERBATIM
	double betah_NaTs2t_value;
	betah_NaTs2t_value = interp4D(_p_betah_NaTs2t_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_NaTs2t_value);
ENDVERBATIM
	fbetah_NaTs2t = betah_NaTs2t_value
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
	gNaTs2_t = gNaTs2_tbar*m*m*m*h
	ina = gNaTs2_t*(Vm-ena)
}

DERIVATIVE states	{
	m' = falpham_NaTs2t() * (1 - m) - fbetam_NaTs2t() * m
	h' = falphah_NaTs2t() * (1 - h) - fbetah_NaTs2t() * h
}

INITIAL{
	update()
	m = falpham_NaTs2t() / (falpham_NaTs2t() + fbetam_NaTs2t())
	h = falphah_NaTs2t() / (falphah_NaTs2t() + fbetah_NaTs2t())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}