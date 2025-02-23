:Comment : mtau deduced from text (said to be 6 times faster than for NaTa)
:Comment : so I used the equations from NaT and multiplied by 6
:Reference : Modeled according to kinetics derived from Magistretti & Alonso 1999
NEURON	{
	SUFFIX Nap_Et22
	USEION na READ ena WRITE ina
	RANGE gNap_Et2bar, gNap_Et2, ina
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_NapEt22_table, betam_NapEt22_table, alphah_NapEt22_table, betah_NapEt22_table, A_1_table, B_1_table
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
	gNap_Et2bar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNap_Et2	(S/cm2)
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

	V_table  alpham_NapEt22_table  betam_NapEt22_table  alphah_NapEt22_table  betah_NapEt22_table  A_1_table  B_1_table  
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

FUNCTION falpham_NapEt22() { 
VERBATIM
	double alpham_NapEt22_value;
	alpham_NapEt22_value = interp4D(_p_alpham_NapEt22_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_NapEt22_value);
ENDVERBATIM
	falpham_NapEt22 = alpham_NapEt22_value
}

FUNCTION fbetam_NapEt22() { 
VERBATIM
	double betam_NapEt22_value;
	betam_NapEt22_value = interp4D(_p_betam_NapEt22_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_NapEt22_value);
ENDVERBATIM
	fbetam_NapEt22 = betam_NapEt22_value
}

FUNCTION falphah_NapEt22() { 
VERBATIM
	double alphah_NapEt22_value;
	alphah_NapEt22_value = interp4D(_p_alphah_NapEt22_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_NapEt22_value);
ENDVERBATIM
	falphah_NapEt22 = alphah_NapEt22_value
}

FUNCTION fbetah_NapEt22() { 
VERBATIM
	double betah_NapEt22_value;
	betah_NapEt22_value = interp4D(_p_betah_NapEt22_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_NapEt22_value);
ENDVERBATIM
	fbetah_NapEt22 = betah_NapEt22_value
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
	gNap_Et2 = gNap_Et2bar*m*m*m*h
	ina = gNap_Et2*(Vm-ena)
}

DERIVATIVE states	{
	m' = falpham_NapEt22() * (1 - m) - fbetam_NapEt22() * m
	h' = falphah_NapEt22() * (1 - h) - fbetah_NapEt22() * h
}

INITIAL{
	update()
	m = falpham_NapEt22() / (falpham_NapEt22() + fbetam_NapEt22())
	h = falphah_NapEt22() / (falphah_NapEt22() + fbetah_NapEt22())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}