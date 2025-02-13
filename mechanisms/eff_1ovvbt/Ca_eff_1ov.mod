:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca
	USEION ca READ eca WRITE ica
	RANGE gCabar, gCa, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_Ca_table, betam_Ca_table, alphah_Ca_table, betah_Ca_table, A_1_table, B_1_table
	RANGE V_val, alpham_Ca_val, betam_Ca_val, alphah_Ca_val, betah_Ca_val, A_1_val, B_1_val
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
	gCabar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa	(S/cm2)
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

	V_table  alpham_Ca_table  betam_Ca_table  alphah_Ca_table  betah_Ca_table  A_1_table  B_1_table  
	V_val (mV)  alpham_Ca_val (/ms)  betam_Ca_val (/ms)  alphah_Ca_val (/ms)  betah_Ca_val (/ms)  A_1_val (nC/cm2)  B_1_val (nC/cm2)  
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

FUNCTION falpham_Ca() { 
VERBATIM
	double alpham_Ca_value;
	alpham_Ca_value = interp4D(_p_alpham_Ca_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_Ca_value);
ENDVERBATIM
	falpham_Ca = alpham_Ca_value
}

FUNCTION fbetam_Ca() { 
VERBATIM
	double betam_Ca_value;
	betam_Ca_value = interp4D(_p_betam_Ca_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_Ca_value);
ENDVERBATIM
	fbetam_Ca = betam_Ca_value
}

FUNCTION falphah_Ca() { 
VERBATIM
	double alphah_Ca_value;
	alphah_Ca_value = interp4D(_p_alphah_Ca_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_Ca_value);
ENDVERBATIM
	falphah_Ca = alphah_Ca_value
}

FUNCTION fbetah_Ca() { 
VERBATIM
	double betah_Ca_value;
	betah_Ca_value = interp4D(_p_betah_Ca_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_Ca_value);
ENDVERBATIM
	fbetah_Ca = betah_Ca_value
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
	gCa = gCabar*m*m*h
	ica = gCa*(Vm-eca)
}

DERIVATIVE states	{
	m' = falpham_Ca() * (1 - m) - fbetam_Ca() * m
	h' = falphah_Ca() * (1 - h) - fbetah_Ca() * h
}

INITIAL{
	update()
	m = falpham_Ca() / (falpham_Ca() + fbetam_Ca())
	h = falphah_Ca() / (falphah_Ca() + fbetah_Ca())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}