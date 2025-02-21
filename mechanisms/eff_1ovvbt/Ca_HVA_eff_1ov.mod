:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca WRITE ica
	RANGE gCa_HVAbar, gCa_HVA, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_CaHVA_table, betam_CaHVA_table, alphah_CaHVA_table, betah_CaHVA_table, A_1_table, B_1_table
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
	gCa_HVAbar = 0.00001 (S/cm2) 
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

	V_table  alpham_CaHVA_table  betam_CaHVA_table  alphah_CaHVA_table  betah_CaHVA_table  A_1_table  B_1_table  
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

FUNCTION falpham_CaHVA() { 
VERBATIM
	double alpham_CaHVA_value;
	alpham_CaHVA_value = interp4D(_p_alpham_CaHVA_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_CaHVA_value);
ENDVERBATIM
	falpham_CaHVA = alpham_CaHVA_value
}

FUNCTION fbetam_CaHVA() { 
VERBATIM
	double betam_CaHVA_value;
	betam_CaHVA_value = interp4D(_p_betam_CaHVA_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_CaHVA_value);
ENDVERBATIM
	fbetam_CaHVA = betam_CaHVA_value
}

FUNCTION falphah_CaHVA() { 
VERBATIM
	double alphah_CaHVA_value;
	alphah_CaHVA_value = interp4D(_p_alphah_CaHVA_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_CaHVA_value);
ENDVERBATIM
	falphah_CaHVA = alphah_CaHVA_value
}

FUNCTION fbetah_CaHVA() { 
VERBATIM
	double betah_CaHVA_value;
	betah_CaHVA_value = interp4D(_p_betah_CaHVA_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_CaHVA_value);
ENDVERBATIM
	fbetah_CaHVA = betah_CaHVA_value
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
	gCa = gCa_HVAbar*m*m*h
	ica = gCa*(Vm-eca)
}

DERIVATIVE states	{
	m' = falpham_CaHVA() * (1 - m) - fbetam_CaHVA() * m
	h' = falphah_CaHVA() * (1 - h) - fbetah_CaHVA() * h
}

INITIAL{
	update()
	m = falpham_CaHVA() / (falpham_CaHVA() + fbetam_CaHVA())
	h = falphah_CaHVA() / (falphah_CaHVA() + fbetah_CaHVA())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}