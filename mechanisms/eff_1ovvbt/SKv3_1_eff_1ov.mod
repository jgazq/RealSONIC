:Comment :
:Reference : :		Characterization of a Shaw-related potassium channel family in rat brain, The EMBO Journal, vol.11, no.7,2473-2486 (1992)

NEURON	{
	SUFFIX SKv3_1
	USEION k READ ek WRITE ik
	RANGE gSKv3_1bar, gSKv3_1, ik 
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_SKv31_table, betam_SKv31_table, A_1_table, B_1_table
	RANGE V_val, alpham_SKv31_val, betam_SKv31_val, A_1_val, B_1_val
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
	gSKv3_1bar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gSKv3_1	(S/cm2)
	mInf
	mTau
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)

	V_table  alpham_SKv31_table  betam_SKv31_table  A_1_table  B_1_table  
	V_val (mV)  alpham_SKv31_val (/ms)  betam_SKv31_val (/ms)  A_1_val (nC/cm2)  B_1_val (nC/cm2)  
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

FUNCTION falpham_SKv31() { 
VERBATIM
	double alpham_SKv31_value;
	alpham_SKv31_value = interp4D(_p_alpham_SKv31_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_SKv31_value);
ENDVERBATIM
	falpham_SKv31 = alpham_SKv31_value
}

FUNCTION fbetam_SKv31() { 
VERBATIM
	double betam_SKv31_value;
	betam_SKv31_value = interp4D(_p_betam_SKv31_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_SKv31_value);
ENDVERBATIM
	fbetam_SKv31 = betam_SKv31_value
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
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gSKv3_1 = gSKv3_1bar*m
	ik = gSKv3_1*(Vm-ek)
}

DERIVATIVE states	{
	m' = falpham_SKv31() * (1 - m) - fbetam_SKv31() * m
}

INITIAL{
	update()
	m = falpham_SKv31() / (falpham_SKv31() + fbetam_SKv31())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}