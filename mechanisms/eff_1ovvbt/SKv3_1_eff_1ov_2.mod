:Comment :
:Reference : :		Characterization of a Shaw-related potassium channel family in rat brain, The EMBO Journal, vol.11, no.7,2473-2486 (1992)

NEURON	{
	SUFFIX SKv3_12
	USEION k READ ek WRITE ik
	RANGE gSKv3_1bar, gSKv3_1, ik 
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_SKv312_table, betam_SKv312_table, A_1_table, B_1_table
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

	V_table  alpham_SKv312_table  betam_SKv312_table  A_1_table  B_1_table  
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

FUNCTION falpham_SKv312() { 
VERBATIM
	double alpham_SKv312_value;
	alpham_SKv312_value = interp4D(_p_alpham_SKv312_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_SKv312_value);
ENDVERBATIM
	falpham_SKv312 = alpham_SKv312_value
}

FUNCTION fbetam_SKv312() { 
VERBATIM
	double betam_SKv312_value;
	betam_SKv312_value = interp4D(_p_betam_SKv312_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_SKv312_value);
ENDVERBATIM
	fbetam_SKv312 = betam_SKv312_value
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
	m' = falpham_SKv312() * (1 - m) - fbetam_SKv312() * m
}

INITIAL{
	update()
	m = falpham_SKv312() / (falpham_SKv312() + fbetam_SKv312())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}