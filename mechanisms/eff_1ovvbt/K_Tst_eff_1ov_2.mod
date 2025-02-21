:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
NEURON	{
	SUFFIX K_Tst2
	USEION k READ ek WRITE ik
	RANGE gK_Tstbar, gK_Tst, ik
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_KTst2_table, betam_KTst2_table, alphah_KTst2_table, betah_KTst2_table, A_1_table, B_1_table
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
	gK_Tstbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Tst	(S/cm2)
	mInf
	mTau
	hInf
	hTau
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)

	V_table  alpham_KTst2_table  betam_KTst2_table  alphah_KTst2_table  betah_KTst2_table  A_1_table  B_1_table  
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

FUNCTION falpham_KTst2() { 
VERBATIM
	double alpham_KTst2_value;
	alpham_KTst2_value = interp4D(_p_alpham_KTst2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_KTst2_value);
ENDVERBATIM
	falpham_KTst2 = alpham_KTst2_value
}

FUNCTION fbetam_KTst2() { 
VERBATIM
	double betam_KTst2_value;
	betam_KTst2_value = interp4D(_p_betam_KTst2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_KTst2_value);
ENDVERBATIM
	fbetam_KTst2 = betam_KTst2_value
}

FUNCTION falphah_KTst2() { 
VERBATIM
	double alphah_KTst2_value;
	alphah_KTst2_value = interp4D(_p_alphah_KTst2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_KTst2_value);
ENDVERBATIM
	falphah_KTst2 = alphah_KTst2_value
}

FUNCTION fbetah_KTst2() { 
VERBATIM
	double betah_KTst2_value;
	betah_KTst2_value = interp4D(_p_betah_KTst2_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_KTst2_value);
ENDVERBATIM
	fbetah_KTst2 = betah_KTst2_value
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
	gK_Tst = gK_Tstbar*(m^4)*h
	ik = gK_Tst*(Vm-ek)
}

DERIVATIVE states	{
	m' = falpham_KTst2() * (1 - m) - fbetam_KTst2() * m
	h' = falphah_KTst2() * (1 - h) - fbetah_KTst2() * h
}

INITIAL{
	update()
	m = falpham_KTst2() / (falpham_KTst2() + fbetam_KTst2())
	h = falphah_KTst2() / (falphah_KTst2() + fbetah_KTst2())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}