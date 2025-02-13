:Comment : The persistent component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential

NEURON	{
	SUFFIX K_Pst
	USEION k READ ek WRITE ik
	RANGE gK_Pstbar, gK_Pst, ik
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_KPst_table, betam_KPst_table, alphah_KPst_table, betah_KPst_table, A_1_table, B_1_table
	RANGE V_val, alpham_KPst_val, betam_KPst_val, alphah_KPst_val, betah_KPst_val, A_1_val, B_1_val
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
	gK_Pstbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Pst	(S/cm2)
	mInf
	mTau
	hInf
	hTau
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)

	V_table  alpham_KPst_table  betam_KPst_table  alphah_KPst_table  betah_KPst_table  A_1_table  B_1_table  
	V_val (mV)  alpham_KPst_val (/ms)  betam_KPst_val (/ms)  alphah_KPst_val (/ms)  betah_KPst_val (/ms)  A_1_val (nC/cm2)  B_1_val (nC/cm2)  
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

FUNCTION falpham_KPst() { 
VERBATIM
	double alpham_KPst_value;
	alpham_KPst_value = interp4D(_p_alpham_KPst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_KPst_value);
ENDVERBATIM
	falpham_KPst = alpham_KPst_value
}

FUNCTION fbetam_KPst() { 
VERBATIM
	double betam_KPst_value;
	betam_KPst_value = interp4D(_p_betam_KPst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_KPst_value);
ENDVERBATIM
	fbetam_KPst = betam_KPst_value
}

FUNCTION falphah_KPst() { 
VERBATIM
	double alphah_KPst_value;
	alphah_KPst_value = interp4D(_p_alphah_KPst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_KPst_value);
ENDVERBATIM
	falphah_KPst = alphah_KPst_value
}

FUNCTION fbetah_KPst() { 
VERBATIM
	double betah_KPst_value;
	betah_KPst_value = interp4D(_p_betah_KPst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_KPst_value);
ENDVERBATIM
	fbetah_KPst = betah_KPst_value
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
	gK_Pst = gK_Pstbar*m*m*h
	ik = gK_Pst*(Vm-ek)
}

DERIVATIVE states	{
	m' = falpham_KPst() * (1 - m) - fbetam_KPst() * m
	h' = falphah_KPst() * (1 - h) - fbetah_KPst() * h
}

INITIAL{
	update()
	m = falpham_KPst() / (falpham_KPst() + fbetam_KPst())
	h = falphah_KPst() / (falphah_KPst() + fbetah_KPst())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}