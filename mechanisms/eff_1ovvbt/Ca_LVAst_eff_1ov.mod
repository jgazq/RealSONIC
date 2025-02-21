:Comment : LVA ca channel. Note: mtau is an approximation from the plots
:Reference : :		Avery and Johnston 1996, tau from Randall 1997
:Comment: shifted by -10 mv to correct for junction potential
NEURON	{
	SUFFIX Ca_LVAst
	USEION ca READ eca WRITE ica
	RANGE gCa_LVAstbar, gCa_LVAst, ica
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)

	POINTER V_table, alpham_CaLVAst_table, betam_CaLVAst_table, alphah_CaLVAst_table, betah_CaLVAst_table, A_1_table, B_1_table
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
	gCa_LVAstbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_LVAst	(S/cm2)
	mInf
	mTau
	hInf
	hTau
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)

	V_table  alpham_CaLVAst_table  betam_CaLVAst_table  alphah_CaLVAst_table  betah_CaLVAst_table  A_1_table  B_1_table  
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

FUNCTION falpham_CaLVAst() { 
VERBATIM
	double alpham_CaLVAst_value;
	alpham_CaLVAst_value = interp4D(_p_alpham_CaLVAst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alpham_CaLVAst_value);
ENDVERBATIM
	falpham_CaLVAst = alpham_CaLVAst_value
}

FUNCTION fbetam_CaLVAst() { 
VERBATIM
	double betam_CaLVAst_value;
	betam_CaLVAst_value = interp4D(_p_betam_CaLVAst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betam_CaLVAst_value);
ENDVERBATIM
	fbetam_CaLVAst = betam_CaLVAst_value
}

FUNCTION falphah_CaLVAst() { 
VERBATIM
	double alphah_CaLVAst_value;
	alphah_CaLVAst_value = interp4D(_p_alphah_CaLVAst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(alphah_CaLVAst_value);
ENDVERBATIM
	falphah_CaLVAst = alphah_CaLVAst_value
}

FUNCTION fbetah_CaLVAst() { 
VERBATIM
	double betah_CaLVAst_value;
	betah_CaLVAst_value = interp4D(_p_betah_CaLVAst_table, _p_A_arr, _p_Q_arr, _p_A1_arr, _p_B1_arr, A_s, Q_s, A1_s, B1_s, A_t, v, a1, b1);
	return(betah_CaLVAst_value);
ENDVERBATIM
	fbetah_CaLVAst = betah_CaLVAst_value
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
	gCa_LVAst = gCa_LVAstbar*m*m*h
	ica = gCa_LVAst*(Vm-eca)
}

DERIVATIVE states	{
	m' = falpham_CaLVAst() * (1 - m) - fbetam_CaLVAst() * m
	h' = falphah_CaLVAst() * (1 - h) - fbetah_CaLVAst() * h
}

INITIAL{
	update()
	m = falpham_CaLVAst() / (falpham_CaLVAst() + fbetam_CaLVAst())
	h = falphah_CaLVAst() / (falphah_CaLVAst() + fbetah_CaLVAst())
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}