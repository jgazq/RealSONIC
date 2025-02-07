:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA
	USEION ca READ eca WRITE ica
	RANGE gCa_HVAbar, gCa_HVA, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t, q1, f1 : section specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
	PI = (pi) (1)
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
	q1 (nC/cm2)
	f1 (rad)
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
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE psi1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE alpham_CaHVA(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betam_CaHVA(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE alphah_CaHVA(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betah_CaHVA(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)

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
	m' = alpham_CaHVA(A_t, y, q1, f1) * (1 - m) - betam_CaHVA(A_t, y, q1, f1) * m
	h' = alphah_CaHVA(A_t, y, q1, f1) * (1 - h) - betah_CaHVA(A_t, y, q1, f1) * h
}

INITIAL{
	update()
	m = alpham_CaHVA(A_t, y, q1, f1) / (alpham_CaHVA(A_t, y, q1, f1) + betam_CaHVA(A_t, y, q1, f1))
	h = alphah_CaHVA(A_t, y, q1, f1) / (alphah_CaHVA(A_t, y, q1, f1) + betah_CaHVA(A_t, y, q1, f1))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}