:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA2
	USEION ca READ eca WRITE ica
	RANGE gCa_HVAbar, gCa_HVA, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t, q1, f1, q2, f2 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
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
	q1  (nC/cm2)
	f1  (rad)
	q2  (nC/cm2)
	f2  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (mV)
FUNCTION_TABLE A_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (mV)
FUNCTION_TABLE phi_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (rad)
FUNCTION_TABLE alpham_CaHVA2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betam_CaHVA2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE alphah_CaHVA2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betah_CaHVA2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)

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
	m' = alpham_CaHVA2(A_t, y, q1, f1, q2, f2) * (1 - m) - betam_CaHVA2(A_t, y, q1, f1, q2, f2) * m
	h' = alphah_CaHVA2(A_t, y, q1, f1, q2, f2) * (1 - h) - betah_CaHVA2(A_t, y, q1, f1, q2, f2) * h
}

INITIAL{
	update()
	m = alpham_CaHVA2(A_t, y, q1, f1, q2, f2) / (alpham_CaHVA2(A_t, y, q1, f1, q2, f2) + betam_CaHVA2(A_t, y, q1, f1, q2, f2))
	h = alphah_CaHVA2(A_t, y, q1, f1, q2, f2) / (alphah_CaHVA2(A_t, y, q1, f1, q2, f2) + betah_CaHVA2(A_t, y, q1, f1, q2, f2))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}