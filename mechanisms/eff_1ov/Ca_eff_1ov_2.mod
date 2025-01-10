:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca2
	USEION ca READ eca WRITE ica
	RANGE gCabar, gCa, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t, q1, f1 : section (even segment) specific
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
	q1  (nC/cm2)
	f1  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE A_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE phi_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (rad)
FUNCTION_TABLE alpham_Ca2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betam_Ca2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE alphah_Ca2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betah_Ca2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)

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
	m' = alpham_Ca2(A_t, y, q1, f1) * (1 - m) - betam_Ca2(A_t, y, q1, f1) * m
	h' = alphah_Ca2(A_t, y, q1, f1) * (1 - h) - betah_Ca2(A_t, y, q1, f1) * h
}

INITIAL{
	update()
	m = alpham_Ca2(A_t, y, q1, f1) / (alpham_Ca2(A_t, y, q1, f1) + betam_Ca2(A_t, y, q1, f1))
	h = alphah_Ca2(A_t, y, q1, f1) / (alphah_Ca2(A_t, y, q1, f1) + betah_Ca2(A_t, y, q1, f1))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}