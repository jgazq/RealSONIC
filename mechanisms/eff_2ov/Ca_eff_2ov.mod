:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca
	USEION ca READ eca WRITE ica
	RANGE gCabar, gCa, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1, a2, b2 : section (even segment) specific
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
	a1  (nC/cm2)
	b1  (rad)
	a2  (nC/cm2)
	b2  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (mV)
FUNCTION_TABLE alpham_Ca(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE betam_Ca(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_Ca(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE betah_Ca(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)

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
	m' = alpham_Ca(A_t, y, a1, b1, a2, b2) * (1 - m) - betam_Ca(A_t, y, a1, b1, a2, b2) * m
	h' = alphah_Ca(A_t, y, a1, b1, a2, b2) * (1 - h) - betah_Ca(A_t, y, a1, b1, a2, b2) * h
}

INITIAL{
	update()
	m = alpham_Ca(A_t, y, a1, b1, a2, b2) / (alpham_Ca(A_t, y, a1, b1, a2, b2) + betam_Ca(A_t, y, a1, b1, a2, b2))
	h = alphah_Ca(A_t, y, a1, b1, a2, b2) / (alphah_Ca(A_t, y, a1, b1, a2, b2) + betah_Ca(A_t, y, a1, b1, a2, b2))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}