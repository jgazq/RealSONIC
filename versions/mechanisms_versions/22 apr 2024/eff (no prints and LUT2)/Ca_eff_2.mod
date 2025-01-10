:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca2
	USEION ca READ eca WRITE ica
	RANGE gCabar, gCa, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t : section specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
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
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_Ca2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_Ca2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_Ca2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betah_Ca2(A(kPa), Q(nC/cm2)) (/ms)

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
	m' = alpham_Ca2(A_t, y) * (1 - m) - betam_Ca2(A_t, y) * m
	h' = alphah_Ca2(A_t, y) * (1 - h) - betah_Ca2(A_t, y) * h
}

INITIAL{
	update()
	m = alpham_Ca2(A_t, y) / (alpham_Ca2(A_t, y) + betam_Ca2(A_t, y))
	h = alphah_Ca2(A_t, y) / (alphah_Ca2(A_t, y) + betah_Ca2(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}