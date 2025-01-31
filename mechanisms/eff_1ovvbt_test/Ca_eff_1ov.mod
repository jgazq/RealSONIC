:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca
	USEION ca READ eca WRITE ica
	RANGE gCabar, gCa, ica 
	RANGE Adrive, Vm, y, Fdrive, A_t, q1, f1 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
	POINTER V
	RANGE V_val
	POINTER A_V1, phi_V1
	POINTER alpham_Ca, betam_Ca, alphah_Ca, betah_Ca
	RANGE alpham_Ca_val, betam_Ca_val, alphah_Ca_val, betah_Ca_val
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
	V_val (mV)
	alpham_Ca_val (/ms)
	betam_Ca_val (/ms)
	alphah_Ca_val (/ms)
	betah_Ca_val (/ms)
}

INCLUDE "update.inc"

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
	VERBATIM
	alpham_Ca_val = interp(alpham_Ca)
	betam_Ca_val = interp(betam_Ca)
	alphah_Ca_val = interp(alphah_Ca)
	betah_Ca_val = interp(betah_Ca)
	ENDVERBATIM
	m' = alpham_Ca_val * (1 - m) - betam_Ca_val * m
	h' = alphah_Ca_val * (1 - h) - betah_Ca_val * h
}

INITIAL{
	update()
	VERBATIM
	alpham_Ca_val = interp(alpham_Ca)
	betam_Ca_val = interp(betam_Ca)
	alphah_Ca_val = interp(alphah_Ca)
	betah_Ca_val = interp(betah_Ca)
	ENDVERBATIM
	m = alpham_Ca_val / (alpham_Ca_val + betam_Ca_val)
	h = alphah_Ca_val / (alphah_Ca_val + betah_Ca_val)
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}