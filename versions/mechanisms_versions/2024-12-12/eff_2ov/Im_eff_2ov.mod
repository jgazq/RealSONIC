:Reference : :		Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones
NEURON	{
	SUFFIX Im
	USEION k READ ek WRITE ik
	RANGE gImbar, gIm, ik
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
	gImbar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gIm	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
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
FUNCTION_TABLE alpham_Im(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betam_Im(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)

STATE	{ 
	m
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gIm = gImbar*m
	ik = gIm*(Vm-ek)
}

DERIVATIVE states	{
	m' = alpham_Im(A_t, y) * (1 - m) - betam_Im(A_t, y) * m
}

INITIAL{
	update()
	m = alpham_Im(A_t, y) / (alpham_Im(A_t, y) + betam_Im(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}