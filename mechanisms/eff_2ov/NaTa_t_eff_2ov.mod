:Reference :Colbert and Pan 2002
NEURON	{
	SUFFIX NaTa_t
	USEION na READ ena WRITE ina
	RANGE gNaTa_tbar, gNaTa_t, ina
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
	gNaTa_tbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTa_t	(S/cm2)
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
FUNCTION_TABLE alpham_NaTat(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betam_NaTat(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE alphah_NaTat(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betah_NaTat(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)

STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gNaTa_t = gNaTa_tbar*m*m*m*h
	ina = gNaTa_t*(Vm-ena)
}

DERIVATIVE states	{
	m' = alpham_NaTat(A_t, y, q1, f1, q2, f2) * (1 - m) - betam_NaTat(A_t, y, q1, f1, q2, f2) * m
	h' = alphah_NaTat(A_t, y, q1, f1, q2, f2) * (1 - h) - betah_NaTat(A_t, y, q1, f1, q2, f2) * h
}

INITIAL{
	update()
	m = alpham_NaTat(A_t, y, q1, f1, q2, f2) / (alpham_NaTat(A_t, y, q1, f1, q2, f2) + betam_NaTat(A_t, y, q1, f1, q2, f2))
	h = alphah_NaTat(A_t, y, q1, f1, q2, f2) / (alphah_NaTat(A_t, y, q1, f1, q2, f2) + betah_NaTat(A_t, y, q1, f1, q2, f2))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}