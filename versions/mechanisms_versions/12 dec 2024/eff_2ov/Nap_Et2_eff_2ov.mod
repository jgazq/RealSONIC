:Comment : mtau deduced from text (said to be 6 times faster than for NaTa)
:Comment : so I used the equations from NaT and multiplied by 6
:Reference : Modeled according to kinetics derived from Magistretti & Alonso 1999
NEURON	{
	SUFFIX Nap_Et2
	USEION na READ ena WRITE ina
	RANGE gNap_Et2bar, gNap_Et2, ina
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
	gNap_Et2bar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNap_Et2	(S/cm2)
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
FUNCTION_TABLE alpham_NapEt2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betam_NapEt2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE alphah_NapEt2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betah_NapEt2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)

STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gNap_Et2 = gNap_Et2bar*m*m*m*h
	ina = gNap_Et2*(Vm-ena)
}

DERIVATIVE states	{
	m' = alpham_NapEt2(A_t, y) * (1 - m) - betam_NapEt2(A_t, y) * m
	h' = alphah_NapEt2(A_t, y) * (1 - h) - betah_NapEt2(A_t, y) * h
}

INITIAL{
	update()
	m = alpham_NapEt2(A_t, y) / (alpham_NapEt2(A_t, y) + betam_NapEt2(A_t, y))
	h = alphah_NapEt2(A_t, y) / (alphah_NapEt2(A_t, y) + betah_NapEt2(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}