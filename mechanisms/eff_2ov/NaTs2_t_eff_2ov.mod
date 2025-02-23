:Reference :Colbert and Pan 2002
:comment: took the NaTa and shifted both activation/inactivation by 6 mv
NEURON	{
	SUFFIX NaTs2_t
	USEION na READ ena WRITE ina
	RANGE gNaTs2_tbar, gNaTs2_t, ina
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
	gNaTs2_tbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTs2_t	(S/cm2)
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
FUNCTION_TABLE alpham_NaTs2t(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE betam_NaTs2t(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_NaTs2t(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)
FUNCTION_TABLE betah_NaTs2t(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2), A2(nC/cm2), B2(nC/cm2)) (/ms)

STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gNaTs2_t = gNaTs2_tbar*m*m*m*h
	ina = gNaTs2_t*(Vm-ena)
}

DERIVATIVE states	{
	m' = alpham_NaTs2t(A_t, y, a1, b1, a2, b2) * (1 - m) - betam_NaTs2t(A_t, y, a1, b1, a2, b2) * m
	h' = alphah_NaTs2t(A_t, y, a1, b1, a2, b2) * (1 - h) - betah_NaTs2t(A_t, y, a1, b1, a2, b2) * h
}

INITIAL{
	update()
	m = alpham_NaTs2t(A_t, y, a1, b1, a2, b2) / (alpham_NaTs2t(A_t, y, a1, b1, a2, b2) + betam_NaTs2t(A_t, y, a1, b1, a2, b2))
	h = alphah_NaTs2t(A_t, y, a1, b1, a2, b2) / (alphah_NaTs2t(A_t, y, a1, b1, a2, b2) + betah_NaTs2t(A_t, y, a1, b1, a2, b2))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}