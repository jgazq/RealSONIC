:Comment : mtau deduced from text (said to be 6 times faster than for NaTa)
:Comment : so I used the equations from NaT and multiplied by 6
:Reference : Modeled according to kinetics derived from Magistretti & Alonso 1999
NEURON	{
	SUFFIX Nap_Et22
	USEION na READ ena WRITE ina
	RANGE gNap_Et2bar, gNap_Et2, ina
	RANGE Adrive, Vm, y, Fdrive, A_t, a1, b1 : section (even segment) specific
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
	a1  (nC/cm2)
	b1  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)
FUNCTION_TABLE alpham_NapEt22(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betam_NapEt22(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_NapEt22(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betah_NapEt22(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE A_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)
FUNCTION_TABLE B_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)

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
	m' = alpham_NapEt22(A_t, y, a1, b1) * (1 - m) - betam_NapEt22(A_t, y, a1, b1) * m
	h' = alphah_NapEt22(A_t, y, a1, b1) * (1 - h) - betah_NapEt22(A_t, y, a1, b1) * h
}

INITIAL{
	update()
	m = alpham_NapEt22(A_t, y, a1, b1) / (alpham_NapEt22(A_t, y, a1, b1) + betam_NapEt22(A_t, y, a1, b1))
	h = alphah_NapEt22(A_t, y, a1, b1) / (alphah_NapEt22(A_t, y, a1, b1) + betah_NapEt22(A_t, y, a1, b1))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}