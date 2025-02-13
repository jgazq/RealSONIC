:Comment : The persistent component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential

NEURON	{
	SUFFIX K_Pst
	USEION k READ ek WRITE ik
	RANGE gK_Pstbar, gK_Pst, ik
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
	gK_Pstbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Pst	(S/cm2)
	mInf
	mTau
	hInf
	hTau
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)
FUNCTION_TABLE alpham_KPst(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betam_KPst(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_KPst(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betah_KPst(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE A_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (nC/cm2)
FUNCTION_TABLE B_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (nC/cm2)

STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gK_Pst = gK_Pstbar*m*m*h
	ik = gK_Pst*(Vm-ek)
}

DERIVATIVE states	{
	m' = alpham_KPst(A_t, y, a1, b1) * (1 - m) - betam_KPst(A_t, y, a1, b1) * m
	h' = alphah_KPst(A_t, y, a1, b1) * (1 - h) - betah_KPst(A_t, y, a1, b1) * h
}

INITIAL{
	update()
	m = alpham_KPst(A_t, y, a1, b1) / (alpham_KPst(A_t, y, a1, b1) + betam_KPst(A_t, y, a1, b1))
	h = alphah_KPst(A_t, y, a1, b1) / (alphah_KPst(A_t, y, a1, b1) + betah_KPst(A_t, y, a1, b1))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}