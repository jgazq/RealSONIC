:Comment : The persistent component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential

NEURON	{
	SUFFIX K_Pst2
	USEION k READ ek WRITE ik
	RANGE gK_Pstbar, gK_Pst, ik
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
	q1  (nC/cm2)
	f1  (rad)
	q2  (nC/cm2)
	f2  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (mV)
FUNCTION_TABLE A_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (mV)
FUNCTION_TABLE phi_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (rad)
FUNCTION_TABLE alpham_KPst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betam_KPst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE alphah_KPst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)
FUNCTION_TABLE betah_KPst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (/ms)

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
	m' = alpham_KPst2(A_t, y, q1, f1, q2, f2) * (1 - m) - betam_KPst2(A_t, y, q1, f1, q2, f2) * m
	h' = alphah_KPst2(A_t, y, q1, f1, q2, f2) * (1 - h) - betah_KPst2(A_t, y, q1, f1, q2, f2) * h
}

INITIAL{
	update()
	m = alpham_KPst2(A_t, y, q1, f1, q2, f2) / (alpham_KPst2(A_t, y, q1, f1, q2, f2) + betam_KPst2(A_t, y, q1, f1, q2, f2))
	h = alphah_KPst2(A_t, y, q1, f1, q2, f2) / (alphah_KPst2(A_t, y, q1, f1, q2, f2) + betah_KPst2(A_t, y, q1, f1, q2, f2))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}