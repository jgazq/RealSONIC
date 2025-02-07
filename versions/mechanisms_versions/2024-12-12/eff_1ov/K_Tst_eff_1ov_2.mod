:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
NEURON	{
	SUFFIX K_Tst2
	USEION k READ ek WRITE ik
	RANGE gK_Tst2bar, gK_Tst2, ik
	RANGE Adrive, Vm, y, Fdrive, A_t, q1, f1 : section (even segment) specific
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
	gK_Tst2bar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Tst2	(S/cm2)
	mInf
	mTau
	hInf
	hTau
	A_t  (kPa)
	y
	q1  (nC/cm2)
	f1  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE A_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (mV)
FUNCTION_TABLE phi_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (rad)
FUNCTION_TABLE alpham_KTst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betam_KTst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE alphah_KTst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)
FUNCTION_TABLE betah_KTst2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad)) (/ms)

STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gK_Tst2 = gK_Tst2bar*(m^4)*h
	ik = gK_Tst2*(Vm-ek)
}

DERIVATIVE states	{
	m' = alpham_KTst2(A_t, y, q1, f1) * (1 - m) - betam_KTst2(A_t, y, q1, f1) * m
	h' = alphah_KTst2(A_t, y, q1, f1) * (1 - h) - betah_KTst2(A_t, y, q1, f1) * h
}

INITIAL{
	update()
	m = alpham_KTst2(A_t, y, q1, f1) / (alpham_KTst2(A_t, y, q1, f1) + betam_KTst2(A_t, y, q1, f1))
	h = alphah_KTst2(A_t, y, q1, f1) / (alphah_KTst2(A_t, y, q1, f1) + betah_KTst2(A_t, y, q1, f1))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}