:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
NEURON	{
	SUFFIX K_Tst2
	USEION k READ ek WRITE ik
	RANGE gK_Tstbar, gK_Tst, ik
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
	gK_Tstbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Tst	(S/cm2)
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

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_KTst2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_KTst2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_KTst2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betah_KTst2(A(kPa), Q(nC/cm2)) (/ms)

STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gK_Tst = gK_Tstbar*(m^4)*h
	ik = gK_Tst*(Vm-ek)
}

DERIVATIVE states	{
	m' = alpham_KTst2(A_t, y) * (1 - m) - betam_KTst2(A_t, y) * m
	h' = alphah_KTst2(A_t, y) * (1 - h) - betah_KTst2(A_t, y) * h
}

INITIAL{
	update()
	m = alpham_KTst2(A_t, y) / (alpham_KTst2(A_t, y) + betam_KTst2(A_t, y))
	h = alphah_KTst2(A_t, y) / (alphah_KTst2(A_t, y) + betah_KTst2(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}