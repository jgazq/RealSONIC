:Comment : The persistent component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential

NEURON	{
	SUFFIX K_Pst2
	USEION k READ ek WRITE ik
	RANGE gK_Pstbar, gK_Pst, ik
	RANGE Adrive, Vm, y, Fdrive, A_t : section specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
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
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_KPst(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_KPst(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_KPst(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betah_KPst(A(kPa), Q(nC/cm2)) (/ms)

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
	m' = alpham_KPst(A_t, y) * (1 - m) - betam_KPst(A_t, y) * m
	h' = alphah_KPst(A_t, y) * (1 - h) - betah_KPst(A_t, y) * h
}

INITIAL{
	update()
	m = alpham_KPst(A_t, y) / (alpham_KPst(A_t, y) + betam_KPst(A_t, y))
	h = alphah_KPst(A_t, y) / (alphah_KPst(A_t, y) + betah_KPst(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}