:Comment :
:Reference : :		Characterization of a Shaw-related potassium channel family in rat brain, The EMBO Journal, vol.11, no.7,2473-2486 (1992)

NEURON	{
	SUFFIX SKv3_12
	USEION k READ ek WRITE ik
	RANGE gSKv3_12bar, gSKv3_12, ik 
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
	gSKv3_12bar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gSKv3_12	(S/cm2)
	mInf
	mTau
	A_t  (kPa)
	y
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_SKv312(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_SKv312(A(kPa), Q(nC/cm2)) (/ms)

STATE	{ 
	m
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gSKv3_12 = gSKv3_12bar*m
	ik = gSKv3_12*(Vm-ek)
}

DERIVATIVE states	{
	m' = alpham_SKv312(A_t, y) * (1 - m) - betam_SKv312(A_t, y) * m
}

INITIAL{
	update()
	m = alpham_SKv312(A_t, y) / (alpham_SKv312(A_t, y) + betam_SKv312(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}