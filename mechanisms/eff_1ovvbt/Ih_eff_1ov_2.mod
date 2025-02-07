:Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih2
	NONSPECIFIC_CURRENT ihcn
	RANGE gIhbar, gIh, ihcn 
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
	gIhbar = 0.00001 (S/cm2) 
	ehcn =  -45.0 (mV)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ihcn	(mA/cm2)
	gIh	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	A_t  (kPa)
	y
	a1  (nC/cm2)
	b1  (rad)
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_Ih2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_Ih2(A(kPa), Q(nC/cm2)) (/ms)

STATE	{ 
	m
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gIh = gIhbar*m
	ihcn = gIh*(Vm-ehcn)
}

DERIVATIVE states	{
	m' = alpham_Ih2(A_t, y) * (1 - m) - betam_Ih2(A_t, y) * m
}

INITIAL{
	update()
	m = alpham_Ih2(A_t, y) / (alpham_Ih2(A_t, y) + betam_Ih2(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}