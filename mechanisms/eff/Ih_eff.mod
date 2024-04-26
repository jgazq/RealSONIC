:Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih
	NONSPECIFIC_CURRENT ihcn
	RANGE gIhbar, gIh, ihcn 
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
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_Ih(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_Ih(A(kPa), Q(nC/cm2)) (/ms)

STATE	{ 
	m
}

BREAKPOINT	{
printf("Ih.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
	SOLVE states METHOD cnexp
	gIh = gIhbar*m
	ihcn = gIh*(Vm-ehcn)
}

DERIVATIVE states	{
	m' = alpham_Ih(A_t, y) * (1 - m) - betam_Ih(A_t, y) * m
}

INITIAL{
printf("Ih.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
printf("Ih.mod: \n")
printf("V = %g\t",V(A_t,y))
printf("alpha = %g\t" ,alpham_Ih(A_t, y))
printf("beta = %g\t" ,betam_Ih(A_t, y))
	m = alpham_Ih(A_t, y) / (alpham_Ih(A_t, y) + betam_Ih(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}