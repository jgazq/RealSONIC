:Reference : :		Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones
NEURON	{
	SUFFIX Im2
	USEION k READ ek WRITE ik
	RANGE gIm2bar, gIm2, ik
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
	gIm2bar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ek	(mV)
	ik	(mA/cm2)
	gIm2	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
	A_t  (kPa)
	y
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_Im2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_Im2(A(kPa), Q(nC/cm2)) (/ms)

STATE	{ 
	m
}

BREAKPOINT	{
printf("Im2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
	SOLVE states METHOD cnexp
	gIm2 = gIm2bar*m
	ik = gIm2*(Vm-ek)
}

DERIVATIVE states	{
	m' = alpham_Im2(A_t, y) * (1 - m) - betam_Im2(A_t, y) * m
}

INITIAL{
printf("Im2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
printf("Im2.mod: \n")
printf("V = %g\t",V(A_t,y))
printf("alpha = %g\t" ,alpham_Im2(A_t, y))
printf("beta = %g\t" ,betam_Im2(A_t, y))
	m = alpham_Im2(A_t, y) / (alpham_Im2(A_t, y) + betam_Im2(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}