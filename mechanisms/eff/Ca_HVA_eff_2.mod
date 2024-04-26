:Comment :
:Reference : :		Reuveni, Friedman, Amitai, and Gutnick, J.Neurosci. 1993

NEURON	{
	SUFFIX Ca_HVA2
	USEION ca READ eca WRITE ica
	RANGE gCa_HVA2bar, gCa_HVA2, ica 
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
	gCa_HVA2bar = 0.00001 (S/cm2) 
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa	(S/cm2)
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
}

INCLUDE "update.inc"

FUNCTION_TABLE V(A(kPa), Q(nC/cm2)) (mV)
FUNCTION_TABLE alpham_CaHVA2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_CaHVA2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_CaHVA2(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betah_CaHVA2(A(kPa), Q(nC/cm2)) (/ms)

STATE	{ 
	m
	h
}

BREAKPOINT	{
printf("Ca_HVA2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
	SOLVE states METHOD cnexp
	gCa = gCa_HVA2bar*m*m*h
	ica = gCa*(Vm-eca)
}

DERIVATIVE states	{
	m' = alpham_CaHVA2(A_t, y) * (1 - m) - betam_CaHVA2(A_t, y) * m
	h' = alphah_CaHVA2(A_t, y) * (1 - h) - betah_CaHVA2(A_t, y) * h
}

INITIAL{
printf("Ca_HVA2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
printf("Ca_HVA2.mod: \n")
printf("V = %g\t",V(A_t,y))
printf("alpha = %g\t" ,alpham_CaHVA2(A_t, y))
printf("beta = %g\t" ,betam_CaHVA2(A_t, y))
	m = alpham_CaHVA2(A_t, y) / (alpham_CaHVA2(A_t, y) + betam_CaHVA2(A_t, y))
printf("Ca_HVA2.mod: \n")
printf("V = %g\t",V(A_t,y))
printf("alpha = %g\t" ,alphah_CaHVA2(A_t, y))
printf("beta = %g\t" ,betah_CaHVA2(A_t, y))
	h = alphah_CaHVA2(A_t, y) / (alphah_CaHVA2(A_t, y) + betah_CaHVA2(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}