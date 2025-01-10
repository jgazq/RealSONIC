:Reference :Colbert and Pan 2002
NEURON	{
	SUFFIX NaTa_t
	USEION na READ ena WRITE ina
	RANGE gNaTa_tbar, gNaTa_t, ina
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
	gNaTa_tbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTa_t	(S/cm2)
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
FUNCTION_TABLE alpham_NaTat(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_NaTat(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_NaTat(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betah_NaTat(A(kPa), Q(nC/cm2)) (/ms)

STATE	{
	m
	h
}

BREAKPOINT	{
printf("NaTa_t.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
	SOLVE states METHOD cnexp
	gNaTa_t = gNaTa_tbar*m*m*m*h
	ina = gNaTa_t*(Vm-ena)
}

DERIVATIVE states	{
	m' = alpham_NaTat(A_t, y) * (1 - m) - betam_NaTat(A_t, y) * m
	h' = alphah_NaTat(A_t, y) * (1 - h) - betah_NaTat(A_t, y) * h
}

INITIAL{
printf("NaTa_t.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
printf("NaTa_t.mod: \n")
printf("V = %g, alpha = %g, beta = %g\n",V(A_t,y), alpham_NaTat(A_t, y), betam_NaTat(A_t, y))
	m = alpham_NaTat(A_t, y) / (alpham_NaTat(A_t, y) + betam_NaTat(A_t, y))
printf("NaTa_t.mod: \n")
printf("V = %g, alpha = %g, beta = %g\n",V(A_t,y), alphah_NaTat(A_t, y), betah_NaTat(A_t, y))
	h = alphah_NaTat(A_t, y) / (alphah_NaTat(A_t, y) + betah_NaTat(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}