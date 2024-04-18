:Reference :Colbert and Pan 2002
:comment: took the NaTa and shifted both activation/inactivation by 6 mv
NEURON	{
	SUFFIX NaTs2_t2
	USEION na READ ena WRITE ina
	RANGE gNaTs2_tbar, gNaTs2_t, ina
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
	gNaTs2_tbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	ena	(mV)
	ina	(mA/cm2)
	gNaTs2_t	(S/cm2)
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
FUNCTION_TABLE alpham_NaTs2t(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betam_NaTs2t(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_NaTs2t(A(kPa), Q(nC/cm2)) (/ms)
FUNCTION_TABLE betah_NaTs2t(A(kPa), Q(nC/cm2)) (/ms)

STATE	{
	m
	h
}

BREAKPOINT	{
printf("NaTs2_t2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
	SOLVE states METHOD cnexp
	gNaTs2_t = gNaTs2_tbar*m*m*m*h
	ina = gNaTs2_t*(Vm-ena)
}

DERIVATIVE states	{
	m' = alpham_NaTs2t(A_t, y) * (1 - m) - betam_NaTs2t(A_t, y) * m
	h' = alphah_NaTs2t(A_t, y) * (1 - h) - betah_NaTs2t(A_t, y) * h
}

INITIAL{
printf("NaTs2_t2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
printf("NaTs2_t2.mod: \n")
printf("V = %g, alpha = %g, beta = %g\n",V(A_t,y), alpham_NaTs2t(A_t, y), betam_NaTs2t(A_t, y))
	m = alpham_NaTs2t(A_t, y) / (alpham_NaTs2t(A_t, y) + betam_NaTs2t(A_t, y))
printf("NaTs2_t2.mod: \n")
printf("V = %g, alpha = %g, beta = %g\n",V(A_t,y), alphah_NaTs2t(A_t, y), betah_NaTs2t(A_t, y))
	h = alphah_NaTs2t(A_t, y) / (alphah_NaTs2t(A_t, y) + betah_NaTs2t(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}