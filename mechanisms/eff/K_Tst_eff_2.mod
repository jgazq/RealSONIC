:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
NEURON	{
	SUFFIX K_Tst2
	USEION k READ ek WRITE ik
	RANGE gK_Tst2bar, gK_Tst2, ik
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
printf("K_Tst2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
	SOLVE states METHOD cnexp
	gK_Tst2 = gK_Tst2bar*(m^4)*h
	ik = gK_Tst2*(Vm-ek)
}

DERIVATIVE states	{
	m' = alpham_KTst2(A_t, y) * (1 - m) - betam_KTst2(A_t, y) * m
	h' = alphah_KTst2(A_t, y) * (1 - h) - betah_KTst2(A_t, y) * h
}

INITIAL{
printf("K_Tst2.mod: \n")
printf("V = %g\n",V(A_t,y))
	update()
printf("K_Tst2.mod: \n")
printf("V = %g\t",V(A_t,y))
printf("alpha = %g\t" ,alpham_KTst2(A_t, y))
printf("beta = %g\t" ,betam_KTst2(A_t, y))
	m = alpham_KTst2(A_t, y) / (alpham_KTst2(A_t, y) + betam_KTst2(A_t, y))
printf("K_Tst2.mod: \n")
printf("V = %g\t",V(A_t,y))
printf("alpha = %g\t" ,alphah_KTst2(A_t, y))
printf("beta = %g\t" ,betah_KTst2(A_t, y))
	h = alphah_KTst2(A_t, y) / (alphah_KTst2(A_t, y) + betah_KTst2(A_t, y))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}