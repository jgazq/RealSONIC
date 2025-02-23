:Comment : LVA ca channel. Note: mtau is an approximation from the plots
:Reference : :		Avery and Johnston 1996, tau from Randall 1997
:Comment: shifted by -10 mv to correct for junction potential
NEURON	{
	SUFFIX Ca_LVAst2
	USEION ca READ eca WRITE ica
	RANGE gCa_LVAstbar, gCa_LVAst, ica
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
	gCa_LVAstbar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v (nC/cm2)
	Vm (mV)
	eca	(mV)
	ica	(mA/cm2)
	gCa_LVAst	(S/cm2)
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

FUNCTION_TABLE V(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (mV)
FUNCTION_TABLE alpham_CaLVAst2(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betam_CaLVAst2(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE alphah_CaLVAst2(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE betah_CaLVAst2(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (/ms)
FUNCTION_TABLE A_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (nC/cm2)
FUNCTION_TABLE B_1(A(kPa), Q(nC/cm2), A1(nC/cm2), B1(nC/cm2)) (nC/cm2)

STATE	{
	m
	h
}

BREAKPOINT	{
	update()
	SOLVE states METHOD cnexp
	gCa_LVAst = gCa_LVAstbar*m*m*h
	ica = gCa_LVAst*(Vm-eca)
}

DERIVATIVE states	{
	m' = alpham_CaLVAst2(A_t, y, a1, b1) * (1 - m) - betam_CaLVAst2(A_t, y, a1, b1) * m
	h' = alphah_CaLVAst2(A_t, y, a1, b1) * (1 - h) - betah_CaLVAst2(A_t, y, a1, b1) * h
}

INITIAL{
	update()
	m = alpham_CaLVAst2(A_t, y, a1, b1) / (alpham_CaLVAst2(A_t, y, a1, b1) + betam_CaLVAst2(A_t, y, a1, b1))
	h = alphah_CaLVAst2(A_t, y, a1, b1) / (alphah_CaLVAst2(A_t, y, a1, b1) + betah_CaLVAst2(A_t, y, a1, b1))
}

INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}