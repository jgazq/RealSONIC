: Dynamics that track inside calcium concentration
: modified from Destexhe et al. 1994

NEURON	{
	SUFFIX CaDynamics_E2
	USEION ca READ ica WRITE cai
	RANGE decay, gamma, minCai, depth
	RANGE Adrive, Vm, y, Fdrive, A_t : section specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
}

UNITS	{
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um)	= (micron)
}

PARAMETER	{
	stimon       : Stimulation state
	Fdrive (kHz) : Stimulation frequency
	Adrive (kPa) : Stimulation amplitude
	detailed     : Simulation type
	gamma = 0.05 : percent of free calcium (not buffered)
	decay = 80 (ms) : rate of removal of calcium
	depth = 0.1 (um) : depth of shell
	minCai = 1e-4 (mM)
}

ASSIGNED	{
	ica (mA/cm2)
	A_t  (kPa)
	y
}
STATE	{
	cai (mM)
	}

BREAKPOINT	{
	 SOLVE states METHOD cnexp 
}
DERIVATIVE states	{
	cai' = -(10000)*(ica*gamma/(2*FARADAY*depth)) - (cai - minCai)/decay
}
INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}