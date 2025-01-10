: SK-type calcium-activated potassium current
: Reference : Kohler et al. 1996

NEURON {
       SUFFIX SK_E22
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gSK_E2bar, gSK_E2, ik
	RANGE Adrive, Vm, y, Fdrive, A_t, q1, f1, q2, f2 : section (even segment) specific
	RANGE stimon, detailed    : common to all sections (but set as RANGE to be accessible from caller)
}

UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      (mM) = (milli/liter)
	PI = (pi) (1) :in order to use the constant pi
}

PARAMETER {
	stimon       : Stimulation state
	Fdrive (kHz) : Stimulation frequency
	Adrive (kPa) : Stimulation amplitude
	detailed     : Simulation type
	v (nC/cm2)
	Vm (mV)
          gSK_E2bar = .000001 (mho/cm2)
          zTau = 1              (ms)
          ek           (mV)
          cai          (mM)
}

ASSIGNED {
         zInf
         ik            (mA/cm2)
         gSK_E2	       (S/cm2)
	A_t  (kPa)
	y
	q1  (nC/cm2)
	f1  (rad)
	q2  (nC/cm2)
	f2  (rad)
}

INCLUDE "update.inc"
FUNCTION_TABLE V(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (mV)
FUNCTION_TABLE A_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (mV)
FUNCTION_TABLE A_V2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (mV)
FUNCTION_TABLE phi_V1(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (rad)
FUNCTION_TABLE phi_V2(A(kPa), Q(nC/cm2), Q1(nC/cm2), phi1(rad), Q2(nC/cm2), phi2(rad)) (rad)

STATE {
      z   FROM 0 TO 1
}

BREAKPOINT {
	update()
           SOLVE states METHOD cnexp
           gSK_E2  = gSK_E2bar * z
           ik   =  gSK_E2 * (Vm - ek)
}

DERIVATIVE states {
        rates(cai)
        z' = (zInf - z) / zTau
}

PROCEDURE rates(ca(mM)) {
          if(ca < 1e-7){
	              ca = ca + 1e-07
          }
          zInf = 1/(1 + (0.00043 / ca)^4.8)
}

INITIAL {
        rates(cai)
        z = zInf
}
INDEPENDENT {
	t FROM 0 TO 1 WITH 1 (ms)
}