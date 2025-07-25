TITLE Low voltage activated calcium channel 
: it calculates I_Ca using channel permeability instead of conductance
: Based on Bruehl and Wadman 1999 and Pastoll et al 2012.
: alphas and betas fitted .
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degK)
	KTOMV = .0853 (mV/degC)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {           :parameters that can be entered when function is called in cell-setup 
	dt            (ms)
	v             (mV)
    tBase = 23.5  (degC)
	celsius (degC) := 22  
	glvabar = 0   (mho/cm2)  : initialized conductance
	ki = 0.001    (mM)
	cai = 10.0e-5 (mM)       : initial internal Ca++ concentration
	cao = 2       (mM)       : initial external Ca++ concentration
    tfa = 1                  : activation time constant scaling factor
    tfi = 1.2                : inactivation time constant scaling factor
    eca = 140                : Ca++ reversal potential
	
	mvhalf = -52.431 (mV)
	mslope = 8.2023  (mV)
	
	hvhalf = -88.165 (mV)
	hslope = 6.6772	 (mV)
	
	mfactor = 1
	hfactor = 1
}

NEURON {
	SUFFIX LVA
	USEION ca READ cai,cao WRITE ica
	RANGE glvabar, hinf, minf, taum, tauh, mvhalf, mslope, hvhalf, hslope, mfactor, hfactor
}

STATE {	m h }  : unknown activation and inactivation parameters to be solved in the DEs 

ASSIGNED {     : parameters needed to solve DE
	ica (mA/cm2)
	glva  (mho/cm2) 
	minf
	hinf
	taum
	tauh
}

INITIAL {
:   tadj = 3^((celsius-tBase)/10)   : assume Q10 of 3
	rates(v)
        m = minf
        h = hinf
	glva = glvabar*m*m*h*h2(cai)
	ica = glva*ghk(v,cai,cao)    : dummy calcium current induced by this channel

}

BREAKPOINT {
	SOLVE states METHOD cnexp
	glva = glvabar*m*m*h*h2(cai) : maximum channel permeability
	ica = glva*ghk(v,cai,cao)    : dummy calcium current induced by this channel

}

UNITSOFF
FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) { LOCAL nu,f
        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {   : temperature-dependent adjustment factor
        KTF = ((25./293.15)*(celsius + 273.15))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alph(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alph = (1.6e-4)*exp(-(v+79.5)/20) : +89.5
}

FUNCTION beth(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
	beth = 1/(exp((-v-5)/10)+1.0)  : +5
}

FUNCTION alpm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	alpm = 0.8967*(-1.0*v-7.88)/(exp((-1.0*v-7.88)/10.0)-1.0) :+3.88 
}

FUNCTION betm(v(mV)) {
	TABLE FROM -150 TO 150 WITH 200
	betm = 0.046*exp(-v/22.73)
}

UNITSON
LOCAL facm,fach

:if state_cagk is called from hoc, garbage or segmentation violation will
:result because range variables won't have correct pointer.  This is because
: only BREAKPOINT sets up the correct pointers to range variables.

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        m' = (minf - m)/taum
        h' = (hinf - h)/tauh
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a
        a = alpm(v)
        taum = mfactor*(1/(tfa*(a + betm(v)))) : estimation of activation tau
        minf =  1/(1+exp((mvhalf-v)/mslope))        : estimation of activation steady state
:        facm = (1 - exp(-dt/taum))
        a = alph(v)
        tauh = hfactor*(1/(tfi*(a + beth(v)))) : estimation of inactivation tau
        hinf = 1+(-1/(1+exp((hvhalf-v)/hslope)))         : estimation of inactivation steady state
 :       fach = (1 - exp(-dt/tauh))
}
