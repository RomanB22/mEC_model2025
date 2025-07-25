TITLE Guillem's active synapses with depression

NEURON {
	POINT_PROCESS synactdep
	RANGE  tau_fall, tau_rise, Es, f, U_SE, tau_d
	    NONSPECIFIC_CURRENT Isyn
}

PARAMETER {

	tau_fall=2.3 (ms)
	tau_rise=0.4 (ms)
	tau_d=1 (ms) :will re read in
	Es=-75 (mV)
	f=1
	U_SE=1

	v 		(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	Isyn 		(uS)
	c_rise
	c_fall
	rs

}
 

STATE { g_fall g_rise xs }

BREAKPOINT {
        SOLVE states METHOD cnexp
        Isyn = (v-Es)*f*(g_fall-g_rise) : amp 
                    
} 

INITIAL {
	g_fall=0 
	g_rise=0 

}

DERIVATIVE states {        
		g_fall' = -g_fall/tau_fall : siemens
        g_rise' = -g_rise/tau_rise : siemens
		xs'= (1 - xs)/tau_d
}


NET_RECEIVE(weight (uS)) {
	rs=U_SE*xs
	xs=xs-rs
	g_fall = g_fall+weight*rs
	g_rise = g_rise+weight*rs
}
