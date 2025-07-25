TITLE Guillem's active synapses

NEURON {
	POINT_PROCESS synact
	RANGE  tau_fall, tau_rise, Es, f
	NONSPECIFIC_CURRENT Isyn
}

PARAMETER {

	tau_fall=2.3 (ms)
	tau_rise=0.4 (ms)
	Es=-75 (mV)
	f=1

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

}
 

STATE { g_fall g_rise }

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
}


NET_RECEIVE(weight (uS)) {
	g_fall = g_fall+weight
	g_rise = g_rise+weight
}
