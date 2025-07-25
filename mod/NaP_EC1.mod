TITLE Persistent sodium current
: Implemented in Rubin and Cleland (2006) J Neurophysiology
: Adapted from Fransen et al (2004) Hippocampus

INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX NaP
	USEION na READ ena WRITE ina
	RANGE gbar, ina, mvhalf , mslope, hvhalf, hslope, mfactor, hfactor
}

PARAMETER { 
	gbar = 0.0 	(mho/cm2)
	v ena 		(mV)
	
	mvhalf = 48.7
	mslope = 4.4
	
	hvhalf = 48.8
	hslope = 9.98
	
	mfactor = 1
	hfactor = 1
} 

ASSIGNED { 
	ina 		(mA/cm2) 
	minf 		(1)
	mtau 		(ms)
	hinf
	htau		(ms)
} 
STATE {
	m h
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ina = gbar * m * h * ( v - ena )
} 

INITIAL { 
	settables(v) 
	m = minf
	h = hinf
	ina = gbar * m * h * ( v - ena )

} 

DERIVATIVE states { 
	settables(v) 
	m' = ( minf - m ) / mtau
	h' = ( hinf - h ) / htau 
}
UNITSOFF
 
PROCEDURE settables(v) { 
	TABLE minf, mtau, hinf, htau FROM -120 TO 40 WITH 641
	minf  = 1 / ( 1 + exp( -(v + mvhalf) / mslope ) ) 
	if( v == -38.0 ) {
		mtau = .0013071895424837 :limit as v --> -38, a discontinuity in the mtau function
	}else{
		mtau = mfactor*(1 / ((.091 * 1000 * (v + 38))/(1 - exp(-(v + 38)/5)) + (-.062 * 1000 * (v + 38))/(1 - exp((v + 38)/5))))
	}
	hinf = 1 / ( 1 + exp((v + hvhalf) / hslope ))
	htau = hfactor*(1 / ((-2.88*.001*(v + 17.049))/(1 - exp((v - 49.1)/4.63)) + (6.94*.001*(v + 64.409))/(1 - exp(-(v + 447)/2.63))))
}
UNITSON