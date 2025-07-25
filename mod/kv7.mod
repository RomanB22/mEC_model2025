TITLE Kv1

NEURON {
	SUFFIX kv7
	USEION k WRITE ik
	RANGE  gbar, tha1, siga1, ka1, ka2
}

PARAMETER {
	gbar = 0.010   	(mho/cm2)	:this should be read in
	tha1= 1 (mV) :this should be read in
	siga1=12 (mV)
	siga2=-80 (mV)
	ka2=0.02
        ka1=1
	ek=-90		(mV)            : must be explicitly def. in hoc
	celsius

	v 		(mV)
}


UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
} 

ASSIGNED {
	ik 		(mA/cm2)
	thegk		(S/cm2)
	alphaa
	betaa

}
 

STATE { a }

BREAKPOINT {
        SOLVE states METHOD euler
        thegk = gbar*a*a*a*a
	ik = thegk * (v - ek)
} 

INITIAL {  
	trates(v)  
}

DERIVATIVE states {   
        trates(v)      
        a' = alphaa*(1-a) - betaa*a
}

PROCEDURE trates(vm) {  
        
	alphaa = -ka1*(v-tha1)/(exp(-(v-tha1)/siga1)-1)
	betaa = ka2/exp(-v/siga2)

}
