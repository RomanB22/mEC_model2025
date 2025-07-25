TITLE sinusoidal input with v

NEURON {
	POINT_PROCESS optodrive
	RANGE  gsin, Ese, f
	NONSPECIFIC_CURRENT i
}

PARAMETER {

	gsin	(uS)
	Ese=0	(mV)
        f=8 (1/s)
	v 		(mV)
}


UNITS {

	(nA) = (nanoamp)
	(mV) = (millivolt)
        PI = (pi) (1)
} 

ASSIGNED {
	i 		(nA)


}
 

BREAKPOINT {
        
        i = gsin/2*(1+sin(2*PI*f*t*1e-3-PI/2))*(v-Ese) 
                    
} 

INITIAL {
	i=0

}

