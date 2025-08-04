COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
ENDCOMMENT

NEURON {
	POINT_PROCESS optodrive
	RANGE del, dur, gsin, Ese, f
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
    PI =   (pi) (1)
}

PARAMETER {
	del=0    (ms)
	dur    (ms)	<0,1e9>
	gsin   (uS)
	Ese=0  (mV)
    f=8    (1/s)
	v 	   (mV)
}
ASSIGNED { i (nA) }

INITIAL {
	i = 0
}

BREAKPOINT {
	at_time(del)
	at_time(del+dur)
	if (t < del + dur && t >= del) {
		if (f == 0) {
			i = gsin*(v-Ese)
		} else {
			i = gsin/2*(1+sin(2*PI*f*t*1e-3-PI/2))*(v-Ese)
			}
		}
	else{
			i = 0	
		}
}