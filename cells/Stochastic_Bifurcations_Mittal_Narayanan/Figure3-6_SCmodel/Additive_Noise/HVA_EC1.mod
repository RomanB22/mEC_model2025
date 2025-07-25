TITLE Ca-High voltage activated channel for stellate cells in layer II of Medial entorhinal cortex


COMMENT 
  Kinetics taken from Castelli and Magistretti 2006, Bruehl and Wadman 1999 and Pastoll et al 2012
  Tiime constants are modified from R-type Ca channel from CA1 based on fitting of data from above papers

ENDCOMMENT

NEURON {
    SUFFIX HVA
    USEION ca READ cai, cao WRITE ica
    RANGE gmax, m, h, mvhalf, hvhalf, mslope, hslope, mfactor, hfactor
    RANGE minf, hinf, taum, tauh
    GLOBAL q10, taum_fit, tauh_fit, z
}

UNITS {
    (molar) = (/liter)
    (mA) = (milliamp)
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (mM) = (millimolar)
    (S) = (siemens)
    (uS) = (microsiemens)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {   
    gmax = 0      (S/cm2) <0,1e9> 
    q10  = 3  
    taum_fit = 0.92  (ms)
	tauh_fit = 250   (ms)       : fitting pof taum and tauh
    z = 2                       : valency of Ca ions
	
	mvhalf = 11.1	(mV)
	hvhalf = 37		(mV)
	mslope = 8.4	(mV)
	hslope = 9		(mV)
	mfactor = 1		(mV)
	hfactor = 1		(mV)
}  

STATE {	mO mC hO hC }    

ASSIGNED {               : parameters needed to solve DE
    v       (mV)
    celsius (degC)
    cai     (mM)
    cao     (mM)
	ica     (mA/cm2)
    minf
    hinf
	taum    (ms)
    tauh    (ms)
}

BREAKPOINT {
    SOLVE kin METHOD sparse
	ica = gmax*mO*mO*mO*hO*ghkg(v,cai,cao,z)
}

INITIAL { 
    :taum = q10^(-(celsius-22(degC))/10(degC))*taum_fit
	taum = mfactor*taum_fit
    :tauh = q10^(-(celsius-22(degC))/10(degC))*tauh_fit
	tauh = hfactor*tauh_fit
    SOLVE kin STEADYSTATE sparse    
    ica = gmax*mO*mO*mO*hO*ghkg(v,cai,cao,z)
}

KINETIC kin {
    minf = 1/(1+exp(-(v+ mvhalf)/mslope))
    hinf = 1/(1+exp( (v+hvhalf)/hslope))
    ~ mC <-> mO (minf/taum, (1-minf)/taum)
    ~ hC <-> hO (hinf/tauh, (1-hinf)/tauh)
    CONSERVE mC + mO = 1
    CONSERVE hC + hO = 1
}

FUNCTION ghkg(v(mV), ci(mM), co(mM), z) (mV) {
    LOCAL xi, f, exi, fxi
    f = R*(celsius+273.15)/(z*(1e-3)*FARADAY)
    xi = v/f
    exi = exp(xi)
    if (fabs(xi) < 1e-4) {
        fxi = 1 - xi/2
    }else{
        fxi = xi/(exi - 1)
    }
    ghkg = f*((ci/co)*exi - 1)*fxi
}

