
:  This is a NEURON mod file generated from a ChannelML file

:  Unit system of original ChannelML file: Physiological Units

COMMENT
    ChannelML file containing a single Channel description
ENDCOMMENT

TITLE Channel: NaT

COMMENT
   Provided by Nolan based on Pastoll et al 2012
ENDCOMMENT


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (um) = (micrometer)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (l) = (liter)
}


    
NEURON {

    SUFFIX NaT
    USEION na READ ena WRITE ina VALENCE 1  : reversal potential of ion is read, outgoing current is written
           
        
    RANGE gmax, gion, mvhalf, mslope, hvhalf, hslope, mfactor, hfactor
    
    RANGE minf, mtau
    
    RANGE hinf, htau
    
}

PARAMETER { 

    gmax = 0.025 (S/cm2)  : default value, should be overwritten when conductance placed on cell
	
	mvhalf = -26.138
	mslope = 9.3829
	
	hvhalf = -23.774
	hslope = 6.108
    
	mfactor = 1
	hfactor = 1
}



ASSIGNED {

    v (mV)
    
    celsius (degC)
    
    : Reversal potential of na
    ena (mV)
    : The outward flow of ion: na calculated by rate equations...
    ina (mA/cm2)
    
    
    gion (S/cm2)
    minf
    mtau (ms)
    hinf
    htau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
        
    gion = gmax * (m^3) * (h^1)
    ina = gion*(v - ena)
            

}



INITIAL {
    
        
    rates(v)
    m = minf
	h = hinf
    gion = gmax * (m^3) * (h^1)
    ina = gion*(v - ena)       
    
}
    
STATE {
    m
    h
    
}



DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
	h' = (hinf - h)/htau
            

}

PROCEDURE rates(v(mV)) {  
    
    : Note: not all of these may be used, depending on the form of rate equations
    LOCAL  alpha, beta, tau, inf, gamma, zeta, temp_adj_m,
         A_alpha_m, B_alpha_m, Vhalf_alpha_m,
         A_beta_m, B_beta_m, Vhalf_beta_m, temp_adj_h,
         A_alpha_h, B_alpha_h, Vhalf_alpha_h,
         A_beta_h, B_beta_h, Vhalf_beta_h
    
    TABLE minf, mtau,hinf, htau DEPEND celsius FROM -100 TO 100 WITH 400
    
    UNITSOFF
    temp_adj_m = 1
    temp_adj_h = 1
    
        
    :      ***  Adding rate equations for gate: m  ***
        
    : Found a parameterised form of rate equation for alpha, using expression: A*((v-Vhalf)/B) / (1 - exp(-((v-Vhalf)/B)))
    A_alpha_m = 4
    B_alpha_m = 9
    Vhalf_alpha_m = -33 
    alpha = A_alpha_m * vtrap((v - Vhalf_alpha_m), B_alpha_m)
    
    
    : Found a parameterised form of rate equation for beta, using expression: A*((v-Vhalf)/B) / (1 - exp(-((v-Vhalf)/B)))
    A_beta_m = 27.6
    B_beta_m = -12
    Vhalf_beta_m = -58 
    beta = A_beta_m * vtrap((v - Vhalf_beta_m), B_beta_m)
    
    mtau = mfactor*(1/(temp_adj_m*(alpha + beta)))
	
    minf = 1/(1+exp((mvhalf-v)/mslope))
    


    :     *** Finished rate equations for gate: m ***
    

    
        
    :      ***  Adding rate equations for gate: h  ***
        
    : Found a parameterised form of rate equation for alpha, using expression: A*((v-Vhalf)/B) / (1 - exp(-((v-Vhalf)/B)))
    A_alpha_h = 0.36
    B_alpha_h = -12
    Vhalf_alpha_h = -48 
    alpha = A_alpha_h * vtrap((v - Vhalf_alpha_h), B_alpha_h)
    
    
    : Found a parameterised form of rate equation for beta, using expression: A*((v-Vhalf)/B) / (1 - exp(-((v-Vhalf)/B)))
    A_beta_h = 0.40
    B_beta_h = 6
    Vhalf_beta_h = -11 
    beta = A_beta_h * vtrap((v - Vhalf_beta_h), B_beta_h)
    
    htau = hfactor*(1/(temp_adj_h*(alpha + beta)))
    hinf = 1+(-1/(1+exp((hvhalf-v)/hslope)))
    


    :     *** Finished rate equations for gate: h ***
    

    
}


: Function to assist with parameterised expressions of type linoid/exp_linear

FUNCTION vtrap(VminV0, B) {
    if (fabs(VminV0/B) < 1e-6) {
    vtrap = (1 + VminV0/B/2)
}else{
    vtrap = (VminV0 / B) /(1 - exp((-1 *VminV0)/B))
    }
}

UNITSON


