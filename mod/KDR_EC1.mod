
:  This is a NEURON mod file generated from a ChannelML file

:  Unit system of original ChannelML file: Physiological Units

COMMENT
    ChannelML file containing a single Channel description
ENDCOMMENT

TITLE Channel: KDR

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

    SUFFIX KDR
    USEION k READ ek WRITE ik VALENCE 1  : reversal potential of ion is read, outgoing current is written
           
        
    RANGE gmax, gion, nvhalf, nslope, nfactor
    
    RANGE ninf, ntau
    
}

PARAMETER { 

    gmax = 0.025 (S/cm2)  : default value, should be overwritten when conductance placed on cell
	
	nvhalf = -17.655
	nslope = 19.655
	nfactor = 1
	
}



ASSIGNED {

    v (mV)
    
    celsius (degC)
    
    : Reversal potential of k
    ek (mV)
    : The outward flow of ion: k calculated by rate equations...
    ik (mA/cm2)
    
    
    gion (S/cm2)
    ninf
    ntau (ms)
    
}

BREAKPOINT { 
                        
    SOLVE states METHOD cnexp
        
    gion = gmax * (n^4)
    ik = gion*(v - ek)
            

}



INITIAL {
        
    rates(v)
    n = ninf
    gion = gmax * (n^4)
    ik = gion*(v - ek)        
    
}
    
STATE {
    n
    
}



DERIVATIVE states {
    rates(v)
    n' = (ninf - n)/ntau
            

}

PROCEDURE rates(v(mV)) {  
    
    : Note: not all of these may be used, depending on the form of rate equations
    LOCAL  alpha, beta, tau, inf, gamma, zeta, temp_adj_n,
         A_alpha_n, B_alpha_n, Vhalf_alpha_n,
         A_beta_n, B_beta_n, Vhalf_beta_n
    
    TABLE ninf, ntau DEPEND celsius FROM -100 TO 100 WITH 400
    
    UNITSOFF
    temp_adj_n = 1
    
        
    :      ***  Adding rate equations for gate: n  ***
        
    : Found a parameterised form of rate equation for alpha, using expression: A*((v-Vhalf)/B) / (1 - exp(-((v-Vhalf)/B)))
    A_alpha_n = 0.2
    B_alpha_n = 10
    Vhalf_alpha_n = -38 
    alpha = A_alpha_n * vtrap((v - Vhalf_alpha_n), B_alpha_n)
    
    
    : Found a parameterised form of rate equation for beta, using expression: A*((v-Vhalf)/B) / (1 - exp(-((v-Vhalf)/B)))
    A_beta_n = 0.6294
    B_beta_n = -35
    Vhalf_beta_n = -47 
    beta = A_beta_n * vtrap((v - Vhalf_beta_n), B_beta_n)
    
    ntau = nfactor*(1/(temp_adj_n*(alpha + beta)))
	
    ninf = 1/(1+exp((nvhalf-v)/nslope))
    


    :     *** Finished rate equations for gate: n ***
    

    
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


