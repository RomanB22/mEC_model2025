: Sah P, Clements JD (1999) Photolytic manipulation of [Ca2+]i
: reveals slow kinetics of potassium channels underlying the
: afterhyperpolarization in hippocampal pyramidal neurons.
: J  Neurosci 19:3657-3664.

: Sah P, Isaacson JS (1995) Channels underlying the slow
: afterhyperpolarization in hippocampal pyramidal neurons:
: neurotransmitters modulate the open probability. Neuron 15:435-441.
NEURON {                                                                                                                        
    SUFFIX skkin                                                                                                            
    USEION k READ ek WRITE ik                                                                                           
    USEION ca READ cai                                                                                                          
    RANGE g, gbar                                                                                                   
}                                                                                                                                                                                                                                               
UNITS {                                                                                                                         
    (mA) = (milliamp)                                                                                                       
    (mV) = (millivolt)                                                                                                      
    (molar) = (1/liter)                                                                                                     
    (uM)    = (micromolar)                                                                                                  
    (mM)    = (millimolar)                                                                                          
}                                                                                                                                                                                                                                               
PARAMETER {                                                                                                                    
    cai0 = 50e-6  (mM)                                                                                                      
    gbar = 33     (mho/cm2)                                                                                                 
    kf   = 10e-6  (mM/ms) : 10 microM/s = 10e-6 mM/ms                                                                       
    kb   = 0.5e-3 (/ms)   : 0.5 /s = 0.5e-3 /ms                                                                             
    kfo  = 600e-3 (/ms)   : 600 /s = 600e-3 /ms                                                                             
    kbo  = 400e-3 (/ms)   : 400 /s = 400e-3 /ms                                                                     
}                                                                                                                                                                                                                                               
ASSIGNED {                                                                                                                      
    v    (mV)                                                                                                               
    ek   (mV)                                                                                                               
    g    (mho/cm2)                                                                                                          
    ik   (amp/cm2)                                                                                                          
    cai  (mM)                                                                                                       
}                                                                                                                                                                                                                                               
STATE {                                                                                                                         
    c cac ca2c ca3c ca4c o calcm                                                                                    
}                                                                                                                                                                                                                                               
INITIAL {                                                                                                                       
    calcm=cai0 * 1e6                                                                                                        
    SOLVE kin STEADYSTATE sparse                                                                                            
    g = gbar*o                                                                                                              
    ik = g*(v - ek)                                                                                                 
}                                                                                                                                                                                                                                               
BREAKPOINT {                                                                                                                    
    calcm=cai * 1e6                                                                                                         
    SOLVE kin METHOD sparse                                                                                                 
    g = gbar*o                                                                                                              
    ik = g*(v - ek)                                                                                                 
}                                                                                                                                                                                                                                                                                                                                                                       
KINETIC kin {                                                                                                                   
    ~ c + calcm     <-> cac  (4*kf, kb)                                                                                     
    ~ cac + calcm   <-> ca2c (3*kf, 2*kb)                                                                                   
    ~ ca2c + calcm  <-> ca3c (2*kf, 3*kb)                                                                                   
    ~ ca3c + calcm  <-> ca4c (kf, 4*kb)                                                                                     
    ~ ca4c          <-> o    (kfo, kbo)                                                                                     
    CONSERVE c + cac + ca2c + ca3c + ca4c + o = 1                                                                   
}     
