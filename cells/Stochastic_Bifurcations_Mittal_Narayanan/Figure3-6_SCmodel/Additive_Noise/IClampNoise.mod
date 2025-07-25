NEURON {
 POINT_PROCESS IClampnoise
GLOBAL i
 RANGE del,dur,f0,f1,r,torn,std,bias,on
 ELECTRODE_CURRENT i
}

UNITS {
 (nA) = (nanoamp)
}

PARAMETER {
 del=50  (ms)
 dur=10000 (ms)
 torn=50 (ms)
 std=0.5 (nA)
 f0=0.2  (nA)
 f1=0.8  (nA)
 r =60
 bias = 0 (nA)
}

ASSIGNED {
 ival (nA)
 i (nA)
 amp (nA)
 noise (nA)
 on (1)
}

INITIAL {
 i = 0
 on = 0
 seed(10)
 net_send(del, 1)
}

PROCEDURE seed(x) {
 set_seed(x)
}

BEFORE BREAKPOINT {
 if (on) {
  noise = normrand(1,std*1(/nA))*1(nA)
  ival = noise + bias 
 } else {
  ival = 0
 }
}

BREAKPOINT {
 i = ival
}

NET_RECEIVE (w) {
 if (flag == 1) {
  if (on == 0) {
   : turn it on
   on = 1
   : prepare to turn it off
   net_send(dur, 1)
  } else {
   : turn it off
   on = 0
  }
 }
}

