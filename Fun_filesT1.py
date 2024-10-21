#!/usr/bin/env python

import scipy
from numpy import *

p_Alb_dia = [9.46*10**9.0, 6.03*10**(-4.0), 0.832, 65.7]
p_Hb_dia = [4.22*10**15.0, 4.55*10**(-4.0), 1.97, 42.9]
p_saline = [3.46*10**9.0, 7.45*10**(-11.0), 1.22, 1.90*10**(-3.0)]
p_MetHb_para = [6.40*10**6.0, 9.72*10**(-6.0), 3.07*10**(-10.0)];
p_Hb_para = [8.55*10**3.0, 8.68*10**(-9.0), 8.97*10**(-12.0)];


def Fun_koenig(p, B0):
	return p[0]*p[1]/(1.0+(267.513*10**6.0*B0*p[1])**p[2])+p[3]

def jw(w, tau, tau_s):
	comb1 = (1j*w*tau + tau/tau_s)**0.5;
	y1 = (1+0.25*comb1)/(1+comb1+4.0/9*comb1**2.0+1.0/9*comb1**3);
	return y1.real
	
def Fun_Freed(p, B0):
	w = B0*267.513*10**6.0;
	we = B0*1.760860*10**11.0; 
	return p[0]*(3*jw(w, p[1], p[2])+7*jw(we, p[1], p[2]));	


def Main_T1_cal(B0T, Hct, Y):
    ctHb = 33.2 # g/dL
    b = ctHb/((0.377*ctHb+100)-ctHb)/64.5 # mol/kg Equation 10, 100 is used rather than 1000 in equation 10 becuase g/dL is used rather than g/L
    b_alb = (0.75/1000)/(1.025-0.75/1000*66.5)  # Different than Equation 8, b_alb is calculated directly using n_protein/m_water, Alb concentration is 0.75mM, Denominator is the water mass. The results for two equation is 0.0007691 versus 0.0007695 
    MetHb=0.004 # 0.4%
    Ye = (1-MetHb)*Y

    if B0T == 3.0:
        r_MetHb = 1260
        r_dia = 62.9
        r_para = 65.6
        r_R1_water = 0.287
        r_alb = 172
    elif B0T == 7.0:
        r_MetHb = 1080
        r_dia = 45.1
        r_para = 56.2
        r_R1_water = 0.263
        r_alb = 114
    elif B0T == 9.4:
        r_para = 64.9
        r_R1_water = 0.232
        r_MetHb = 854.3
        r_dia = 47.2
        r_alb = 55.8
    elif B0T == 11.7: 
        r_para = 58.2
        r_R1_water = 0.235
        r_MetHb = 671.9
        r_dia = 43.3
        r_alb = 39.0
    elif B0T == 1.5:
        r_MetHb = 1340
        r_dia = 121
        r_para = 80.35
        r_R1_water = 0.253
        r_alb = 253
    elif B0T == 4.7: 
        r_MetHb = 1190
        r_dia = 51.1
        r_para = 61.68
        r_R1_water = 0.23
        r_alb = 139
    else:
        r_MetHb = Fun_Freed(p_MetHb_para, B0T)
        r_dia = Fun_koenig(p_Hb_dia, B0T)
        r_para = Fun_Freed(p_Hb_para, B0T)
        r_R1_water = Fun_koenig(p_saline, B0T)
        r_alb = Fun_koenig(p_Alb_dia, B0T)

    frac = 0.7*Hct/(0.7*Hct+0.95*(1-Hct)) #water fraction corrected with protein  Equation 2
    return 1000/ ( ((r_dia+r_para*(1-Ye-MetHb)+r_MetHb*MetHb)*b+r_R1_water)*frac+(r_R1_water+ b_alb*r_alb)*(1-frac))  #Equation 14
                  #  according to the definition above Ye = (1-MetHb)*Y, so 1-Ye - MetHb = 1 - (1-MetHb)*Y - MetHb = (1-Y)*(1-MetHb) same as Equation 14