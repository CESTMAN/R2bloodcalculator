#!/usr/bin/env python

import scipy
from scipy.integrate import quad
from numpy import *
from scipy.optimize import least_squares

import Fun_filesT1

def Fun_Rex_G_lysed(dw1, B0, R21, R22, frac, tau_cpmg, tau):
	dw = 2.0 * math.pi * B0 * dw1
	kex = 1/tau/(1-frac) #unit per second
	p1 = (R21-R22-frac*kex+(1-frac)*kex)**2.0 - dw**2.0 + 4.0*frac*(1-frac)*kex**2
	p2 = 2.0*dw*(R21-R22-frac*kex+(1-frac)*kex)
	pp = tau_cpmg/(8**0.5)*(p1+(p1**2+p2**2)**0.5)**0.5
	pm = tau_cpmg/(8**0.5)*(-1.0*p1+(p1**2+p2**2)**0.5)**0.5
	dp = 0.5*(1.0+(p1+2*dw**2)/(p1**2+p2**2)**0.5)
	dm = 0.5*(-1.0+(p1+2*dw**2)/(p1**2+p2**2)**0.5)
	kk = cosh(2*pp)*dp-dm*cos(2*pm)
	Rex_general = 0.5*(kex+R21+R22-math.acosh(kk)/tau_cpmg)  
	return Rex_general


def Fun_Rex_G(dw1, B0, R21, R22, Hct, tau_cpmg, tau):
	frac = Hct*0.7/(Hct*0.7+(1-Hct)*0.95)
	dw = 2.0 * math.pi * B0 * dw1 # unit rad/second
	kex = 1./tau/(1-frac)  #unit per second
	p1 = (R21-R22-frac*kex+(1-frac)*kex)**2.0 - dw**2.0 + 4.0*frac*(1-frac)*kex**2
	p2 = 2.0*dw*(R21-R22-frac*kex+(1-frac)*kex)
	pp = tau_cpmg/(2**0.5)*(p1+(p1**2+p2**2)**0.5)**0.5
	pm = tau_cpmg/(2**0.5)*(-1.0*p1+(p1**2+p2**2)**0.5)**0.5
	dp = 0.5*(1.0+(p1+2*dw**2)/(p1**2+p2**2)**0.5)
	dm = 0.5*(-1.0+(p1+2*dw**2)/(p1**2+p2**2)**0.5)
	kk = cosh(pp)*dp-dm*cos(pm)
	Rex_general = 0.5*(kex+R21+R22-math.acosh(kk)/tau_cpmg)
	return Rex_general
	
	
def Fun_Rex_Jensen(G0, B0, tau_diffusion, tau_cpmg, Hct):
	x0 = 4*tau_cpmg/tau_diffusion

	def fun(y, x):
		return (1/pi**0.5)*(exp(-y)/y**0.5*(1.0-tanh(x0*y)/x0/y))
	Integral1 = quad(fun, 0, inf, args=(x0))
	
	Rex_diffusion = Hct*G0**2*B0**2*tau_diffusion*Integral1[0]
	return Rex_diffusion	


def Fun_R2_lysed_predict(B0t, Y, cHb, tau_cpmg):
	if tau_cpmg > 40.0/1000:
	    tau_cpmg = 40./1000

	B0 = B0t * 42.577 
	R2b_o = 94.7*B0t #averaged R2 of the exchangeable protons in oxygenated Hb
	R2b_p = 10.6*B0t**2 + 287  #averaged R2 of the exchangeable protons in deoxygenated Hb

	dw_p = 1.11 #ppm average chemical shift between deoxyHb exchange H and water H
	dw_o = 0.614 #ppm average chemical shift between oxyHb exchange H and water H
	R2w = 0.36  # PBS solution R2
	n_Hb = 499  # number of exchangable proton in one Hb
	frac = n_Hb*cHb/64.5/1000/(2*55.6*(1-0.3*cHb/332))
	tau1 = 1.0/(14.6*1000)

	Rex_p = Fun_Rex_G_lysed(dw_p, B0, R2b_p, R2w, frac, tau_cpmg, tau1)
	Rex_o = Fun_Rex_G_lysed(dw_o, B0, R2b_o, R2w, frac, tau_cpmg, tau1)

	Rex_result = Rex_p*(1-Y) + Rex_o*Y
	return Rex_result

	
def Fun_R2_plasma_predict(B0t, tau_cpmg):
	cAlb = 50 #g/dL
	B0 = B0t * 42.577
	R2b = 55.8*B0t  #averaged R2 of the exchangeable protons in Alb

	if tau_cpmg > 40.0/1000:
	    tau_cpmg = 40./1000

	dw = 1.21 #ppm average chemical shift between Alb exchange H and water H
	R2w = 0.36 # PBS solution R2
	n_Hb = 806 # number of exchangable proton in one Alb
	p1 = n_Hb*cAlb/66.5/1000/(2*55.6*(1-0.05*cAlb/50))
	tau1 = 1/(11.6*1000)
     
	R2 = Fun_Rex_G_lysed(dw, B0, R2b, R2w, p1, tau_cpmg, tau1)
	return R2	

	
def Main_T2_cal(B0, tau_cpmg, Hct, Y):
	# B0: T, tau_cpmg: ms, 
	tau_cpmg = tau_cpmg/1000 #s

	## calculate intrinsic R2 of plasma and cytoplasma 
	R2_p = Fun_R2_plasma_predict(B0, tau_cpmg)
	R2_e = Fun_R2_lysed_predict(B0, Y, 332, tau_cpmg)

	## R2 contribution from membrane
	R2_m = 2.39 # per second

	## exchange and diffusion contribution
	kai = 0.253 #ppm, theoritically calculated susceptibility difference when fully deoxygeanted
	tau_ery = 9.13/1000 #s
	beta_ex = -1.20  # exchange shape factor
	Yoff_ex = 0.889  

	tau_D2 = 3.15/1000 # s  diffusion coeffecient in the plasma
	beta_diff = 0.661 # diffusion shape factor
	Yoff_diff = 0.985

	# chemical shift difference beween inside and outside the erythrocyte
	dw = kai*beta_ex*(Yoff_ex-Y)+0.0298-0.0164*Y  
	# susceptibility difference beween inside and outside the erythrocyte
	kai_d = kai*(Yoff_diff-Y)
	# pre-factor in the Jensen model for diffusion
	g2 = math.pi**2*(128.0/45.0)**0.5 * kai_d * beta_diff 

	# diffusion contribution
	r_d1 = 0 #erythrocyte
	r_d2 = Fun_Rex_Jensen(g2, B0* 42.577, tau_D2, tau_cpmg, Hct) #plasma

	result1 = 1000/Fun_Rex_G(dw, B0* 42.577, (R2_e+r_d1 + R2_m), (r_d2+R2_p + R2_m*Hct), Hct, tau_cpmg, tau_ery)
            
	return result1

def Main_T2corr_cal(B0, tau_cpmg, Hct, Y, tau_p):	    
	T1 = Fun_filesT1.Main_T1_cal(B0, Hct, Y)
	T2 = Main_T2_cal(B0, tau_cpmg, Hct, Y)
	TE_corr = tau_cpmg - tau_p/2*(1-T2/T1)
	T2corr = Main_T2_cal(B0, TE_corr, Hct, Y)
	return T2corr

def Y_cal1(Y, B0, tau_cpmg, Hct, tau_p, T2corr):	
    T1 = Fun_filesT1.Main_T1_cal(B0, Hct, Y)
    T2 = Main_T2corr_cal(B0, tau_cpmg, Hct, Y, tau_p)
    return T2corr - T2

def Main_Y_cal(B0, tau_cpmg, Hct, tau_p, T2):
    Y1 = least_squares(Y_cal1,0.6,bounds=([0,1]),args=(B0, tau_cpmg, Hct, tau_p, T2))
    
    return Y1.x[0]

def Hct_cal1(Hct, B0, Y, T1):	
    T1_temp = Fun_filesT1.Main_T1_cal(B0, Hct, Y)
    return T1 - T1_temp

def Main_Hct_cal(B0, Y, T1):	
    Hct1 = least_squares(Hct_cal1,0.42,bounds=([0.2,0.6]),args=(B0, Y, T1))
    
    return Hct1.x[0]
