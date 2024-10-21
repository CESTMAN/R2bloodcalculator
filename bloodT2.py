#!/usr/bin/env python

# in the terminal run:   python test.py 3 10 0.4 1

import Fun_files

import argparse
import sys

try:
	parser = argparse.ArgumentParser()
	parser.add_argument("B0", help="B0(T)", type=float) 
	parser.add_argument("tau_cpmg", help="echo spacing(ms)", type=float) 
	parser.add_argument("Hct", type=float)
	parser.add_argument("Y", type=float)                
	args = parser.parse_args()
	
	T2 = Fun_files.Main_T2_cal(args.B0, args.tau_cpmg, args.Hct, args.Y)
	print T2
except:
    e = sys.exc_info()[0]
    print e
    
    