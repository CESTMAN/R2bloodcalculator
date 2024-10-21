#!/usr/bin/env python

# in the terminal run:   python test.py 3 10 0.4 1

import Fun_filesT1

import argparse
import sys

try:
	parser = argparse.ArgumentParser()
	parser.add_argument("B0", help="B0(T)", type=float) 
	parser.add_argument("Hct", type=float)
	parser.add_argument("Y", type=float)                
	args = parser.parse_args()
	
	T1 = Fun_filesT1.Main_T1_cal(args.B0, args.Hct, args.Y)
	print T1
except:
    e = sys.exc_info()[0]
    print e
    
    