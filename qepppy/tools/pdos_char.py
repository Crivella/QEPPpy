#! python

import os.path
import re

def pdos_char( fname="", kpt_list=[], bnd_list=[], thr=0):
	if( not os.path.isfile(fname)):
		print( "ERROR: File fname:'{}' does not exist.".format( fname))
		return 1
	if( not kpt_list):
		print( "ERROR: Passing empty 'kpt_list'.")
		return 1
	if( isinstance( kpt_list, str)):
		kl = list( map( int, kpt_list.split(",")))
	if( isinstance( kpt_list, list)):
		kl = kpt_list
	if( isinstance( bnd_list, str)):
		bl = list( map( int, bnd_list.split(",")))
	if( isinstance( bnd_list, list)):
		bl = bnd_list
	if( len( bl) > 0):
		cb = True

	state_l = []

	k = -1
	b = -1
	ck = False
	ce = False
	cp = False
	cb = False
	if( len( bl) > 0):
		cb = True
	for line in open( fname):
		line = line.rstrip()
		if( " state #" in line):
			state_l.append( line.split(":")[1])
		if( " k = " in line or ck):
			if( " k = " in line):
				k+=1
				b=-1
				kpt = list( map( float, list( filter( None, line.split( " ")))[2:5]))
				ck = True
				ce = False
				if( k+1 not in kl):
					ck = False
				else:
					print( "KPT (#{:5d}):\t{}".format( k+1, kpt))
			if( ck and (" e(" in line or " e =" in line)):
				if( " e(" in line):
					off=4
				else:
					off=2
				ce = True
				b+=1
				
				#print( line, line.split( " "), list( filter( None, line.split( " "))))
				el = float( list( filter( None, line.split( " ")))[off])
				if( cb and b+1 not in bl):
					ce = False
				else:
					print("\tBND (#{:3d}):\t{} eV".format( b+1, el))
					
			if( cp and " |psi|^2" in line):
				cp = False
			if( ce and (" psi = " in line or cp)):
				cp = True
				l = list( filter( None, re.split( " +psi = | +\+|\*\[#|\]\+", line)))
				#print( l[1::2], l[0::2])
				for s, wf in zip( l[1::2], l[0::2]):
					wf=float(wf)
					if( wf>thr):
						print( "\t\t{}:\t{:7.3f}%".format( state_l[int(s)-1], wf*100))


	return 0



if __name__ == "__main__":
	import sys
	argc = len( sys.argv)
	if( not 2<=argc<=5):
		print("Incorrect use. Pleas pass arguments:"
			"\n\t'fname',"
			"\n\t'kpt_list\t(comma separated)',"
			"\n\t'bnd_list\t(optional) (comma separated)'"
			"\n\t'thr\t(optional) (threshold for printing percantages)'")
		exit()
	if( argc==2):
		pdos_char( sys.argv[1])
	if( argc==3):
		pdos_char( sys.argv[1], sys.argv[2])
	if( argc==4):
		pdos_char( sys.argv[1], sys.argv[2], sys.argv[3])
	if( argc==5):
		pdos_char( sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]))






