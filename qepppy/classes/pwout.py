import os.path
import xml.etree.ElementTree as ET

import numpy as np

from qepppy.classes.bands     import bands     as bands
from qepppy.classes.structure import structure as structure
#from .eigenv  import egv

class pwout( ):
	__name__ = "pwout"
	def __init__( self, fname=""):
		if not fname:
			raise Exception( "Must initialize class giving the name of the .xml file")
		if not os.path.isfile( fname):
			raise IOError( "File '{}' does not exist".format( fname))
		if ".xml" not in fname:
			raise IOError( "File extension must be '.xml'")
		self.fname = fname
		self.bnd   = bands( fname)
		self.stc   = structure( fname)
		self.smallest_gap = self.smallest_gap
		#self._parse_xml_()

	def __getattr__( self, key):
		if key in self.bnd.__dict__:
			return self.bnd.__getattr__( key)
		if key in self.stc.__dict__:
			return self.stc.__getattr__( key)
		raise AttributeError( "'{}' object does not have attribute '{}'".format( self.__name__, key))

	def __str__( self):
		msg = ""
		msg += self.bnd.__str__( )
		msg += self.stc.__str__( )
		return msg

	def __getitem__( self, key):
		if key in self.__dict__:
			return self.__dict__[key]
		if( isinstance( key, int)):
			return self.bnd.__getitem__( key)
		if( isinstance( key, str)):
			if key in self.bnd.__dict__:
				return self.bnd.__getitem__( key)
			if key in self.stc.__dict__:
				return self.stc.__getitem__( key)
		raise KeyError( "'{}' object does not support key '{}'".format( self.__name__, key))

	def __iter__( self):
		return self.bnd.__iter__( )

	def __enter__( self):
		return self

	def __exit__( self, *args):
		del self

	def band_structure( self, fname="plotted.dat", plot=True):
		kpt   = np.array( self.kpt)
		n_bnd = self.n_bnd
		x = [0]
		for i in range( 1, len( kpt)):
			x.append( np.linalg.norm( kpt[i] - kpt[i-1]))
			x[i] += x[i-1]
		#mod = list( map( np.linalg.norm, self.kpt))
		#x   = list( map( lambda x: sum( mod[0:mod.index(x)]), mod))
		x = np.array( x)
		print( x)
		egv = self.egv

		#print( list( zip ( x, egv)))
		#res = np.array( list( zip ( x, egv)))
		#res = np.concatenate( [x, egv], axis=1)
		res = np.column_stack( ( x, egv))
		#res = hstack( (x, egv))
		print( res[:,0], "\n\n", res[:,1:])

		np.savetxt( fname=fname, X=res, fmt="%13.8f"+"".join('%11.6f'*n_bnd))

		if plot:
			import matplotlib.pyplot as plt
			from matplotlib.ticker import AutoMinorLocator as AML
			fig, ax = plt.subplots()
			plt.plot( res[:,0], res[:,1:])
			plt.ylabel( "Energy( eV)")
			plt.xlabel( "")
			ml1 = AML(5)
			ax.yaxis.set_minor_locator(ml1)
			ax.yaxis.set_tick_params(which='both', right = True)
			plt.legend()
			plt.show()
		


		return 0

	def smallest_gap( self, radius=0., comp_point=[0.,0.,0.]):
		'''
		print( args)
		argc = len( args)
		if( argc > 0):
			if( argc == 1):
				radius     = args[0]
			if( argc == 2):
				radius     = args[0]
				comp_point = args[1]
			if( argc > 2):
				raise Exception( "This function accepts a maximum of 2 positional paramteters\n" \
						 "radius and comp_point")
		'''
		print( "SMALLEST_GAP: radius={}, comp_point={}".format( radius, comp_point))

		if( radius < 0):
			print( "Invalid negative 'radius'")
			return 1
		if( len(comp_point) != 3):
			print( "'comp_point' should be an [x,y,z] vector {}".format( comp_point))
			return 1
		for x in comp_point:
			try:
				x+=1
			except TypeError:
				print( "Invalid type %s for 'comp_point'" % type(x))
				return 1

		print( "E_fermi(from file):\t{:f} eV".format( self.fermi))

		kpt   = self.kpt
		n_kpt = self.n_kpt
		n_el  = self.n_el
		ef    = self.fermi
		if ( self.lsda):
			n_el /= 2
		if ( self.noncolin):
			vb = (n_el - 1)
			print( "spin-orbit correction detected");
		else:
			vb = (n_el/2 - 1)
			print( "No spin-orbit correction")
		vb = int(vb)
		cb = vb + 1
		print( "vb = %d, cb = %d" % (vb+1, cb+1));

		base = np.array( self.egv)
		cp = np.array( comp_point)
		if( radius > 0):
			mod = map( lambda x: radius - np.linalg.norm(cp-x), self.kpt)
		else:
			mod = map( lambda x: 1, self.kpt)
		lcb  = base[:,cb]
		lvb  = base[:,vb]
		lcb1 = base[:,cb+1]
		lg1  = lcb - lvb

		lll  = list( filter( lambda x: x[1]>= 0, zip( range( n_kpt), mod, kpt, lvb, lcb, lcb1)))
		found = len( lll)
		if( found <= 0):
			raise Exception( "No k-point found for the given criteria")
		print( "\nFound {} points with the given criteria.".format( found))

		res = mg1 = max( lll, key = lambda a: a[3])
		print( "\nMax_vb_energy: vb= {:f} eV".format( res[3]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))

		res = mg2 = min( lll, key = lambda a: a[4])
		print( "\nMin_cb_energy: cb= {:f} eV".format( res[4]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))

		if( mg1 == mg2):
			print( "DIRECT GAP %lf eV" % (m2 - m1))
		else:
			if( mg2[4] < mg1[3]):
				print( "METALLIC")
			else:
				print( "INDIRECT GAP %lf eV" % (mg2[4] - mg1[3]))

		if( ef == ef):
			res = min( (i for i in lll if i[3] < ef < i[4]), key = lambda a: a[4] - a[3])
			#mog = l[ lg1.index( m)]
			print( "\nMin_opt_gap: {:f} eV".format( res[4] - res[3]))
			print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1})".format(res[2], res[0]+1))
			print("\t%lf -> %lf   Ef: %lf eV" % (res[3], res[4], ef))		
		else:
			print( "\nCannot calculate min_opt_gap with invalid fermi energy")

		res = min( lll, key = lambda a: a[4] - a[3])
		print( "\nMin_gap_energy (vb->cb): {:f} eV".format( res[4] - res[3]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))
		print("\t%lf -> %lf   Ef: %lf eV" % (res[3], res[4], ef))

		res = min( lll, key = lambda a: a[5] - a[3])
		print( "\nMin_gap_energy (vb->cb+1): {:f} eV".format( res[5] - res[3]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))
		print("\t%lf -> %lf   Ef: %lf eV" % (res[3], res[5], ef))	

		res = min( lll, key = lambda a: a[5] - a[4])
		print( "\nMin_gap_energy (cb->cb+1): {:f} eV".format( res[5] - res[4]))
		print("\tat {0[0]:.6f} {0[1]:.6f} {0[2]:.6f} (2pi/a) (# {1}) ".format( res[2], res[0]+1))
		print("\t%lf -> %lf   Ef: %lf eV" % (res[4], res[5], ef))	

		return 0

		


if __name__ == "__main__":
	import sys
	import ast
	argc = len( sys.argv)
	if( not 3<=argc):
		print("Incorrect use. Pleas pass arguments:"
			"\n\t'filename',"
			"\n\t'function name',"
			"\n\t'function params\t(optional)'")
		exit()

	data = pwout( sys.argv[1])
	func = data[sys.argv[2]]
	if( sys.argv[2] == 'smallest_gap'):
		if( argc == 4):
			func( radius = float( sys.argv[3]))
		if( argc == 5):
			func( radius = float( sys.argv[3]), comp_point=ast.literal_eval( sys.argv[4]))
		if( argc > 5):
			raise Exception( "Invalid number of arguments for function 'smallest_gap'")


















