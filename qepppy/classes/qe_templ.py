from ._qe_templ_base_ import namelist as NL, card as CARD


"""
Template format:
{
	nl:   ['", ... (List of namelists name)]
	card: ['", ... (List of cards name)]
	NAMELIST_NAME1: {
		VAR_NAME: {
			v: (Value of the parameter)
			t: (Type of the parameter)
			d: (Default value)
			c: (List of possible acceppted value for the parameter)
			vec:None/(start,end) (Info for array like variables)
			}
		...
		}
	NAMELIST_NAME2={...}
	...

	CARD_NAME1:{
		v: (Value associated with the card)
		c: (List of possible acceppted value for the card)
		d: (Default value for v)
		r: True/False (is card REQUIRED?)
		u: True/False (True if any values are set in syntax)
		syntax:{
			cond: "..." (Condition on card value)
			l:[ [{n: varname, v: value, t: TYPE}, ..., [{...}, ...]], [...], ([{...}, {...}, ..., ], s, e, kw), ...]
				(Every element of the list represent a line
				A Tuple represent a repeating line ( [line], start, end, keyword)
				A List within a list marks optional arguments)
		}
		syntax1:{...} (if multiple syntaxes are provided)
	}
	CARD_NAME2:{...}
	...
}
"""


class qe_templ( CARD, NL):
	def find( self, *args, up=""):
		"""
		Find a var in all possible namelists and cards.
		Return its value if found otherwise return None
		"""
		def _syntax_find_( el, tof):
			#Recursive find to descend into syntax elements
			if not isinstance( el, list): return None
			for e in el:
				f = None
				if isinstance( e, dict):
					if e['n'] == tof: f = e['v']
				if isinstance( e, list):
					f = _syntax_find_( e, tof=tof)
				if isinstance( e, tuple): 
					f = _syntax_find_( e[0], tof=tof)
				if f != None: return f
			return None

		l=[]
		for name in args:
			n = None
			tof = name
			if "(" in name:
				tof = name.split( "(")[0]
				n = int( name.split( "(")[1].split( ")")[0])

			ret = None
			for nl in self._templ_['nl']:
				if up:
					if up != nl: continue
				f = self._templ_[nl].get( tof)
				if f: 
					if n: 
						try: ret = f['v'][n-1]
						except: pass
					else: ret = f['v']
					break
			for card in self._templ_['card']:
				if up:
					if up != card: continue
				synt = self._get_syntax_( self._templ_[card])
				f = _syntax_find_( synt, tof=tof)
				if f != None: 
					if n: 
						try: ret = f[n-1]
						except: pass
					else: ret = f
					break
			l.append( ret)

		if len( l) == 1: l = l[0]
		else: l = tuple( l)
		return l

	def load_templ( self, fname=""):
		import os
		if os.path.isfile( fname):
			with open( fname) as f:
				file = f.read()
		else:
			from pkg_resources import resource_string, resource_listdir
			#print( resource_listdir('qepppy.data', ''))
			if fname in resource_listdir('qepppy.data', ''):
				file = resource_string( 'qepppy.data', fname).decode('utf-8')

		import ast
		self._templ_ = ast.literal_eval( file)
		return




	



	




