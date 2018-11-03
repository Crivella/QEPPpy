from ._qe_templ_base_ import namelist as NL, card as CARD


"""
Template format:
{
    nl:   ["...", ... (List of namelists name)]
    card: ["...", ... (List of cards name)]
    NAMELIST_NAME1: { (dictionary containing all the namelist parameters)
        VAR_NAME: { (dictionary containing the details of the parameter)
            v: (Value of the parameter)
            t: (Type of the parameter as a string)
            d: (Default value)
            c: (List of possible acceppted value for the parameter)
            vec:None/(start,end) (Info for array like variables e.g. celldm(1/2/3/4/5/6))
            }
        ...
        }
    NAMELIST_NAME2={...}
    ...

    CARD_NAME1:{ (Dictionary containing the data and syntax of a card)
        v: (Value associated with the card)
        c: (List of possible acceppted value for the card)
        d: (Default value for v)
        r: True/False (is card REQUIRED?)
        u: True/False (True if any values are set in syntax)
        syntax:{ (Dictionary defining the syntax that the card should follow)
            cond: "..." (Condition on card value)
            l:[ [{n: varname, v: value, t: TYPE}, ..., [{...}, ...]], [...], ([{...}, {...}, ..., ], s, e, kw), ...]
                Every element of the list represent a line
                A Tuple represent a repeating line ( [line], start, end, keyword)
                A List within a list marks optional arguments
        }
        syntax1:{...} (if multiple syntaxes are provided)
    }
    CARD_NAME2:{...}
    ...
}
"""


class qe_templ( CARD, NL):
	"""
	Instance for handling the QE namelist dict template in the data folder.

	This class provides the following methods:
	    validate(): Check namelist template after reading
	    convert():  Convert the internal dict in a string QE input file
	    check_nl( nl="namelist"): Check if nl is valid (present in the internal namelist)
	    set_nl:( nl="namelist", k="param", v="value to set") Set a namelist value in the namelist template
	    set_card: ( card="", v="", el=[]) Set a card value in the namelist template line by line
	    find: (name) Find a variable (or list ) with name=name in the namelist template

	Template format:
	{ (Template dictionary)
	    nl:   ["...", ... (List of namelists name)]
	    card: ["...", ... (List of cards name)]
	    NAMELIST_NAME1: { (dictionary containing all the namelist parameters)
	        VAR_NAME: { (dictionary containing the details of the parameter)
	            v: (Value of the parameter)
	            t: (Type of the parameter as a string)
	            d: (Default value)
	            c: (List of possible accepted value for the parameter)
	            vec:None/(start,end) (Info for array like variables e.g. celldm(1/2/3/4/5/6))
	            }
	        ...
	        }
	    NAMELIST_NAME2={...}
	    ...
	    CARD_NAME1:{ (Dictionary containing the data and syntax of a card)
	        v: (Value associated with the card)
	        c: (List of possible accepted value for the card)
	        d: (Default value for v)
	        r: True/False (is card REQUIRED?)
	        u: True/False (True if any values are set in syntax)
	        syntax:{ (Dictionary defining the syntax that the card should follow)
	            cond: "..." (Condition on card value)
	            l:[ [{n: varname, v: value, t: TYPE}, ..., [{...}, ...]], [...], ([{...}, {...}, ..., ], s, e, kw), ...]
	                Every element of the list represent a line
	                A Tuple represent a repeating line ( [line], start, end, keyword)
	                A List within a list marks optional arguments
	        }
	        syntax1:{...} (if multiple syntaxes are provided)
	    }
	    CARD_NAME2:{...}
	    ...
	}
	"""
	def find( self, *args, up=""):
		"""
		*args: sequence of names of param to find.
		up: If not set, look for the names in all possible NAMELISTs and CARDs.
		    If set, look only in the NAMELIST or CARD with the specified name.
		Return a tuple of the found values (None if not found)
		NOTE: Looking for array variables (e.g. celldm) without specifying the index (celldm(1)),
		      will return a list of all the values for that param.
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
				if card == tof:
					return self._templ_[card]['v']
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
		"""
		Load a QE template from a specified file (fname) or the internal data
		"""
		import os
		if os.path.isfile( fname):
			with open( fname) as f:
				file = f.read()
		else:
			from pkg_resources import resource_string, resource_listdir
			#print( resource_listdir('qepppy.data', ''))
			if fname in resource_listdir('qepppy.data', ''):
				file = resource_string( 'qepppy.data', fname).decode('utf-8')

		#import json
		#self._templ_ = json.loads( file)
		#return

		import ast
		self._templ_ = ast.literal_eval( file)
		return

	def validate( self):
		"""
		Validate the template.
		Used to check a read input file.
		NAMELIST:
		- checks that all required param are set
		- checks that the TYPE of all set param is conform with the template
		- checks that all set param are within the possible values (if a list is given in the template)
		CARD:
		- checks that all required CARDs are set
		- check that CARD values (the param next to the card name) are within list of possible values
		- check the type of every assigned element
		- check that all the lines are properly filled (groups of param required/optional)
		  have to be all assigned is one of them is.
		"""
		return super().validate()





	



	




