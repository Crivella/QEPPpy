import xml.etree.ElementTree as ET

class bands:
	__name__ = "bands"
	def __init__( self, fname=""):
		self.i = 0
		if not fname:
			raise Exception( "Must initialize class giving the name of the .xml file")
		if not os.path.isfile( fname):
			raise IOError( "File '{}' does not exist".format( fname))
		self.fname = fname
		self._parse_xml_( self)

	def __str__( self):
		bnd = len(self.egv[0])-1
		kpt_fmt = "\nkpt:\t" + "{:8.4f}"*3
		egv_fmt = "\nEigenvalues( eV):\n"+("\t"+"{:10.6f}"*8+"\n")*int(bnd/8)
		egv_fmt += "\t"+"{:10.6f}"*(bnd%8+1)+"\n"
		msg = ""
		for i in range( self.n_kpt):
			msg += kpt_fmt.format( *self.kpt[i])
			msg += egv_fmt.format( *self.egv[i])
		return msg

	def __getitem__( self, key):
		print( key)
		if( 0<=key<self.n_kpt):
			return { 'kpt':self.kpt[key], 'egv':self.egv[key], 'occ':self.occ[key]} 
		else:
			raise Exception( "Index '{}' out of range".format( key))

	def __iter__( self):
		return self

	def __next__( self):
		if( self.i < self.n_kpt):
			i = self.i
			self.i += 1
			return { 'kpt': self.kpt[i], 'egv': self.egv[i], 'occup': self.occ[i]}
		else:
			self.i = 0
			raise StopIteration()

	def __enter__( self):
		return self

	def __exit__( self, *args):
		del self

	def _parse_xml_( self):
		root = ET.parse( self.fname).getroot()
		self.kpt = list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.findall("output//ks_energies/k_point"))))
		self.weight = list( map( float, map( lambda x: x.get('weight'), \
			root.findall("output//ks_energies/k_point"))))
		self.egv = list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.findall("output//ks_energies/eigenvalues"))))
		self.occ = list( map( lambda y: list(map( float, y)), \
			map( lambda x: x.text.split(" "), root.findall("output//ks_energies/occupations"))))
		self.n_el    = float( root.find("output//nelec").text)

		self.lsda    = root.find("output//lsda").text
		self.noncolin = root.find("output//noncolin").text
		if( self.lsda == 'true'):
			self.n_spin = 2
		if( self.noncolin == 'true'):
			self.n_spin = 4








