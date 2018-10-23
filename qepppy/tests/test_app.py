from unittest import TestCase


class Test( TestCase):
	def test_123( self):
		print( "Hello world of tests!!!")

		s=""
		self.assertTrue(isinstance(s, int))
