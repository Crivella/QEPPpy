import sys
import inspect
import functools
import logging
import traceback
import json

__all__ = ['logger', 'debug', 'info','warning', 'error', 'critical', 'function_wrap']

def get_cfg_var( name, c, fallback=None):
	if isinstance( c, dict):
		return c.get( name, fallback)
	return fallback

try:
	with open( "log.json", "r") as f:
		content = json.load( f)
except Exception as e:
	content = None

global_log_enabled = get_cfg_var( 'global_log_enabled', content, True)
global_log_level   = get_cfg_var( 'global_log_level'  , content, 'INFO')
global_log_thr     = get_cfg_var( 'global_log_thr'    , content, 'ERROR')

#print( global_log_enabled)
#print( global_log_level)
#print( global_log_thr)

#logging.basicConfig( format='%(levelname)s:\n%(name)s --- %(message)s.')
#log = logging.getLogger( __name__)
logging.basicConfig( format='%(levelname)s --- %(name)s --- %(message)s.', level=global_log_level)

class critical( Exception):
	level = 'critical'
	@staticmethod
	def print( msg=""):
		logging.critical( msg)

class error( Exception):
	level = 'error'
	@staticmethod
	def print( msg=""):
		logging.error( msg)

class warning( Exception):
	level = 'warning'
	@staticmethod
	def print( msg=""):
		logging.warning( msg)

class info( Exception):
	level = 'info'
	@staticmethod
	def print( msg=""):
		logging.info( msg)

class debug( Exception):
	level = 'debug'
	@staticmethod
	def print( msg=""):
		logging.debug( msg)

def get_original_class_name( method):
	if inspect.ismethod( method):
		#print( "test {}  \n".format( method))
		#print( inspect.getmro( method.__self__.__class__))
		for cls in inspect.getmro( method.__self__.__class__):
			if cls.__name__ == "NewCls":
				continue
			#print( cls, cls.__dict__.get( method.__name__))
			if cls.__dict__.get( method.__name__):
				#print( cls.__name__)
				return cls.__name__
	return "__root__"


def function_wrap( func, msg_lvl=global_log_level, thr_lvl=global_log_thr):
	#print( "Decorating function {}".format( func.__name__))
	@functools.wraps( func)
	def wrapped( *args, **kwargs):
		#print( "I am the wrapped function of {}".format( func.__name__))
		#print( "-----------", msg_name)
		msg_name = func.__name__
		#print( "Inspecting function {}".format( func.__name__))
		try:
			msg_name = get_original_class_name( func) + ":" + msg_name
		except:
			pass
		log = logging.getLogger( msg_name)
		log.setLevel( msg_lvl)
		#handler = logging.StreamHandler( )
		#handler.setLevel( 1)
		#log.addHandler( handler)
		levels ={
			'DEBUG':   log.debug,
			'INFO':    log.info,
			'WARNING': log.warning,
			'ERROR':   log.error,
			'CRITICAL':log.critical
		}
		if logging.getLevelName( msg_lvl) <= logging.getLevelName( 'DEBUG'):
			log.debug( "Calling function '{}' with args:{}, kwargs:{}".format(
				func.__name__, args, kwargs))
		try:
			#print( "Calling function {}({},{})".format( func.__name__, args, kwargs))
			res = func( *args, **kwargs)
			#print( "Returning {}".format( res))
		except Exception as e:
			"""
			print( "Trace before calling handling ---------------")
			import sys
			exc_type, exc_value, exc_traceback = sys.exc_info()
			traceback.print_tb( e.__traceback__)
			print( "")
			traceback.print_tb( exc_traceback)
			"""
			#handling( e, thr_lvl, levels)
			try:
				msg = e.level.upper()
			except:
				if msg_lvl.upper()=='DEBUG':
					tb = e.__traceback__
					tb = tb.tb_next
					e.__traceback__ = tb
					traceback.print_tb( tb)
				log.critical( e)
			else:
				if logging.getLevelName( msg) >= logging.getLevelName( thr_lvl.upper()):
					if msg_lvl.upper()=='DEBUG':
						tb = e.__traceback__
						tb = tb.tb_next
						e.__traceback__ = tb
						traceback.print_tb( tb)
					log.critical( e)
				else:
					levels[msg]( str(e))
		else:
			return res
	return wrapped

def logger( msg_lvl=global_log_level, thr_lvl=global_log_thr, log_enabled=global_log_enabled):
	"""
	Decorator that can be used to log functions or all class methods called thorugh __getatrribute__
	(applied directly to the class)
	"""
	def _log_( elem):
		if not log_enabled:
			return elem
		#print( "Decorating element {}".format( elem))
		if inspect.isclass( elem):
			class NewCls( elem):
				#"""
				@function_wrap
				def __init__( self, *args, **kwargs):
					#self.__name__ = elem.__name__
					#print( "MRO --- ", self.__class__.__mro__)
					try:
						super().__init__( *args, **kwargs)
					except Exception as e:
						#handling( e, thr_lvl)
						raise e
					return
				#"""

				def __getattribute__( self, s):
					#print( "Getting attribute --- ", s)
					x = object.__getattribute__( self, s)

					if inspect.ismethod( x):
						#name = "{}:{}".format( elem.__name__, x.__name__)
						return function_wrap( x, 
							msg_lvl=msg_lvl, 
							thr_lvl=thr_lvl, 
							#msg_name=name
							)
					else:
						return x
			return NewCls
		if inspect.isfunction( elem):
			#name = "{}".format( elem.__name__)
			#if inspect.ismethod( elem):
			#	print( "Inspecting a method of {}".format( elem.__self__.__name__))
			#	name = elem.__self__.__name__ + name
			return function_wrap( elem,
				msg_lvl=msg_lvl,
				thr_lvl=thr_lvl,
			#	msg_name=name
				)
	return _log_
