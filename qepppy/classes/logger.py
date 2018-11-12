import inspect
import functools
import logging

__all__ = ['logger', 'debug', 'info','warning', 'error', 'critical', 'global_msg', 'global_thr']


#logging.basicConfig( format='%(levelname)s:\n%(name)s --- %(message)s.')
#log = logging.getLogger( __name__)
logging.basicConfig( format='%(levelname)s --- %(name)s --- %(message)s.')

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


global_msg = 'INFO'
global_thr = 'ERROR'


def logger_wrap( func, msg_lvl=global_msg, thr_lvl=global_thr, msg_name=__name__):
	#print( "Decorating function {}".format( func.__name__))
	@functools.wraps( func)
	def wrapped( *args, **kwargs):
		log = logging.getLogger( msg_name)
		log.setLevel( msg_lvl)
		#handler = logging.StreamHandler( )
		#handler.setLevel( 1)
		#log.addHandler( handler)
		levels ={
			'DEBUG':log.debug,
			'INFO':log.info,
			'WARNING':log.warning,
			'ERROR':log.error,
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
			try:
				msg = e.level.upper()
			except:
				raise e
			else:
				if logging.getLevelName( msg) >= logging.getLevelName( thr_lvl.upper()):
					raise e
				else:
					levels[msg]( str(e))
		else:
			return res
	return wrapped

def logger( msg_lvl=global_msg, thr_lvl=global_thr):
	def _log_( elem):
		#print( "Decorating element {}".format( elem))
		if inspect.isclass( elem):
			class NewCls( elem):
				def __init__( self, *args, **kwargs):
					self.__name__ = elem.__name__
					#print( "MRO --- ", self.__class__.__mro__)
					try:
						super().__init__( *args, **kwargs)
					except Exception as e:
						pass
					return
				def __getattribute__( self, s):
					#print( "Getting attribute --- ", s)
					x = super().__getattribute__( s)

					if inspect.ismethod( x):
						name = "{}:{}".format( self.__name__, x.__name__)
						return logger_wrap( x, 
							msg_lvl=msg_lvl, 
							thr_lvl=thr_lvl, 
							msg_name=name
							)
					else:
						return x
			return NewCls
		if inspect.isfunction( elem):
			name = "{}",format( elem.__name__)
			return logger_wrap( elem,
				msg_lvl=msg_lvl,
				thr_lvl=thr_lvl,
				msg_name=name
				)
	return _log_
