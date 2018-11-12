import inspect
import functools
import logging

__all__ = ['logger', 'debug', 'info','warning', 'error', 'critical', 'global_msg', 'global_thr']


logging.basicConfig( format='%(levelname)s:\n%(name)s --- %(message)s.')
#log = logging.getLogger( __name__)
#logging.basicConfig( format='%(levelname)s:\n%(name)s --- %(message)s.')

class critical( Exception):
	level = 'critical'
class error( Exception):
	level = 'error'
class warning( Exception):
	level = 'warning'
class info( Exception):
	level = 'info'
class debug( Exception):
	level = 'debug'


global_msg = 'INFO'
global_thr = 'ERROR'


def logger_wrap( func, msg_lvl=global_msg, thr_lvl=global_thr, msg_name=__name__):
	print( "Decorating function {}".format( func.__name__))
	@functools.wraps( func)
	def wrapped( *args, _log_level_=msg_lvl, **kwargs):
		log = logging.getLogger( msg_name)
		logging.basicConfig( level=msg_lvl)
		levels ={
			'DEBUG':log.debug,
			'INFO':log.info,
			'WARNING':log.warning,
			'ERROR':log.error,
			'CRITICAL':log.critical
		}
		try:
			print( "Calling function {}({},{})".format( func.__name__, args, kwargs))
			res = func( *args, **kwargs)
			print( "Returning {}".format( res))
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
		print( "Decorating element {}".format( elem))
		if inspect.isclass( elem):
			class NewCls( object):
				def __init__( self, *args, **kwargs):
					self.oInstance = elem( *args, **kwargs)

				def __getattribute__( self, s):
					try:
						x = super( NewCls, self).__getattribute__(s)
					except AttributeError:
						pass
					else:
						return x
					x = self.oInstance.__getattribute__(s)
					if inspect.ismethod( x):
						name = "{}:{}".format( elem.__name__, x.__name__)
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
