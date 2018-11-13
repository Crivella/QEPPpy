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
		for cls in inspect.getmro( method.__self__.__class__):
			if cls.__name__ == "NewCls":
				continue
			if cls.__dict__.get( method.__name__):
				return cls.__name__
	return "__root__"


def function_wrap( func, msg_lvl=global_log_level, thr_lvl=global_log_thr):
	@functools.wraps( func)
	def wrapped( *args, **kwargs):
		msg_name = func.__name__
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
			res = func( *args, **kwargs)
		except Exception as e:
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
		if inspect.isclass( elem):
			class NewCls( elem):
				@function_wrap
				def __init__( self, *args, **kwargs):
					try:
						super().__init__( *args, **kwargs)
					except Exception as e:
						raise e
					return

				def __getattribute__( self, s):
					x = object.__getattribute__( self, s)

					if inspect.ismethod( x):
						return function_wrap( x, 
							msg_lvl=msg_lvl, 
							thr_lvl=thr_lvl, 
							)
					else:
						return x
			return NewCls
		if inspect.isfunction( elem):
			return function_wrap( elem,
				msg_lvl=msg_lvl,
				thr_lvl=thr_lvl,
				)
	return _log_
