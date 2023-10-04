import functools
import inspect
import json
import logging
import traceback

__all__ = ['logger', 'debug', 'info','warning', 'error', 'critical', 'function_wrap']

try:
    with open("log.json", "r") as f:
        content = json.load(f)
except Exception:
    content = {}

global_log_enabled = content.get('global_log_enabled', True)
global_log_level   = content.get('global_log_level',   'INFO')
global_log_thr     = content.get('global_log_thr',     'ERROR')

logging.basicConfig(
    format=(
        '-'*80 +
        '\n%(levelname)s ---> %(name)s  line:%(lineno)s' + 
        '\n\tMSG: "%(message)s"' + 
        '\n' + '-'*80
        ), 
    level=global_log_level
    )

global_logger    = logging.getLogger("global")
global_logger.propagate = False
ch = logging.StreamHandler()
global_formatter = logging.Formatter(
    '-'*80 + '\n%(levelname)s ---> %(message)s\n' + '-'*80)
ch.setFormatter(global_formatter)
global_logger.addHandler(ch)

def stack_analyzer(lvl=2):
    stack = inspect.stack()[lvl]
    module = stack[1].split("/")[-1]
    lineno = stack[2]
    func = stack[3]

    return "{}:{}:{}".format(module, func, lineno)


class critical(Exception):
    level = 'CRITICAL'
    @staticmethod
    def print(msg=""):
        stack_msg = stack_analyzer()
        global_logger.critical("{}\n{}".format(stack_msg, msg))

class error(Exception):
    level = 'ERROR'
    @staticmethod
    def print(msg=""):
        stack_msg = stack_analyzer()
        global_logger.error("{}\n{}".format(stack_msg, msg))

class warning(Exception):
    level = 'WARNING'
    @staticmethod
    def print(msg=""):
        stack_msg = stack_analyzer()
        global_logger.warning("{}\n{}".format(stack_msg, msg))

class info(Exception):
    level = 'INFO'
    @staticmethod
    def print(msg=""):
        stack_msg = stack_analyzer()
        global_logger.info("{}\n{}".format(stack_msg, msg))

class debug(Exception):
    level = 'DEBUG'
    @staticmethod
    def print(msg=""):
        stack_msg = stack_analyzer()
        global_logger.debug("{}\n{}".format(stack_msg, msg))

def get_original_class_name(method):
    if inspect.ismethod(method):
        for cls in inspect.getmro(method.__self__.__class__):
            if cls.__name__ == "NewCls":
                continue
            if cls.__dict__.get(method.__name__):
                return cls.__module__ + ':' + cls.__name__
    return method.__module__


def function_wrap(func, msg_lvl=global_log_level, thr_lvl=global_log_thr):
    msg_lvl = msg_lvl.upper()
    thr_lvl = thr_lvl.upper()
    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        msg_name = func.__name__
        try:
            msg_name = get_original_class_name(func) + ":" + msg_name
        except:
            pass
        log = logging.getLogger(msg_name)
        log.setLevel(msg_lvl)
        #handler = logging.StreamHandler()
        #handler.setLevel(1)
        #log.addHandler(handler)
        if logging.getLevelName(msg_lvl) <= logging.getLevelName('DEBUG'):
            log.debug("Calling function '{}' with args:{}, kwargs:{}".format(
                func.__name__, args, kwargs))
        try:
            res = func(*args, **kwargs)
        except Exception as e:
            try:
                msg = e.level
            except:
                if msg_lvl == 'DEBUG':
                    tb = e.__traceback__
                    tb = tb.tb_next
                    e.__traceback__ = tb
                    traceback.print_tb(tb)
                log.critical(e)
            else:
                if logging.getLevelName(msg) >= logging.getLevelName(thr_lvl):
                    if msg_lvl == 'DEBUG':
                        tb = e.__traceback__
                        tb = tb.tb_next
                        # Remove decorator signal form trace
                        while 'decorator' in tb.tb_frame.f_code.co_filename:
                            tb = tb.tb_next
                        e.__traceback__ = tb
                        traceback.print_tb(tb)
                    log.critical(e)
                else:
                    log.__getatrribute__(msg.lower())(str(e))
        else:
            return res
    wrapped.__name__=func.__name__
    return wrapped

def logger(msg_lvl=global_log_level, thr_lvl=global_log_thr, log_enabled=global_log_enabled):
    """
    Decorator generator that can be used to log functions or all class methods called through 
    __getatrribute__ (applied directly to the class)
    """
    def _log_(elem):
        if not log_enabled:
            return elem
        if inspect.isclass(elem):
            def new_getattr(cls, attr):
                x = object.__getattribute__(cls, attr)

                if inspect.ismethod(x):
                    return function_wrap(
                        x, 
                        msg_lvl=msg_lvl, 
                        thr_lvl=thr_lvl, 
                        )
                else:
                    return x

            elem.__getattribute__ = new_getattr
            return elem
        if inspect.isfunction(elem):
            return function_wrap(elem,
                msg_lvl=msg_lvl,
                thr_lvl=thr_lvl,
                )
    return _log_
