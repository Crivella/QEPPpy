# from abc import ABCMeta

from ...parsers import fortran_namelist_collection as fnc

tf90_to_py = {
    'INTEGER': int,
    'REAL': float,
    'LOGICAL': bool,
    'CHARACTER': str
    }

def key_getter(key, item=None, attr=None):
    def getter(cls):
        if not item is None:
            items = item.replace(' ', '').split(',')
            res = []
            for i in items:
                res.append(cls[i])
            if len(res) == 1:
                res = res[0]
            return res
        elif not attr is None:
            return getattr(cls, attr)
        else:
            raise ValueError(f'Cannot retrieve {key}')

    return getter

def key_setter(key, item=None, attr=None):
    def setter(cls, value):
        if not item is None:
            items = item.replace(' ', '').split(',')
            try:
                iter(value)
            except TypeError:
                value = [value]
            for i,v in zip(items,value):
                cls[i] = v
        elif not attr is None:
            setattr(cls, attr, value)
        else:
            raise ValueError(f'Cannot set {key}')

    return setter

class VariableLinker(type):
    def __new__(mcs, clsname, base, dct):
        new_dct = {}
        for k,v in dct.items():
            new_dct[k] = v
            if not isinstance(k, str) or not k.startswith('_link_') or not isinstance(v, dict):
                continue

            k       = k[6:]
            attrib  = v.pop('attrib',    None)
            item    = v.pop('item',      None)

            if attrib and item:
                raise ValueError('Must link the variable to only an attibute or an item, not both.')

            getter = key_getter(k, item, attrib)
            setter = key_setter(k, item, attrib)


            new_dct[k] = property(getter, setter)

        res = super().__new__(mcs, clsname, base, new_dct)

        return res



# class MixMeta(ABCMeta, VariableLinker):
#     """Needed for subclassing this togheter with another child that is a a child of abc.ABC"""
#     pass

class qe_input(fnc, metaclass=VariableLinker):
    cards = []

    def parse(self, *args, **kwargs):
        super().parse(*args, **kwargs)

        content = self.unparsed

        for card in self.cards:
            ind = content.find(card)
            if ind == -1:
                continue
            read_func = getattr(self, 'read_' + card.lower())
            read_func(content[ind:])

    @staticmethod
    def get_mode(content):
        l   = content.strip().split('\n')[0]
        val = None
        if '{' in l:
            val = l.split('{')[1].split('}')[0]
        elif ' ' in l:
            val = l.split(' ')[1]

        return val
