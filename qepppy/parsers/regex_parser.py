import re

import numpy as np

matches = {
    int:   r'[\s]*(?P<flag>[\d\-]+)',
    float: r'[\s]*(?P<flag>[\d\.\-EedD]+)',
    str:   r'[\s]*(?P<flag>.*)',
    bool:  r'[\s]*(?P<flag>.)',
}

# def dbg1(*args, **kwargs):
#     # return
#     print(*args, **kwargs)

class Parser_regex():
    """
    Parse a file using regex defined syntax. Meant to be used as a child class.
    The parent class can be passed a `regex_data` on its init to specify the
    syntax:
     - KEYs: name of the attribute/s to be set in the class after parsing the
       file. Can use comma to specify more than one attribute to be set as a
       result of the same parse operation. '_' can be used in order to ignore
       part of the results
       EG:
       - {'test':...} will set an attribute `test` on the parent class
       - {'test1,test2':...} will set attributes `test1` and `test2`.
       - {'p1,_,p2':...} will set attributes `p1` and `p2`, while ignoring the
         second variable resulting from the parse operation.
       For more info on how the result is split from the parsing see the parse
       rules section.
       with the result of the parsing rule given by the value.
     - VALUEs: dictionary specifing a set of rules to be used by the parser:
       - rstring: REGEX string used to parse the file
       - typ: specify the expected resulting type for the attribute.
              Several types are implemented for a 'single_value' or 'iterable'
              reult. If a non-implemented type is used a NotImplemented
              exception will be risen
              - single_value: int/float/str/bool
                This search expect the REGEX string to contain a 'flag' named
                search group, the content of which will be the result of the
                parse. If 'flag' is not present, a default search group will
                be added depending on the specified type.
                EG: rstring='nbnd = ', typ='int
                    final_rstring = 'nbnd = [\\s]*(?P<flag>[\\d\-]+)'
              - iterable: list/np.ndarray
                In this case rstring can specify multiple search groups, the
                result of which will be made into a list.
                - If only an attribute name is given, the result will be a 2D
                  list/array, where the each row is the result of a `finditer`
                  and each column will correspond to a search group.
                - If multiple names are given, every column will be assigned
                  to an attribute in the corresponding order.
                  The number of columns/rgx_groups has to be equal to the number
                  of attribute names, or max_num has to be used to truncate
                  the resulting array.
       - repeatable: True/Fale (default) Allow to read output of calculations
         that print the same output multiple times for different iterations
       - max_num: MAX number of column to consider when parsing iterables
         Can be positive/negative to truncate starting from the left/right
         edge.
       - re_scale_fact: Multiply the final result by this scale factor.
         If this is a string, the parser will look for the corresponding
         attribute and get its value (useful when the factor has to be parsed
         as well)
    """
    def __init__(self, *args, file=None, regex_data={}, **kwargs):
        super().__init__(*args, **kwargs)

        self.regex_data = regex_data

        if file:
            self.parse_file(file)

    @property
    def regex_data(self):
        return self._regex_data

    @regex_data.setter
    def regex_data(self, value):
        if not isinstance(value, dict):
            raise ValueError('Value must be a dictionary.')

        self._regex_data = value

    @staticmethod
    def get_single_val(content, rstring, typ, repeat):
        res = None

        if not 'flag' in rstring:
            rstring += matches[typ]
        val = list(re.finditer(rstring, content))
        if not val:
            return
        if not repeat:
            res = typ(val[-1].group('flag'))
        else:
            res = [typ(_.group('flag')) for _ in val]

        return res

    @staticmethod
    def get_list_val(content, rstring):
        rgx = re.compile(r'\s[a-zA-Z]')

        res = re.finditer(rstring, content)
        res = [x.groupdict() for x in res]
        for n,e in enumerate(res):
            for k,v in e.items():
                if '*' in v or rgx.search(v) or not v.strip():
                    b = str(v).strip()
                else:
                    v = v.replace('-', ' -')
                    b = []
                    try:
                        b = np.fromstring(v, sep=' ')
                    except ValueError:
                        pass
                    if len(b) == 0 :
                        b = str(v).strip()
                    elif len(b) == 1:
                        b = b[0]
                res[n][k] = b

        return res

    # @staticmethod
    def max_num_truncate(self, max_num, val, repeat):
        fact = 1
        if max_num is None:
            return val

        if isinstance(max_num, str):
            if max_num.startswith('-'):
                fact = -1
                max_num = max_num[1:]
            max_num = getattr(self, max_num.strip())
        if not repeat:
            max_num *= fact
            if max_num > 0:
                val = val[:max_num]
            else:
                val = val[max_num:]
        else:
            sx = val.shape[0]
            if sx % max_num != 0:
                raise ValueError(f'Read data cannot be reshaped using given {max_num=}')
            new = sx//max_num
            if new > 1:
                val = val.reshape(sx//max_num, max_num, -1)

        return val

    def scale(self, val, fact):
        if not fact is None:
            if isinstance(fact, str):
                fact = getattr(self, fact)
            if len(val)>0 and isinstance(val.flatten()[0], (int,float,np.number)):
                val *= fact

        return val


    def parse_file(self, file):
        # res = {}

        with open(file, 'r') as f:
            content = f.read()

        for k,v in self.regex_data.items():
            if not 'rstring' in v:
                continue
            # dbg1('-'*40)
            # dbg1('Key:', k)
            # dbg1('val:', v)

            rstring    = v.get('rstring')
            typ        = v.get('typ')
            max_num    = v.get('max_num', None)
            scale_fact = v.get('re_scale_fact', None)
            repeat     = v.get('repeatable', False)

            if typ in matches:
                val = self.get_single_val(content, rstring, typ, repeat)
                setattr(self, k, val)
                # dbg1('found_single:', val)
            elif typ in (list,np.ndarray):
                app = self.get_list_val(content, rstring)
                params = k.split(',')
                # dbg1(params)
                # dbg1(app)
                for num,name in enumerate(params):
                    if name == '_':
                        continue
                    n_elem  = len(app)
                    n_group = len(app[0].values()) if n_elem > 0 else 0
                    if len(params) == 1 and not (n_elem > 0 and n_group == 1):
                        app2 = np.array([list(a.values()) for a in app])
                    else:
                        app2 = np.array([list(a.values())[num] for a in app])

                    val = self.max_num_truncate(max_num, app2, repeat)
                    val = self.scale(val, scale_fact)

                    # dbg1(f'found_mul({name}):', val)
                    setattr(self, name, val)
            else:
                raise NotImplementedError()
