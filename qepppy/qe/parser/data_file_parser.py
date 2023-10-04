import os

import numpy as np

# from ...logger import logger, warning


def _format_(var):
    if var is None:
        return None
    if isinstance(var, np.ndarray):
        return var
    if isinstance(var, str):
        if var.upper() in ['TRUE', 'T', 'YES', 'Y']:
            return True
        elif var.upper() in ['FALSE', 'F', 'NO']:
            return False
    try:
        return int(var)
    except:
        try:
            return float(var)
        except:
            return var

matches = {
    int:   r'[\s]*(?P<flag>[\d\-]+)',
    float: r'[\s]*(?P<flag>[\d\.\-EedD]+)',
    str:   r'[\s]*(?P<flag>.*)',
    bool:  r'[\s]*(?P<flag>.)',
}

def _get_value_(f, search_data, dtype=str):
    import re

    string = search_data['outfile_regex']
    m = search_data.get('modifier', 1)

    a = None
    try:
        if dtype == list:
            a = re.finditer(string, f)
            a = [x.groupdict() for x in a]
            for n,e in enumerate(a):
                for k,v in e.items():
                    b = np.fromstring(v, sep=' ')
                    if len(b) == 0 or re.findall(r'\s[a-zA-Z]', v):
                        b = str(v).strip()
                    elif len(b) == 1:
                        b = b[0]
                    a[n][k] = b*m
        else:
            sstring = string
            if r'(?P<flag>' not in string:
                sstring += matches[dtype]
            a = re.search(sstring, f).group('flag')
    except Exception as e:
        # print(e)
        pass

    return dtype(a)

def _xml_attr_(node, f='', n=''):
    if f:
        node = node.find(f)
    return _format_(node.get(n))

def _xml_text_(node, f='', n=''):
    if f:
        node = node.find(f)
    return _format_(node.text)

def _xml_node_list_(node, f='', n=''):
    if f:
        node = node.findall(f)
    res = []
    for c in node:
        d = c.attrib
        d = {k:_format_(v) for k,v in d.items()}
        add = None
        if not c.text is None:
            add = c.text.strip().replace('\n', ' ')
        if add:
            add = list(filter(None, add.split(' ')))
            if len(add) == 1:
                add = add[0]
            if isinstance(add, list):
                add = np.array(add, dtype=float)
            if n:
                tag = n
            else:
                tag = c.tag
            d[tag] = _format_(add)
        else:
            for e in  _xml_node_list_(c.getchildren()):
                d.update(e)
        res.append(d)

    return res

# @logger()
class data_file_parser(object):
    """
    Parser for QE data'file'schema.xml (QE>=6.2) and pw.x outputs.

    - d:        Rule dictionary defining how the parser will operate.
    - xml:      Name of the "data-file*.xml" to parse or path containing the file
                "data-file-schema.xml" (Parsing will run if the schema is set).
                When schema is given will also genetate 2 attributes:
                - data_path: Path of the folder containing the *.save folder.
                - prefix:    Prefix of calculation
    - outfile:  Name of the pw.x output file to parse (Parsing will run if the
                outfile is set).
    - **kwargs: Overwrite parsed variables with user's one given as a
                var_name=value keyword arg.
                The var_name must already exist after the parsing (cannot set
                unrecognized variables).
    This parser accept a rule dict passed as a keyword argument 'd' to the
    __init__ method.
    The rule dict has to follow the format:
    {
    'varname':{
        'res_type':type, (int/float/ndarray/...)
        'xml_search_string':element to find, (string)
        'xml_ptype':xmlacquisition rule, ('attr'/text/nodelist)
        'extra_name':name
        'outfile_regex': regex string to extract data from outfile
        }
    'varname1':{...}
    }
    The parsing will generate internal variables using as name the keys of the
    rule dict.
    The rule are as follow:
    - attr:     Get the value of the attribute of name "n" of the node given
                by root.find(f)
    - text:     Get the text of the node given by root.find(f)
    - nodelist: Find a list of nodes using root.findall(f) and analyze them.
                The resulting variable will be a list dictionary (one element
                for every node found).
                "extra_name" is used as a key name for the text value of the
                node, otherwise the xml node.tag is used.
                Text that contains arrays of number are automatically converted
                into np.ndarray objects.
                All the node attributes are saved are saved as key:values in the
                dict.
                If the node contains children instead of text, explore them
                recursively:
                - saves all the attributes found in the dict as key:value pairs.
                - saves all the the text values found in the dict as
                   node.tag:text pairs.
                NOTE: No conflict resolution is present for attributes with the
                      same name across children.
                      The attribute are updated at every step, so the value of
                      the last one will be stored.
    """
    __name__ = 'data_file_parser'
    def __init__(self, xml='', outfile='', data={}, **kwargs):
        self._data_ = data
        for i in data:
            self.__dict__[i] = None
        if xml:
            self.xml = os.path.realpath(xml)
            self._set_data_file_()
            self._parse_xml_()
        elif outfile:
            self.outfile = outfile
            self._parse_outfile_()
        if kwargs:
            for k, v in kwargs.items():
                if not k in self.__dict__:
                    continue
                self.__dict__[k] = v

    # def __getitem__(self, key):
    #     return self.__dict__.get(key)

    # def __str__(self):
    #     return ""

    def _set_data_file_(self):
        if os.path.isfile(self.xml):
            self.data_path = os.path.dirname(os.path.realpath(self.xml))
        elif os.path.isdir(self.xml):
            file = os.path.join(self.xml, 'data-file-schema.xml')
            if not os.path.isfile(file):
                raise FileNotFoundError(file)
            self.data_path = self.xml
            self.xml    = file
        else:
            raise FileNotFoundError(f'The file/folder {self.xml} does not exist')
        if 'data-file-schema.xml' in self.xml:
            self.prefix    = '.'.join(os.path.basename(self.data_path).split('.')[:-1])
            self.data_path = os.path.abspath( os.path.join(self.data_path, os.path.pardir))

    def _parse_outfile_(self):
        with open(self.outfile, 'r') as f:
            content = f.read()

        for k, v in self._data_.items():
            t = v['res_type']
            search = v.get('outfile_regex', None)
            if search is None:
                continue
            val = None
            try:
                val = _get_value_(content, v, dtype=t)
            except Exception as e:
                # print("ERROR: ", k, e)
                pass
            self.__dict__[k] = val
        return

    def _parse_xml_(self):
        import xml.etree.ElementTree as ET
        root = ET.parse(self.xml).getroot()

        xml_acq_rule={
            'attr':     _xml_attr_,
            'text':     _xml_text_,
            'nodelist': _xml_node_list_,
        }

        for k, v in self._data_.items():
            res = None
            t = v['res_type']
            try:
                n    = v.get('extra_name', None)
                f    = v['xml_search_string']
                func = xml_acq_rule[v['xml_ptype']]
                res  = func(root, f, n)
            except Exception as e:
                continue

            if res is None:
                self.__dict__[k] = None
            elif t == np.array:
                self.__dict__[k] = np.fromstring(res, sep=' ', dtype=float)
            else:
                self.__dict__[k] = t(res)
        return

    def validate(self):
        return True
        # try:
        #     return True and super().validate()
        # except:
        #     return True
