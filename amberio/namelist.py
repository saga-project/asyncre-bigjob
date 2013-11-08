"""
Support for Fortran namelists (nl).

This module provides classes and routines for reading, writing, and 
manipulating Fortan namelists. The predominant focus is on "new-style" 
F90 conventions, but this should also be compatible with F77 conventions 
in most cases.

A Fortran namelist is comprised of:
    - an ampersand followed by the namelist name
    - zero or more name-value subsequences separated by a value separator
    - a terminating slash

Example:
    &input
     name1 = value1,
     name2 = value2
    /

Exported Classes:
    Namelist           Essentially a named dict with additional 
                       attributes for string formatting.
    NamelistCollection A list-derived class containing only Namelists. 
                       Also has matching utilities using the Namelist 
                       "name" attribute.

Exported Functions:
    read_namelists     Return a NamelistCollection composed of all 
                       namelists in a file.
    read_non_namelists Return non-namelist parts of a file as a list of 
                       lines.
    separate_namelists Return output of read_namelists and 
                       read_non_namelists as a tuple.

Example Usage:

--nl_sample--
 &sample_a
  foo = 1, ! A comment line that will be ignored. 
  bar = 2
 /
 &sample_b
  foo = 3, 
  bar = 4 
 / ! Another comment line
--nl_sample--

>>> import namelist
>>> nlc = namelist.read_namelists('nl_sample')
>>> nlc
[{'foo': 1, 'bar': 2}, {'foo': 3, 'bar': 4}]
>>> nlc[0].name
'sample_a'
>>> nlc.first_match('sample_a').name
'sample_a'
>>> print nlc
 &sample_a
  foo = 1, bar = 2
 /
 &sample_b
  foo = 3, bar = 4
 /

>>> for nl in nlc.matches('sample_a'): print nl
 &sample_a
  foo = 1, bar = 2
 /

"""
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
import re

__author__ = ('Tim Giese (TJG) - <giese@biomaps.rutgers.edu>\n'
              'Brian K. Radak (BKR) - <radakb@biomaps.rutgers.edu>')

__all__ = ['NamelistCollection','Namelist','separate_namelists',
           'read_non_namelists','read_namelists']


class NamelistCollection(list):
    """
    A list of Namelist objects. Elements can be looked up by matching 
    against their (possibly non-unique) "name" attribute. 
    """
    def __init__(self, *nls):
        list.__init__(self)
        self.extend(nls)

    def __str__(self):
        return ''.join([str(nl) for nl in self])

    def append(self, item):
        if not isinstance(item,Namelist):
            raise TypeError('NamelistCollections must contain Namelists!')
        list.append(self,item)

    def extend(self, items):
        if isinstance(items,list) or isinstance(items,tuple):
            for item in items:
                self.append(item)

    def matches(self, *names):
        """
        Return a generator of all namelists whose name matches any of the
        arguments.
        """
        for nl in self:
            for name in names:
                if nl.name == name:
                    yield nl

    def does_not_match(self, *names):
        """
        Return a generator of all namelists whose name does not match any
        of the arguments.
        """
        for nl in self:
            is_a_match = False
            for name in names:
                if nl.name == name:
                    is_a_match = True
                    break
            if not is_a_match:
                yield nl

    def first_match(self, *names):
        """
        Return the first result from matches(). If there are no matches,
        return None.
        """
        try:
            return next(self.matches(*names))
        except StopIteration:
            # This happens if matches() returns an empty generator.
            return None

    def first_non_match(self, *names):
        """
        Return the first result from does_not_match(). If there are no
        matches, return None.
        """
        try:
            return next(self.does_not_match(*names))
        except StopIteration:
            # This happens if does_not_match() returns an empty generator.
            return None


class Namelist(OrderedDict):
    """
    A dict-derived class containing the names and values of a Fortran 
    namelist. Names (keys) are stored as strings while values are cast as
    ints and floats when appropriate (Note: this does not occur for lists
    of ints and/or floats).

    Attributes
    ----------
    name : string, optional
        The name of the namelist.
    line_prefix : string, optional
        Begin all lines with name-value subsequences with this.
    name_value_separator : string, optional 
        Separate names and values in name-value subsequences with this.
    value_separator : string, optional
        Separate name-value subsequences with this.
    max_namevalues_per_line : int, optional
        Start a new line if the number of name-value subsequences on a 
        given line exceeds this.
    max_chars_per_line : int, optional
        Start a new line if the number of characters on a given line 
        exceeds this.
    """
    def __init__(self, name=None, line_prefix = ' ', 
                 name_value_separator = ' = ', value_separator = ', ',
                 max_namevalues_per_line = 72, max_chars_per_line = 72, *args,
                 **kwargs):
        self.name = name
        self.line_prefix = str(line_prefix)
        self.name_value_separator = str(name_value_separator)
        self.value_separator = str(value_separator)
        self.max_namevalues_per_line = int(max_namevalues_per_line)
        self.max_chars_per_line = int(max_chars_per_line)
        OrderedDict.__init__(self,*args,**kwargs)

    def __str__(self):
        txt = ' &%s\n'%self.name
        txt_buf = ' %s'%self.line_prefix
        nvalues = 0
        for name,value in self.iteritems():
            to_add = '%s%s%s%s'%(name,self.name_value_separator,value,
                                 self.value_separator)
            nvalues += 1
            if (len(txt_buf+to_add) >= self.max_chars_per_line or
                nvalues >= self.max_namevalues_per_line):
                txt += txt_buf + '\n'
                txt_buf = ' %s%s'%(self.line_prefix,to_add)
                nvalues = 0
            else:
                txt_buf += to_add
        txt += txt_buf.rstrip(self.value_separator) 
        txt += '\n /\n'
        return txt


def separate_namelists(filename):
    """
    Return a tuple of namelists (as a NamelistCollection) and non-namelists 
    (as a list of lines) from a file.
    """
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r'^(.*?)!.*$',r'\1',fileLines[i])
    bigLine1 = ''
    bigLine2 = ''
    for line in fileLines:
        bigLine1 += line.strip() + ' '
        bigLine2 += line
    return _namelists_from_string(bigLine1),_non_namelists_from_string(bigLine2)

def read_namelists(filename):
    """Return a NamelistCollection of all namelists contained in a file."""
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r'^(.*?)!.*$',r'\1',fileLines[i])
    bigLine = ''
    for line in fileLines:
        bigLine = bigLine + line.strip() + ' '
    return _namelists_from_string(bigLine)

def read_non_namelists(filename):
    """Return a list of lines not containing namelists in a file."""
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r'^(.*?)!.*$',r'\1',fileLines[i])
    bigLine = ''
    for line in fileLines:
        bigLine = bigLine + line
    return _non_namelists_from_string(bigLine)

def _as_common_type(string):
    try:
        return int(string)
    except ValueError:
        try:
            return float(string)
        except ValueError:
            return str(string)

def _namelists_from_string(string):
    # Return a NamelistCollection of all the namelists contained in a string.
    nlObjs = NamelistCollection()
    nls = re.findall(r'(&[a-zA-Z0-9_]+.*?[^\\]\/)',string)
    for nl in nls:
        nlName = None
        nlStr = None
        result = re.match(r'&([a-zA-Z0-9_]+)(.*)[^\\]\/',nl)
        if result is not None:
            nlName = result.group(1).strip()
            nlStr  = result.group(2).strip()
        keyvals = re.findall(r'(.*?)=([^=]+)',nlStr)
        for i in range(len(keyvals)):
            k,v = keyvals[i]
            if len(k) == 0:
                pk,pv = keyvals[i-1]
                cols = pv.split()
                k = cols.pop()
                pv = " ".join(cols)
                keyvals[i-1] = pk,pv
                keyvals[i] = k,v
        for i in range(len(keyvals)):
            k,v = keyvals[i]
            v = re.sub(r"\,$",'',v)
            keyvals[i] = k,v
        nlMap = OrderedDict()
        for k,v in keyvals:
            nlMap[k.strip()] = _as_common_type(v.strip())
        nlObjs.append(Namelist(nlName,**nlMap))
    return nlObjs

def _non_namelists_from_string(string):
    # Return all those parts of a string that are not a namelist.
    nls = re.findall(r'(&[a-zA-Z0-9_]+.*?[^\\]\/)',string,re.S)
    for nl in nls:
        string = string.replace(nl,'')
    lines = re.findall(r'(.*)\n*',string)
    newlines = []
    for line in lines:
        if len(line.strip()):
            newlines.append(line)
    return newlines
