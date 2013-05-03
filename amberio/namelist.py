"""
FILE: namelist.py - Python routines for parsing files with Fortran namelists

DESCRIPTION:

AUTHOR: Tim Giese (TJG) - <giese@biomaps.rutgers.edu>
        Brian K. Radak (BKR) - <radakb@biomaps.rutgers.edu>
"""
import re

__all__ = ['NamelistCollection','Namelist','separate_namelists',
           'read_non_namelists','read_namelists']


class NamelistCollection(list):
    """
    A NamelistCollection is a list of Namelist objects. Namelist objects can be
    looked up by matching against their (possibly non-unique) names. 
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
        """Return all namelists whose name matches any of the arguments.
        """
        for nl in self:
            for name in names:
                if nl.name == name:
                    yield nl

    def does_not_match(self, *names):
        """Return all namelists whose name does not match any of the arguments.
        """
        for nl in self:
            is_a_match = False
            for name in names:
                if nl.name == name:
                    is_a_match = True
                    break
            if not is_a_match:
                yield nl

    def first_match(self, name):
        try:
            return next(self.matches(name))
        except StopIteration:
            # This happens if matches returns an empty generator.
            return None


class Namelist(dict):
    """
    A Namelist is a dictionary derived class that contains the names and values 
    of a Fortran namelist. Names (keys) are stored as strings while values are 
    cast as ints and floats when appropriate (Note: this does not occur for 
    lists of ints and/or floats).

    A Fortran namelist is comprised of:
    - an ampersand followed by the namelist name
    - zero or more name-value subsequences separated by a value separator
    - a terminating slash

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
        Start a new line if the number of name-value subsequences on a given
        line exceeds this.
    max_chars_per_line : int, optional
        Start a new line if the number of characters on a given line exceeds
        this.
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
        dict.__init__(self,*args,**kwargs)

    def __str__(self):
        txt = ' &%s\n'%self.name
        txt_buf = ' %s'%self.line_prefix
        nvalues = 0
        for name,value in self.iteritems():
            txt_buf += '%s%s%s%s'%(name,self.name_value_separator,value,
                                   self.value_separator)
            nvalues += 1
            if (len(txt_buf) >= self.max_chars_per_line or
                nvalues >= self.max_namevalues_per_line):
                txt += txt_buf + '\n'
                txt_buf = ' %s'%self.line_prefix
                nvalues = 0
        txt += txt_buf.rstrip(self.value_separator) 
        txt += '\n /\n'
        return txt


def ValueIsInt(e):
    isInt = False
    isFloat = False
    isString = False
    try:
        a = int(e)
        isInt = True
    except ValueError:
        try:
            a = float(e)
            isFloat = True
        except ValueError:
            isString = True
    return isInt

def ValueIsFloat(e):
    isInt = False
    isFloat = False
    isString = False
    try:
        a = int(e)
        isInt = True
    except ValueError:
        try:
            a = float(e)
            isFloat = True
        except ValueError:
            isString = True
    return isFloat

def ValueIsString(e):
    isInt = False
    isFloat = False
    isString = False
    try:
        a = int(e)
        isInt = True
    except ValueError:
        try:
            a = float(e)
            isFloat = True
        except ValueError:
            isString = True
    return isString

def AllElementsAreInt(eles):
    isInt = False
    for e in eles:
        isInt = ValueIsInt(e)
        if not isInt:
            break
    return isInt

def AllElementsAreFloat(eles):
    isFloat = False
    for e in eles:
        isFloat = ValueIsFloat(e)
        if not isFloat:
            break
    return isFloat

def ConvertToInts(eles):
    for i in range(len(eles)):
        eles[i] = int(eles[i])
    return eles

def ConvertToFloats(eles):
    for i in range(len(eles)):
        eles[i] = float(eles[i])
    return eles

def ConvertToCommonType(e):
    if ValueIsInt(e): 
        return int(e)
    elif ValueIsFloat(e):
        return float(e)
    else:
        return str(e)

def IsLikeAList(v):
   """Return True if v is a non-string sequence and is iterable. Note that
   not all objects with getitem() have the iterable attribute
   Taken from http://stackoverflow.com/questions/836387/how-can-i-tell-if-a-python-variable-is-a-string-or-a-list
   """
   if hasattr(v, '__iter__') and not isinstance(v, basestring):
       return True
   else:
       #This will happen for most atomic types like numbers and strings
       return False

def ConvertToString(v):
    line = ""
    if IsLikeAList(v):
        for item in iter(v):
            if ValueIsString(item):
                line += "'" + str(item) + "' "
            else:
                line += str(item) + " "
    else:
        line = str(v)
    return line

def separate_namelists(filename):
    """
    Return a tuple of namelists (as a NamelistCollection) and non-namelists 
    (as a list of lines) from a file.
    """
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r"^(.*?)!.*$",r"\1",fileLines[i])
    bigLine1 = ""
    bigLine2 = ""
    for line in fileLines:
        bigLine1 += line.strip() + " "
        bigLine2 += line
    return _namelists_from_string(bigLine1),_non_namelists_from_string(bigLine2)

def read_namelists(filename):
    """Return a NamelistCollection of all namelists contained in a file."""
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r"^(.*?)!.*$",r"\1",fileLines[i])
    bigLine = ""
    for line in fileLines:
        bigLine = bigLine + line.strip() + " "
    return _namelists_from_string(bigLine)

def read_non_namelists(filename):
    """Return a list of lines not containing namelists in a file."""
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r"^(.*?)!.*$",r"\1",fileLines[i])
    bigLine = ""
    for line in fileLines:
        bigLine = bigLine + line
    return _non_namelists_from_string(bigLine)

def _namelists_from_string(string):
    # Return a NamelistCollection of all the namelists contained in a string.
    nlObjs = NamelistCollection()
    nls = re.findall(r"(&[a-zA-Z]+.*?[^\\]\/)",string)
    for nl in nls:
        nlName = None
        nlStr  = None
        result = re.match(r"&([a-zA-Z]+)(.*)[^\\]\/",nl)
        if result is not None:
            nlName = result.group(1).strip()
            nlStr  = result.group(2).strip()
        keyvals = re.findall(r"(.*?)=([^=]+)",nlStr)
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
            v = re.sub(r"\,$","",v)
            keyvals[i] = k,v
        nlMap = {}
        for key,val in keyvals:
            nlMap[key.strip()] = ConvertToCommonType(val.strip())
        nlObjs.append(Namelist(nlName,**nlMap))
    return nlObjs

def _non_namelists_from_string(string):
    # Return all those parts of a string that are not a namelist.
    nls = re.findall(r"(&[a-zA-Z]+.*?[^\\]\/)",string,re.S)
    for nl in nls:
        string = string.replace(nl,"")
    lines = re.findall(r"(.*)\n*",string)
    newlines = []
    for line in lines:
        if len(line.strip()):
            newlines.append(line)
    return newlines

