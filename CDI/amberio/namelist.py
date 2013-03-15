"""
FILE: namelist.py - Python routines for parsing files with Fortran namelists

DESCRIPTION:

AUTHOR: Tim Giese (TJG) - <giese@biomaps.rutgers.edu>
        Brian K. Radak (BKR) - <radakb@biomaps.rutgers.edu>
"""
__all__ = ['NamelistCollection','Namelist','ReadEverythingExceptNamelists',
           'ReadNamelists']

import re

class NamelistCollection(list):
    """A list of namelists with easy access to namelist keys and values.
    """
    def __init__(self,*nls):
        list.__init__(self)
        for nl in nls: self.extend(nl)

    def append(self, item):
        if not isinstance(item, Namelist):
            raise TypeError('NamelistCollections must contain Namelists!')
        list.append(self, item)

    def extend(self, items):
        if hasattr(items, '__iter__'):
            for item in items:
                self.append(item)

    def GetAllMatches(self, name):
        return (nl for nl in self if nl.name == name)

    def GetFirstMatch(self, name):
        return next(self.GetAllMatches(name))

class Namelist(dict):
    """
    A dictionary derived class that contains the keys and values of a Fortran
    namelist. Keys are stored as strings while values are cast as ints and 
    floats when appropriate (Note: this does not occur for lists of ints and/or 
    floats).
    """
    def __init__(self,name,*args,**kwargs):
        self.name = name
        super(Namelist, self).__init__(*args,**kwargs)

    def SprintFortran(self):
        line = "&" + self.name.strip() + "\n"
        for k,v in self.iteritems():
            line += k + " = " + ConvertToString(v) + "\n"
        line += "/\n"
        return line

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

# Replaced this when it was decided to not return lists anymore
# def ConvertToCommonType(eles):
#     if AllElementsAreInt(eles):
#         eles = ConvertToInts(eles)
#     elif AllElementsAreFloat(eles):
#         eles = ConvertToFloats(eles)
#     return eles

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

def ReadEverythingExceptNamelists(filename):
    """
    Read a file containing Fortran namelists. Return a list of the lines that do
    NOT contain anything inside a namelist.
    """
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r"^(.*?)!.*$",r"\1",fileLines[i])
    bigLine = ""
    for line in fileLines:
        bigLine = bigLine + line

    nls = re.findall(r"(&[a-zA-Z]+.*?[^\\]\/)",bigLine,re.S)
    for nl in nls:
        bigLine = bigLine.replace(nl,"")
    lines = re.findall(r"(.*)\n*",bigLine)
    newlines = []
    for line in lines:
        if len(line.strip()):
            newlines.append(line)
    return newlines

def ReadNamelists(filename):
    """
    Read a file containing Fortran namelists. Return a list of (dict-derived)
    Namelist objects for each namelist found.
    """
    fileLines = open(filename,'r').readlines()
    for i in range(len(fileLines)):
        fileLines[i] = re.sub(r"^(.*?)!.*$",r"\1",fileLines[i])
    bigLine = ""
    for line in fileLines:
        bigLine = bigLine + line.strip() + " "

    nlObjs = NamelistCollection()
    nls = re.findall(r"(&[a-zA-Z]+.*?[^\\]\/)",bigLine)
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

        nlObjs.append(Namelist(nlName,nlMap))

    return nlObjs
