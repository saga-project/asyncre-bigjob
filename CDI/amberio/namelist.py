#! /usr/bin/env python
################################################################################
#                                                                              
# FILE: namelist.py - Python routines for parsing files with Fortran namelists
#
# DESCRIPTION: 
#
# AUTHOR: Tim Giese (TJG) - <giese@biomaps.rutgers.edu>
#         Brian K. Radak (BKR) - <radakb@biomaps.rutgers.edu>
#
################################################################################
import re
class Namelist(dict):
    """
    A dictionary derived class that contains the keys and values of a Fortran
    namelist. The keys are generally strings while the values (which are comma
    delimited in Fortran files) are lists, even if only one value is present.
    """
    def __init__(self,name,*args,**kwargs):
        self.name = name
        super(Namelist, self).__init__(*args,**kwargs)

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

def ConvertToCommonType(eles):
    if AllElementsAreInt(eles):
        eles = ConvertToInts(eles)
    elif AllElementsAreFloat(eles):
        eles = ConvertToFloats(eles)
    return eles

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

    nlObjs = []
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
            key = key.strip()
            val = val.strip()
            values = val.split(",")
            values = ConvertToCommonType(values)
            nlMap[key] = values

        nlObjs.append(Namelist(nlName,nlMap))

    return nlObjs
