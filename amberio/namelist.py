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
>>> nlc = namelist.Namelist.from_file('nl_sample')
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
import re
try:
    from collections import OrderedDict
except ImportError:
    # for earlier than python2.5?
    from ordereddict import OrderedDict

__author__ = ('Tim Giese (TJG) - <giese@biomaps.rutgers.edu>\n'
              'Brian K. Radak (BKR) - <radakb@biomaps.rutgers.edu>')

__all__ = ['NamelistCollection','Namelist']


class NamelistCollection(list):
    """
    A list of Namelist objects. Elements can be looked up by matching 
    against their (possibly non-unique) "name" attribute. 
    """
    def __init__(self, *nls):
        list.__init__(self)
        self.extend(nls)

    @classmethod
    def separate_nls(cls, filename):
        """
        Return a NamelistCollection of all namelists contained in a file AND a
        list of all remaining lines not containing namelists.
        """
        lines = open(filename,'r').readlines()
        for i in range(len(lines)):
            lines[i] = re.sub(r'^(.*?)!.*$',r'\1',lines[i])
        bigline1 = ''
        bigline2 = ''
        for line in lines:
            bigline1 += line.strip() + ' '
            bigline2 += line
        return cls.from_str(bigline1),cls.non_nls_from_str(cls,bigline2)

    @classmethod
    def from_file(cls, filename):
        """Return a NamelistCollection of all namelists contained in a file."""
        lines = open(filename,'r').readlines()
        for i in range(len(lines)):
            lines[i] = re.sub(r'^(.*?)!.*$',r'\1',lines[i])
        bigline = ''
        for line in lines:
            bigline = bigline + line.strip() + ' '
        return cls.from_str(bigline) 
    
    @classmethod
    def from_str(cls, string):
        """Return a NamelistCollection of all namelists in a string."""
        inst = cls()
        nls = re.findall(r'(&[a-zA-Z0-9_]+.*?[^\\]\/)',string)
        for nl in nls:
            nl_name = None
            nl_str = None
            result = re.match(r'&([a-zA-Z0-9_]+)(.*)[^\\]\/',nl)
            if result is not None:
                nl_name = result.group(1).strip()
                nl_str  = result.group(2).strip()
            keyvals = re.findall(r'(.*?)=([^=]+)',nl_str)
            for i in range(len(keyvals)):
                k,v = keyvals[i]
                if len(k) == 0:
                    pk,pv = keyvals[i-1]
                    cols = pv.split()
                    k = cols.pop()
                    pv = ' '.join(cols)
                    keyvals[i-1] = pk,pv
                    keyvals[i] = k,v
            for i in range(len(keyvals)):
                k,v = keyvals[i]
                v = re.sub(r"\,$",'',v)
                keyvals[i] = k,v
            nl_map = Namelist(nl_name)
            for k,v in keyvals:
                v = v.strip()
                try:
                    v = int(v)
                except ValueError:
                    try:
                        v = float(v)
                    except ValueError:
                        pass
                nl_map[k.strip()] = v
            inst.append(nl_map)
        return inst

    @staticmethod
    def non_nls_from_file(self, filename):
        """Return a list of lines from a file that do not contain namelists."""
        lines = open(filename,'r').readlines()
        for i in range(len(lines)):
            lines[i] = re.sub(r'^(.*?)!.*$',r'\1',lines[i])
        bigline = ''
        for line in lines:
            bigline = bigline + line
        return self.non_nls_from_str(self,bigline)

    @staticmethod
    def non_nls_from_str(self, string):
        """Return all those parts of a string that are not a namelist."""
        nls = re.findall(r'(&[a-zA-Z0-9_]+.*?[^\\]\/)',string,re.S)
        for nl in nls:
            string = string.replace(nl,'')
        lines = re.findall(r'(.*)\n*',string)
        newlines = []
        for line in lines:
            if len(line.strip()):
                newlines.append(line)
        return newlines

    def __str__(self):
        return ''.join([str(nl) for nl in self])

    def append(self, item):
        if not isinstance(item,Namelist):
            raise TypeError('NamelistCollections must contain Namelists!')
        list.append(self,item)

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
        return None. This is similar to list.index() except that it uses the
        namelist name as the value.
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
        Begin all lines with name-value sub-sequences with this.
    name_value_separator : string, optional 
        Separate names and values in name-value sub-sequences with this.
    value_separator : string, optional
        Separate name-value subsequences with this.
    max_namevalues_per_line : int, optional
        Start a new line if the number of name-value sub-sequences on a 
        given line exceeds this.
    max_chars_per_line : int, optional
        Start a new line if the number of characters on a given line 
        exceeds this.
    """
    def __init__(self, name=None, line_prefix = ' ', 
                 name_value_separator = ' = ', value_separator = ', ',
                 max_namevalues_per_line = 72, max_chars_per_line = 72, *args,
                 **kwargs):
        if name is None:
            self.name = name
        else:
            self.name = str(name)
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
            if (len(txt_buf+to_add) >= self.max_chars_per_line
                or nvalues >= self.max_namevalues_per_line):
                txt += txt_buf + '\n'
                txt_buf = ' %s%s'%(self.line_prefix,to_add)
                nvalues = 0
            else:
                txt_buf += to_add
        txt += txt_buf.rstrip(self.value_separator) 
        txt += '\n /\n'
        return txt
