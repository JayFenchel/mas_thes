#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = 'jayf'

import numpy as np

# from cvxopt.modeling import op
#
# lp = op()
# lp.fromfile('QPTEST.QPS')
#
# print(lp.solve())

def fromfile(filename):

    '''
    Reads LP from file 'filename' assuming it is a fixed format
    ascii MPS file.

    Does not include serious error checking.

    MPS features that are not allowed: comments preceded by
    dollar signs, linear combinations of rows, multiple righthand
    sides, ranges columns or bounds columns.
    '''

    _inequalities = []
    _equalities = []
    objective = []
    name = ''

    f = open(filename, 'r')

    s = f.readline()
    while s[:4] != 'NAME':
        s = f.readline()
        if not s:
            raise SyntaxError("EOF reached before 'NAME' section "\
                "was found")
    name = s[14:22].strip()

    s = f.readline()
    while s[:4] != 'ROWS':
        if not s:
            raise SyntaxError("EOF reached before 'ROWS' section "\
                "was found")
        s = f.readline()
    s = f.readline()

    # ROWS section
    functions = dict()   # {MPS row label: affine function}
    rowtypes = dict()    # {MPS row label: 'E', 'G' or 'L'}
    foundobj = False     # first occurrence of 'N' counts
    while s[:7] != 'COLUMNS':
        if not s: raise SyntaxError("file has no 'COLUMNS' section")
        if len(s.strip()) == 0 or s[0] == '*':
            pass
        elif s[1:3].strip() in ['E', 'L', 'G']:
            rowlabel = s[4:12].strip()
            functions[rowlabel] = np.array([[]])
            rowtypes[rowlabel] = s[1:3].strip()
        elif s[1:3].strip() == 'N':
            rowlabel = s[4:12].strip()
            if not foundobj:
                functions[rowlabel] = objective
                foundobj = True
        else:
            raise ValueError("unknown row type '%s'" %s[1:3].strip())
        s = f.readline()
    s = f.readline()


    # COLUMNS section
    variables = dict()   # {MPS column label: variable}
    while s[:3] != 'RHS':
        if not s:
            raise SyntaxError("EOF reached before 'RHS' section "\
                "was found")
        if len(s.strip()) == 0 or s[0] == '*':
            pass
        else:
            if s[4:12].strip(): collabel = s[4:12].strip()
            if collabel not in variables:
                variables[collabel] = variable(1, collabel)
            v = variables[collabel]
            rowlabel = s[14:22].strip()
            if rowlabel not in functions:
                raise KeyError("no row label '%s'" %rowlabel)
            functions[rowlabel]._linear._coeff[v] = \
                matrix(float(s[24:36]), tc='d')
            rowlabel = s[39:47].strip()
            if rowlabel:
                if rowlabel not in functions:
                    raise KeyError("no row label '%s'" %rowlabel)
                functions[rowlabel]._linear._coeff[v] =  \
                    matrix(float(s[49:61]), tc='d')
        s = f.readline()
    s = f.readline()

    f.close()
    return functions



    # # RHS section
    # # The RHS section may contain multiple right hand sides,
    # # identified with different labels in s[4:12].
    # # We read in only one rhs, the one with the first rhs label
    # # encountered.
    # rhslabel = None
    # while s[:6] != 'RANGES' and s[:6] != 'BOUNDS' and \
    #     s[:6] != 'ENDATA':
    #     if not s: raise SyntaxError( \
    #          "EOF reached before 'ENDATA' was found")
    #     if len(s.strip()) == 0 or s[0] == '*':
    #         pass
    #     else:
    #         if None != rhslabel != s[4:12].strip():
    #             # skip if rhslabel is different from 1st rhs label
    #             # encountered
    #             pass
    #         else:
    #             if rhslabel is None: rhslabel = s[4:12].strip()
    #             rowlabel = s[14:22].strip()
    #             if rowlabel not in functions:
    #                 raise KeyError("no row label '%s'" %rowlabel)
    #             functions[rowlabel]._constant = \
    #                 matrix(-float(s[24:36]), tc='d')
    #             rowlabel = s[39:47].strip()
    #             if rowlabel:
    #                 if rowlabel not in functions:
    #                     raise KeyError("no row label '%s'" \
    #                         %rowlabel)
    #                 functions[rowlabel]._constant = \
    #                     matrix(-float(s[49:61]), tc='d')
    #     s = f.readline()
    #
    #
    # # RANGES section
    # # The RANGES section may contain multiple range vectors,
    # # identified with different labels in s[4:12].
    # # We read in only one vector, the one with the first range label
    # # encountered.
    # ranges = dict()
    # for l in iter(rowtypes.keys()):
    #     ranges[l] = None   # {rowlabel: range value}
    # rangeslabel = None
    # if s[:6] == 'RANGES':
    #     s = f.readline()
    #     while s[:6] != 'BOUNDS' and s[:6] != 'ENDATA':
    #         if not s: raise SyntaxError( \
    #             "EOF reached before 'ENDATA' was found")
    #         if len(s.strip()) == 0 or s[0] == '*':
    #             pass
    #         else:
    #             if None != rangeslabel != s[4:12].strip():
    #                 pass
    #             else:
    #                 if rangeslabel == None:
    #                     rangeslabel = s[4:12].strip()
    #                 rowlabel = s[14:22].strip()
    #                 if rowlabel not in rowtypes:
    #                     raise KeyError("no row label '%s'"%rowlabel)
    #                 ranges[rowlabel] = float(s[24:36])
    #                 rowlabel = s[39:47].strip()
    #                 if rowlabel != '':
    #                     if rowlabel not in functions:
    #                         raise KeyError("no row label '%s'" \
    #                             %rowlabel)
    #                     ranges[rowlabel] =  float(s[49:61])
    #         s = f.readline()
    #
    #
    # # BOUNDS section
    # # The BOUNDS section may contain bounds vectors, identified
    # # with different labels in s[4:12].
    # # We read in only one bounds vector, the one with the first
    # # label encountered.
    # boundslabel = None
    # bounds = dict()
    # for v in iter(variables.keys()):
    #     bounds[v] = [0.0, None] #{column label: [l.bound, u. bound]}
    # if s[:6] == 'BOUNDS':
    #     s = f.readline()
    #     while s[:6] != 'ENDATA':
    #         if not s: raise SyntaxError( \
    #             "EOF reached before 'ENDATA' was found")
    #         if len(s.strip()) == 0 or s[0] == '*':
    #             pass
    #         else:
    #             if None != boundslabel != s[4:12].strip():
    #                 pass
    #             else:
    #                 if boundslabel is None:
    #                     boundslabel = s[4:12].strip()
    #                 collabel = s[14:22].strip()
    #                 if collabel not in variables:
    #                     raise ValueError('unknown column label ' \
    #                         + "'%s'" %collabel)
    #                 if s[1:3].strip() == 'LO':
    #                     if bounds[collabel][0] != 0.0:
    #                         raise ValueError("repeated lower "\
    #                             "bound for variable '%s'" %collabel)
    #                     bounds[collabel][0] = float(s[24:36])
    #                 elif s[1:3].strip() == 'UP':
    #                     if bounds[collabel][1] != None:
    #                         raise ValueError("repeated upper "\
    #                             "bound for variable '%s'" %collabel)
    #                     bounds[collabel][1] = float(s[24:36])
    #                 elif s[1:3].strip() == 'FX':
    #                     if bounds[collabel] != [0, None]:
    #                         raise ValueError("repeated bounds "\
    #                             "for variable '%s'" %collabel)
    #                     bounds[collabel][0] = float(s[24:36])
    #                     bounds[collabel][1] = float(s[24:36])
    #                 elif s[1:3].strip() == 'FR':
    #                     if bounds[collabel] != [0, None]:
    #                         raise ValueError("repeated bounds "\
    #                             "for variable '%s'" %collabel)
    #                     bounds[collabel][0] = None
    #                     bounds[collabel][1] = None
    #                 elif s[1:3].strip() == 'MI':
    #                     if bounds[collabel][0] != 0.0:
    #                         raise ValueError("repeated lower " \
    #                             "bound for variable '%s'" %collabel)
    #                     bounds[collabel][0] = None
    #                 elif s[1:3].strip() == 'PL':
    #                     if bounds[collabel][1] != None:
    #                         raise ValueError("repeated upper " \
    #                             "bound for variable '%s'" %collabel)
    #                 else:
    #                     raise ValueError("unknown bound type '%s'"\
    #                         %s[1:3].strip())
    #         s = f.readline()
    #
    # for l, type in iter(rowtypes.items()):
    #
    #     if type == 'L':
    #         c = functions[l] <= 0.0
    #         c.name = l
    #         _inequalities += [c]
    #         if ranges[l] != None:
    #             c = functions[l] >= -abs(ranges[l])
    #             c.name = l + '_lb'
    #             _inequalities += [c]
    #     if type == 'G':
    #         c = functions[l] >= 0.0
    #         c.name = l
    #         _inequalities += [c]
    #         if ranges[l] != None:
    #             c = functions[l] <= abs(ranges[l])
    #             c.name = l + '_ub'
    #             _inequalities += [c]
    #     if type == 'E':
    #         if ranges[l] is None or ranges[l] == 0.0:
    #             c = functions[l] == 0.0
    #             c.name = l
    #             _equalities += [c]
    #         elif ranges[l] > 0.0:
    #             c = functions[l] >= 0.0
    #             c.name = l + '_lb'
    #             _inequalities += [c]
    #             c = functions[l] <= ranges[l]
    #             c.name = l + '_ub'
    #             _inequalities += [c]
    #         else:
    #             c = functions[l] <= 0.0
    #             c.name = l + '_ub'
    #             _inequalities += [c]
    #             c = functions[l] >= ranges[l]
    #             c.name = l + '_lb'
    #             _inequalities += [c]
    #
    # for l,bnds in iter(bounds.items()):
    #     v = variables[l]
    #     if None != bnds[0] != bnds[1]:
    #         c = v >= bnds[0]
    #         _inequalities += [c]
    #     if bnds[0] != bnds[1] != None:
    #         c = v  <= bnds[1]
    #         _inequalities += [c]
    #     if None != bnds[0] == bnds[1]:
    #         c = v == bnds[0]
    #         _equalities += [c]
    #
    # # Eliminate constraints with no variables
    # for c in _inequalities + _equalities:
    #     if len(c._f._linear._coeff) == 0:
    #         if c.type() == '=' and c._f._constant[0] != 0.0:
    #             raise ValueError("equality constraint '%s' "\
    #                "has no variables and a nonzero righthand side"\
    #                %c.name)
    #         elif c.type() == '<' and c._f._constant[0] > 0.0:
    #             raise ValueError("inequality constraint '%s' "\
    #                "has no variables and a negative righthand side"\
    #                %c.name)
    #         else:
    #             print("removing redundant constraint '%s'" %c.name)
    #             if c.type() == '<': _inequalities.remove(c)
    #             if c.type() == '=': _equalities.remove(c)
    #
    #
    # _variables = dict()
    # for v in objective._linear._coeff.keys():
    #     _variables[v] = {'o': True, 'i': [], 'e': []}
    # for c in _inequalities:
    #     for v in c._f._linear._coeff.keys():
    #         if v in _variables:
    #             _variables[v]['i'] += [c]
    #         else:
    #             _variables[v] = {'o': False, 'i': [c], 'e': []}
    # for c in _equalities:
    #     for v in c._f._linear._coeff.keys():
    #         if v in _variables:
    #             _variables[v]['e'] += [c]
    #         else:
    #             _variables[v] = {'o': False, 'i': [], 'e': [c]}
    #
    # status = None
    #
    # f.close()