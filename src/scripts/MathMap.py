#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script performs mathematical operations with Healpix maps.

HOW TO WRITE THE EXPRESSIONS:
-- It is recommended to separate all files, numbers and operators in the expression
   by whitespaces, although a few cases might work without them (e.g. powers and 
   parenthesis).
-- Parenthesis must be input with the scape character '\\', e.g. \\( 1 + 2 \\)'  
-- Multiplication is described by ' x ' (with whitespace) and not by '*' (wich is 
   a wildcard).
-- Powers are represented by '^', not '**' (for the same reason);

If the expression starts with '<FILE> = ', the result of the operation is saved 
to <FILE> and not shown.
Else, the result is shown and not saved anywhere.

USAGE:   MathMap.py <EXPRESSION> 
EXAMPLE: MathMap.py sim01-poisson1-f1z1.fits + contam01-f1z1.fits
EXAMPLE: MathMap.py sim01+contam01.fits = sim01-poisson1-f1z1.fits + contam01-f1z1.fits

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 11/jun/2018
"""

import healpy as hp
import numpy as np
import sys
import re
import parser
import matplotlib.pyplot as pl

# Docstring output:
if len(sys.argv) < 1 + 1: 
    print(__doc__)
    sys.exit(1)

# Load input expression:
# Check if there is demand for saving the result:
if sys.argv[1][-1]=='=':
    WillSave = True
    expression = ' '.join(sys.argv[2:])
    destiny = sys.argv[1][:-1]
elif sys.argv[2][0]=='=':
    WillSave = True
    destiny = sys.argv[1]
    if len(sys.argv[2])>1:
        expression = sys.argv[2][1:]+' '.join(sys.argv[3:])
    else:
        expression = ' '.join(sys.argv[3:])
else:
    WillSave = False
    expression = ' '.join(sys.argv[1:])


# Check if it is a float number:
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# Check if it is an operator:
op = ['+','-','x','/','^','(',')']
def is_operator(s):
    if s in op:
        return True
    return False

# Find maps, load them:
print ''
Map = []
MapFile = []
MapVar = []
MapCount = 0
term = re.split(' \+ | - | / | x |\^| \^|\^ | \^ |\(| \(|\( | \( |\)| \)|\) | \) | == ', expression)
for i in range(0,len(term)):
    if is_number(term[i])==False and is_operator(term[i])==False and len(term[i])>0:
        MapFile.append(term[i])
        MapVar.append('Map[{}]'.format(MapCount))
        Map.append(hp.read_map(MapFile[-1], verbose=False))
        MapCount = MapCount + 1

# Prepare expression for compiling:
expression = expression.replace(' x ','*')
expression = expression.replace('^','**')
FileLength = [len(mf) for mf in MapFile]
# Replace filename by variables, starting from the longest filename to avoid confusion:
zipped = zip(FileLength,MapFile,MapVar)
zipped.sort(key = lambda t: t[0], reverse=True)
for z in zipped:
    expression = expression.replace(z[1],z[2])

# Print a summary of the operations:
print 'Will do:', expression
print 'With:'
for i in range(0,MapCount):
    print 'Map[{}] = {}'.format(i,MapFile[i])
if WillSave==True:
    print 'Result will be saved to:', destiny
print ''

# Perform the operation:
code = parser.expr(expression).compile()
result = eval(code)

# Output:
if WillSave==True:
    hp.write_map(destiny,result)
else:
    hp.mollview(result)
    pl.show()
