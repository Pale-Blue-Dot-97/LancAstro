##############################################################
# Set of definitions to read columns from a whitespace 
# delimited file with '#' comments (e.g. Sextractor file)
#
# By Jim Geach
# Modified by Philip Best (August 2007)
# Modified by David Sobral (November 2007)
##############################################################


from sys import *
from string import *
from numpy import *



################################################################
# Definition "count": counts number of lines in a file
################################################################

def count(f):
    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()
    j=0
    return len(alldata_lines)


################################################################
# Definition "skip": counts commented lines at start of file
################################################################

def skip(f):
    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()
    n=0
    while alldata_lines[n][0]=='#':
        n+=1
    return n	
	

################################################################
# Definition "read_col": reads a column "c" from a file,       #
# skipping the commented lines at start of file                #
################################################################
	
def read_col(c, f):
    L = count(f)
    sk = skip(f)
    L-=sk
	
    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()
	
    D = array(range(L),float32)
    for i in range(L):
        line = alldata_lines[i+sk]
        sp = line.split()
        D[i] = float(sp[c-1])
    return D

################################################################
# Definition "mread_col": reads multiple column "cs" from a 
# file, skipping the commented lines at start of file
################################################################

def mread_col(cs, f):
    L = count(f)
    sk = skip(f)
    L-=sk
	
    file_object = open(f)
    alldata_lines = file_object.readlines()
    file_object.close()

    O = []

    for j in cs:
        D = array(range(L),float32)
        for i in range(L):
            line = alldata_lines[i+sk]
            sp = line.split()
            D[i] = float(sp[j-1])
        O.append(D)
    return O


