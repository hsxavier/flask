#! /usr/bin/env python

from contextlib import contextmanager
from itertools import imap, izip
from glob import iglob
from math import sqrt
from sys import exit, argv

# Colours:
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


@contextmanager
def multi_file_manager(files, mode='rt'):
    files = [open(file, mode) for file in files]
    yield files
    for file in files:
        file.close()

# generator function to read and yield each value from a file
def read_values(file):
    for line in file:
        if not line.startswith("#"):                # skip comments.  
            for value in imap(float, line.split()): # might only need 'int' here
                yield value

# enumerate multiple egual length iterables simultaneously as (i. n0, n1, ...)
def multi_enumerate(start, *iterables):
    # returns generator
    return ((n,)+t for n, t in enumerate(izip(*iterables), start))


# Beggining of the main code:
argv.pop(0)                   # First argument is the code itself.
meanfile = argv.pop(0)        # Second is the file output for the means.
devfile  = argv.pop(0)        # Third is the file output for the standard deviations.
#                               Remaining are the files that will be averaged over.
 

# Safeguard against wrong input:
#print ""
#print "Will write means to"+bcolors.FAIL, meanfile +bcolors.ENDC, "and std. dev. to"+bcolors.FAIL, devfile+bcolors.ENDC
#go = raw_input("Is this what you want (y/n)? ")
#if go!="y" and go!="Y":
#    exit(0)
#print ""

# Run code:
with multi_file_manager(argv) as datfiles:
    num_files = len(datfiles)
    if num_files < 2:
        print 'Less than 2 .dat files were found to process, terminating.'
        sys.exit(1)

    # determine number of rows and cols from first file
    temp = [line.split() for line in datfiles[0]]
    while (temp[0][0].startswith("#")==True): temp.pop(0)
    num_rows = len(temp)
    num_cols = len(temp[0])
    del temp  # no longer needed
    datfiles[0].seek(0)  # rewind first file
    print 'Found {} .dat files, each {} rows x {} cols'.format(num_files, num_rows, num_cols)
    print 'Will write their means to', meanfile, 'and their deviations to', devfile
    means  = []  # initialize
    sigmas = []  # standard deviations
    generators = [read_values(file) for file in datfiles]
    for j in xrange(num_rows):  # main loop
        for i in xrange(num_cols):
            values = map(next, generators)  # next cell value from each file
            mean = float(sum(values)) / num_files
            means.append(mean)
            means_diff_sq = imap(lambda value: (value-mean)**2, values)
            sigma = sqrt(sum(means_diff_sq) / (num_files-1)) # Corrected to give unbiased sample standard deviation.
            sigmas.append(sigma)

print 'Calculating means and deviations...'
with open(meanfile, 'wt') as averages, open(devfile, 'wt') as deviations:
    for i, mean, sigma in multi_enumerate(0, means, sigmas):
        averages.write('{:.8e}'.format(mean))
        deviations.write('{:.8e}'.format(sigma))
        if i % num_cols != num_cols-1:  # not last column?
             averages.write(' ')        # delimiter between values on line
             deviations.write(' ')  
        else:
            averages.write('\n')
            deviations.write('\n')       
print 'done.'
print

