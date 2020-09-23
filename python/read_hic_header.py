from __future__ import print_function
#Reads the genome name from the hic header
#Can take in a .hic file or URL that points to .hic file
import sys
import struct
import requests
import io

from straw import straw_module

if (len(sys.argv) != 2 and len(sys.argv) != 3):
    sys.stderr.write('Usage: '+ sys.argv[0]+' <hic file or URL> [verbose]\n')
    sys.exit(1)
verbose=0

if (len(sys.argv) == 3):
    verbose=1

_=straw_module.read_metadata(sys.argv[1],verbose=verbose)
