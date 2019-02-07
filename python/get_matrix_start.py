from __future__ import print_function
#Reads the genome name from the hic header
#Can take in a .hic file or URL that points to .hic file
import sys
import struct
import io

def readcstr(f):
    buf = ""
    while True:
        b = f.read(1)
        b = b.decode('utf-8', 'backslashreplace')
        if b is None or b == '\0':
            return str(buf)
        else:
            buf = buf + b

if (len(sys.argv) != 2):
  sys.stderr.write('Usage: '+ sys.argv[0]+' <hic file or URL>\n')
  sys.exit(1)

infile = sys.argv[1]
req=open(infile, 'rb')
magic_string = struct.unpack('<3s', req.read(3))[0]
req.read(1)
if (magic_string != b"HIC"):
    print('This does not appear to be a HiC file magic string is incorrect')
    sys.exit(1)
version = struct.unpack('<i',req.read(4))[0]
masterindex = struct.unpack('<q',req.read(8))[0]
genome = ""
c=req.read(1).decode("utf-8") 
while (c != '\0'):
    genome += c
    c=req.read(1).decode("utf-8") 

nattributes = struct.unpack('<i',req.read(4))[0]
for x in range(0, nattributes):
  key = readcstr(req)
  value = readcstr(req)
nChrs = struct.unpack('<i',req.read(4))[0]
for x in range(0, nChrs):
    name = readcstr(req)
    length = struct.unpack('<i',req.read(4))[0]
nBpRes = struct.unpack('<i',req.read(4))[0]
for x in range(0, nBpRes):
    res = struct.unpack('<i',req.read(4))[0]
nFrag = struct.unpack('<i',req.read(4))[0]
for x in range(0, nFrag):
    res = struct.unpack('<i',req.read(4))[0]
if (nFrag > 0):
    for x in range(0, nChrs):
        nSites = struct.unpack('<i',req.read(4))[0]
        for y in range(0, nSites):
            site = struct.unpack('<i',req.read(4))[0]
print("{}".format(req.tell()))
