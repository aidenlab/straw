from __future__ import print_function
#Reads the genome name from the hic header
#Can take in a .hic file or URL that points to .hic file
import sys
import struct
import requests
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
#            buf.append(b.decode('utf'))

if (len(sys.argv) != 2):
  sys.stderr.write('Usage: '+ sys.argv[0]+' <hic file or URL>\n')
  sys.exit(1)

infile = sys.argv[1]
magic_string = ""

if (infile.startswith("http")):
    # try URL first. 100K should be sufficient for header
    headers={'range' : 'bytes=0-100000', 'x-amz-meta-requester' : 'straw'}
    s = requests.Session()
    r=s.get(infile, headers=headers)
    if (r.status_code >=400):
        print("Error accessing " + infile) 
        print("HTTP status code " + str(r.status_code))
        sys.exit(1)
    req=io.BytesIO(r.content)        
    myrange=r.headers['content-range'].split('/')
    totalbytes=myrange[1]
else:
    req=open(infile, 'rb')

magic_string = struct.unpack('<3s', req.read(3))[0]
req.read(1)
if (magic_string != b"HIC"):
    print('This does not appear to be a HiC file magic string is incorrect')
    sys.exit(1)
version = struct.unpack('<i',req.read(4))[0]
print('HiC version:')
print('  {0}'.format(str(version))) 
masterindex = struct.unpack('<q',req.read(8))[0]
print('Master index:')
print('  {0}'.format(str(masterindex)))
genome = ""
c=req.read(1).decode("utf-8") 
while (c != '\0'):
    genome += c
    c=req.read(1).decode("utf-8") 
print('Genome ID:')
print('  {0}'.format(str(genome))) 
# read and throw away attribute dictionary (stats+graphs)
print('Attribute dictionary:')
nattributes = struct.unpack('<i',req.read(4))[0]
for x in range(0, nattributes):
  key = readcstr(req)
  value = readcstr(req)
#  print('   Key:{0}'.format(key))
#  print('   Value:{0}'.format(value))
nChrs = struct.unpack('<i',req.read(4))[0]
print("Chromosomes: ")
for x in range(0, nChrs):
  name = readcstr(req)
  length = struct.unpack('<i',req.read(4))[0]
  print('  {0}  {1}'.format(name, length))
nBpRes = struct.unpack('<i',req.read(4))[0]
print("Base pair-delimited resolutions: ")
for x in range(0, nBpRes):
  res = struct.unpack('<i',req.read(4))[0]
  print('   {0}'.format(res))
nFrag = struct.unpack('<i',req.read(4))[0]
print("Fragment-delimited resolutions: ")
for x in range(0, nFrag):
  res = struct.unpack('<i',req.read(4))[0]
  print('   {0}'.format(res))
