#Reads the genome name from the hic header
#Can take in a .hic file or URL that points to .hic file
import sys
import struct
import urllib

def readcstr(f):
    buf = bytearray()
    while True:
        b = f.read(1)
        if b is None or b == '\0':
            return str(buf)
        else:
            buf.append(b)

if (len(sys.argv) != 2):
  sys.stderr.write('Usage: '+ sys.argv[0]+' <hic file or URL>\n')
  sys.exit(1)

infile = sys.argv[1]
magic_string = ""
try: req = urllib.urlopen(infile)
except:
  req=open(infile, 'r')
c=req.read(1)
while (c != '\0'):
    magic_string += c
    c=req.read(1)
if (magic_string != "HIC"):
    print 'This does not appear to be a HiC file; magic string is incorrect'
    sys.exit(1)
version = struct.unpack('<i',req.read(4))[0]
print 'HiC version:'
print '  {0}'.format(str(version)) 
masterindex = struct.unpack('<q',req.read(8))[0]
genome = ""
c=req.read(1)
while (c != '\0'):
    genome += c
    c=req.read(1)
print 'Genome ID:'
print '  {0}'.format(str(genome)) 
# read and throw away attribute dictionary (stats+graphs)
print 'Attribute dictionary:'
nattributes = struct.unpack('<i',req.read(4))[0]
for x in range(0, nattributes):
  key = readcstr(req)
  value = readcstr(req)
#  print '   Key:{0}'.format(key)
#  print '   Value:{0}'.format(value)
nChrs = struct.unpack('<i',req.read(4))[0]
print "Chromosomes: "
for x in range(0, nChrs):
  name = readcstr(req)
  length = struct.unpack('<i',req.read(4))[0]
  print '  {0}  {1}'.format(name, length)
nBpRes = struct.unpack('<i',req.read(4))[0]
print "Base pair-delimited resolutions: "
for x in range(0, nBpRes):
  res = struct.unpack('<i',req.read(4))[0]
  print '   {0}'.format(res)
nFrag = struct.unpack('<i',req.read(4))[0]
print "Fragment-delimited resolutions: "
for x in range(0, nFrag):
  res = struct.unpack('<i',req.read(4))[0]
  print '   {0}'.format(res)
