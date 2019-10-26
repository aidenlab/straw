from __future__ import print_function
#Reads the genome name from the hic header
#Can take in a .hic file or URL that points to .hic file
from typing import List
import sys
import struct
import requests
import io
import argparse

def readcstr(f):
  """
  Function to parse the requested .hic file.
  ----------
  Parameters:
    f: io.BytesIO object representation of the requested .hic file.
  ---------- 
  Returns:
    string version of the decoded .hic file.
  """
  buf = ""
  while True:
      b = f.read(1)
      b = b.decode('utf-8', 'backslashreplace')
      if b is None or b == '\0':
        return str(buf)
      else:
        buf = buf + b
#            buf.append(b.decode('utf'))

def main(
  infile:str, 
  verbose:bool=False,
  noprint_attributes:List[str]=None):
  """
  Main function used to read and print the header of a .hic file.
  ----------
  Parameters:
    infile: string path to .hic file (can also be an URL).
    verbose: boolean indicator of verbosity, default is False.
    noprint_attributes: list of strings representing keys in the .hic archive's attribute dictionary
                        which should not be printed.
                        None is the "sentinel default value", gets overwritten to ['graphs']
  ----------
  Returns:
    True if the function ran to completion.
  ----------
  Raises:
    requests.HTTPError: if infile is a URL, but requests.Session.get() returns an invalid (>=400) status code.
    ValueError: if the header of the .hic file does not contain the 'HIC' keyword.
  """
  # sentinel value replacement
  if noprint_attributes is None:
    noprint_attributes = ['graphs']

  magic_string = ""
  if (infile.startswith("http")):
      # try URL first. 100K should be sufficient for header
      headers={'range' : 'bytes=0-100000', 'x-amz-meta-requester' : 'straw'}
      s = requests.Session()
      r=s.get(infile, headers=headers)
      if verbose: 
        print("requested access to:\n" + infile)
      if (r.status_code >= 400):
        error_message = "Error accessing\n" + infile + "\nHTTP status code " + str(r.status_code)
        raise requests.HTTPError(error_message)
      req=io.BytesIO(r.content)        
      myrange=r.headers['content-range'].split('/')
      totalbytes=myrange[1]
  else:
      req=open(infile, 'rb')

  magic_string = struct.unpack('<3s', req.read(3))[0]
  req.read(1)
  if (magic_string != b"HIC"):
    error_message = 'This does not appear to be a .hic file; magic string is incorrect'
    raise ValueError(error_message)

  version = struct.unpack('<i',req.read(4))[0]
  if verbose:
    print('HiC version:')
    print('  {0}'.format(str(version))) 
  masterindex = struct.unpack('<q',req.read(8))[0]
  if verbose:
    print('Master index:')
    print('  {0}'.format(str(masterindex)))
  genome = ""
  c=req.read(1).decode("utf-8") 
  while (c != '\0'):
    genome += c
    c=req.read(1).decode("utf-8") 
  if verbose:
    print('Genome ID:')
    print('  {0}'.format(str(genome))) 
  
  # read and throw away attribute dictionary (stats+graphs)
  if verbose:
    print('Attribute dictionary:')
  nattributes = struct.unpack('<i',req.read(4))[0]
  for _ in range(0, nattributes):
    key = readcstr(req)
    value = readcstr(req)
    if verbose:
      print('   Key:{0}'.format(key))

      if key in noprint_attributes:
        print('   (not displaying this attribute field as specified)')
        continue 
      
      print('   Value:{0}'.format(value))
  nChrs = struct.unpack('<i',req.read(4))[0]
  if verbose:
    print("Chromosomes: ")
  for _ in range(0, nChrs):
    name = readcstr(req)
    length = struct.unpack('<i',req.read(4))[0]
    if verbose:
      print('  {0}  {1}'.format(name, length))
  nBpRes = struct.unpack('<i',req.read(4))[0]
  if verbose:
    print("Base pair-delimited resolutions: ")
  for _ in range(0, nBpRes):
    res = struct.unpack('<i',req.read(4))[0]
    if verbose:
      print('   {0}'.format(res))
  nFrag = struct.unpack('<i',req.read(4))[0]
  if verbose:
    print("Fragment-delimited resolutions: ")
  for _ in range(0, nFrag):
    res = struct.unpack('<i',req.read(4))[0]
    if verbose:
      print('   {0}'.format(res))
  req.close()
  return True 

if __name__ == '__main__':

  parser = argparse.ArgumentParser(
    prog="read_hic_header", 
    #usage="python read_hic_header.py <file path/URL> [--silent]", 
    description="Python script used to read and display the header of .hic files.\nSilent mode (--silent/-s) can be used to parse the file header as a check for .hic file integrity.",
    formatter_class=argparse.RawTextHelpFormatter
  )
  filepath_help_string = '''string path to .hic file (can also be an URL).
e.g.  read_hic_header.py /home/usr/xyz/hic_files/mm10_wt.hic
      read_hic_header.py https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/HIC001.hic'''

  silent_help_string = '''boolean switch for silent mode (minimal stdout output, checks file header for file integrity.).
e.g.  read_hic_header.py /home/usr/xyz/hic_files/mm10_wt.hic -s
      read_hic_header.py https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/HIC001.hic --silent'''

  parser.add_argument( 'filepath', type=str, metavar='filepath',
                      help=filepath_help_string)
  parser.add_argument( '--silent', '-s', action='store_false', required=False, default=True,
                      help=silent_help_string)
  args = parser.parse_args()
  main(args.filepath, args.silent)  