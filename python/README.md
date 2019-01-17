Straw enables programmatic access to .hic files.
.hic files store the contact matrices from Hi-C experiments and the
normalization and expected vectors, along with meta-data in the header.

The main function, straw, takes in the normalization, the filename or URL,
chromosome1 (and optional range), chromosome2 (and optional range),
whether the bins desired are fragment or base pair delimited, and bin size.

It then reads the header, follows the various pointers to the desired matrix
and normalization vector, and stores as [x, y, count]

Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <\
BP/FRAG> <binsize>

Example:
>>>import straw
>>>result = straw.straw('NONE', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
>>>for i in range(len(result[0])):
...   print("{0}\t{1}\t{2}".format(result[0][i], result[1][i], result[2][i]))

See https://github.com/theaidenlab/straw/wiki/Python for more documentation

# Install
Straw uses the [requests library](http://docs.python-requests.org/en/master/user/install/#install) for support of URLs.  Be sure it is installed.

# Examples

* Extract all reads on chromosome X at 1MB resolution with no normalization in local file "HIC001.hic" 
   ```python
   import straw
   result = straw.straw('NONE', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
   # the values returned are in x / y / counts
   for i in range(len(result[0])):
      print("{0}\t{1}\t{2}".format(result[0][i], result[1][i], result[2][i]))
   ```

* Extract reads between 1MB and 7.5MB on chromosome 1 at 25KB resolution with KR (balanced) normalization and write to a file:
   ```python
   import straw
   straw.printme("KR", "HIC001.hic", "1:1000000:7500000", "1:1000000:7500000", "BP", 25000, 'out.txt')
   ```
To get in an array:
   ```python
   import straw
   result = straw.straw("KR", "HIC001.hic", "1:1000000:7500000", "1:1000000:7500000", "BP", 25000)
   ```

* Extract all interchromosomal reads between chromosome 5 and chromosome 12 at 500 fragment resolution with VC (vanilla coverage) normalization:
   ```python
   import straw

   result = straw.straw("VC", "HIC001.hic", "5", "12", "FRAG", 500)
   # the values returned are in results
   for i in range(len(result[0])):
      print("{0}\t{1}\t{2}".format(result[0][i], result[1][i], result[2][i]))
   ```

See the script [straw.py](https://github.com/theaidenlab/straw/blob/master/python/straw.py) for an example of how to print the results to a file.  
# Read header
See the file [read_hic_header.py](https://github.com/theaidenlab/straw/blob/master/python/read_hic_header.py) for a Python script that reads the header of a hic file and outputs the information (including resolutions).
