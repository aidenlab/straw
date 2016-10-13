# straw: Extract data quickly from Juicebox via straw
Straw is a fast C++ implementation of reading/dump for .hic files.
 
Please see [the wiki](wiki) for more documentation.

For questions, please use
[the Google Group](https://groups.google.com/forum/#!forum/3d-genomics).

If you use this tool in your work, please cite 

**Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov, Eric S. Lander, and Erez Lieberman Aiden. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell Systems 3(1), 2016.**

## C++
Straw must be compiled with the -lz flag to include the zlib.h library:
`g++ -lz -o straw main.cpp straw.cpp`

Usage: `straw <NONE/VC/VC_SQRT/KR> <hicFile> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize> `

Users have also reported needing this flag: -std=c++11
## R
Straw is compatible with R via the Rcpp library.  Usage is very similar to the above.  Future versions might use a different signature to call the function.

```
library(Rcpp)
sourceCpp("straw-R.cpp")
A<-straw_R("NONE /Users/nchernia/Downloads/drosophila.hic arm_2L arm_2L BP 100000")
```
In the above example, A is a data frame containing the counts information:

```
> head(A)
        x        y counts
1  800000  4700000     77
2       0        0  41105
3 5500000  8100000     82
4 2700000  8300000     85
5 5700000 17700000     11
6 3000000 16500000     27
```

## Python 
Straw can be called from Python by using Boost <http://www.boost.org/>

The files [Jamroot](Jamroot), [straw.py](straw.py), and 
[main-python.cpp](main-python.cpp) are used to hook the C++ to Python. 
[main-python.cpp](main-python.cpp) is essentially the same as 
[main.cpp](main.cpp)

1. Get Boost:  
  http://www.boost.org/doc/libs/1_61_0/more/getting_started/unix-variants.html#get-boost
2. Compile bjam:  
  http://www.boost.org/doc/libs/1_49_0/more/getting_started/unix-variants.html#prepare-to-use-a-boost-library-binary
(just do --with-libraries python)
3. Type "bjam" in the directory with the files (assuming you've put boost and bjam in your path).
