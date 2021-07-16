# Straw: rapidly stream data from .hic files
Straw is library which allows rapid streaming of contact data from .hic files. 
This repository contains source code for the C++, R, Python, and MATLAB versions of Straw.

For the Java version of Straw, see: https://github.com/aidenlab/java-straw/

For the Javascript version of Straw, see: https://github.com/igvteam/hic-straw/

There are two Python versions - a pure Python one and a version where the C++ code is bound with pybind11. The former version is deprecated in favor of using the pybind11 version, which is much faster.

A Jupyter notebook example can be found here: https://aidenlab.gitbook.io/juicebox/accessing-raw-data

## Quick Start Python

For the fastest version, you must have pybind11 installed.

Clone the library and cd into the `straw/` directory.
```
pip install ./pybind11_python
```
Then run via `import strawC` and `strawC.strawC` 

```
    Example:
    >>>import strawC
    >>>result = strawC.strawC('NONE', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
    >>>for i in range(len(result)):
    ...   print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
```

## Compile on Linux

         g++ -std=c++0x -o straw main.cpp straw.cpp -lcurl -lz
 
Please see [the wiki](https://github.com/theaidenlab/straw/wiki) for more documentation.

For questions, please use
[the Google Group](https://groups.google.com/forum/#!forum/3d-genomics).

Ongoing development work is carried out by <a href="http://mshamim.com">Muhammad S. Shamim</a>, <a href="https://github.com/cwenger">Craig Wenger</a>, and <a href="http://www.cherniavsky.net/neva/">Neva C. Durand</a>.

If you use this tool in your work, please cite 

**Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov, Eric S. Lander, and Erez Lieberman Aiden. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell Systems 3(1), 2016.**

