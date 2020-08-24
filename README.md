# straw: Extract data quickly from Juicebox via straw
Straw is a fast implementation of reading/dump for .hic files. Available in C++, R, and Python.

There are two Python versions - a pure Python one and a version where the C++ code is bound with pybind11. The latter is much faster.

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

For the pure Python version:
```
pip install ./python
```
Then run via `import straw` and `straw.straw`

Some more information on running: https://github.com/aidenlab/straw/tree/master/python

## Compile on Linux

         g++ -std=c++0x -o straw main.cpp straw.cpp -lcurl -lz
 
Please see [the wiki](https://github.com/theaidenlab/straw/wiki) for more documentation.

For questions, please use
[the Google Group](https://groups.google.com/forum/#!forum/3d-genomics).

Ongoing development work is carried out by <a href="http://www.cherniavsky.net/neva/">Neva C. Durand</a> and <a href="https://mikeaalv.github.io/">Yue Wu</a>.

If you use this tool in your work, please cite 

**Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov, Eric S. Lander, and Erez Lieberman Aiden. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell Systems 3(1), 2016.**

