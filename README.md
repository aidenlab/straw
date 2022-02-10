# Straw: rapidly stream data from .hic files
Straw is library which allows rapid streaming of contact data from .hic files. 
This repository contains source code for the C++, R, Python, and MATLAB versions of Straw.

- For the Java version, see: https://github.com/aidenlab/java-straw/
- For the Javascript version, see: https://github.com/igvteam/hic-straw/

There are two Python versions - a pure Python flavor and one which wraps the C++ code with pybind11. 
The former version has been deprecated in favor of using the pybind11 version, which is much faster.

- For archival purposes, the pure python version has been moved to: https://github.com/aidenlab/pystraw

A Jupyter notebook example of using straw can be found here: https://aidenlab.gitbook.io/juicebox/accessing-raw-data

## Install straw for python

Use `pip install hic-straw`. 
If you want to build from the source code, you must have pybind11 installed. 
Clone the library and `cd` into the `straw/` directory. Then `pip install ./pybind11_python`.

## Compile straw for C++

```bash
g++ -std=c++0x -o straw main.cpp straw.cpp -lcurl -lz
```

You must have cURL installed.
Please see [the wiki](https://github.com/aidenlab/straw/wiki) for more documentation.

For questions, please use
[the Google Group](https://groups.google.com/forum/#!forum/3d-genomics).

Ongoing development work is carried out by <a href="http://mshamim.com">Muhammad S. Shamim</a>.
Past contributors include <a href="http://www.cherniavsky.net/neva/">Neva C. Durand</a> and many others.

If you use this tool in your work, please cite 

**Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov, Eric S. Lander, and Erez Lieberman Aiden. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell Systems 3(1), 2016.**

