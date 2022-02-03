# Straw: rapidly stream data from .hic files
Straw is library which allows rapid streaming of contact data from .hic files. 
This repository contains source code for the C++, R, Python, and MATLAB versions of Straw.

For the Java version of Straw, see: https://github.com/aidenlab/java-straw/

For the Javascript version of Straw, see: https://github.com/igvteam/hic-straw/

There are two Python versions - a pure Python one and a version where the C++ code is bound with pybind11. The former version is deprecated in favor of using the pybind11 version, which is much faster.

A Jupyter notebook example can be found here: https://aidenlab.gitbook.io/juicebox/accessing-raw-data

## Quick Start Python

Use `pip install hic-straw`. Or if you want to build from the source code, you must have pybind11 installed. Clone the library and `cd` into the `straw/` directory. Then `pip install ./pybind11_python`.

You can run your code via `import strawC` (or `hic-straw`) and `strawC.strawC`, for example:

```python
import strawC
result = strawC.strawC('observed', 'NONE', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
for i in range(len(result)):
    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
```
To query observed/expected data:
```python
import strawC
result = strawC.strawC('oe', 'NONE', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
for i in range(len(result)):
    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
```

### Usage
```
strawC.strawC(data_type, normalization, file, region_x, region_y, 'BP', resolution)
```

`data_type`: `'observed'` (previous default / "main" data) or `'oe'` (observed/expected)<br>
`normalization`: `NONE`, `VC`, `VC_SQRT`, `KR`, `SCALE`, etc.<br>
`file`: filepath (local or URL)<br>
`region_x/y`: provide the `chromosome` or utilize the syntax `chromosome:start_position:end_position` if using a smaller window within the chromosome<br>
`resolution`: typically `2500000`, `1000000`, `500000`, `100000`, `50000`, `25000`, `10000`, `5000`, etc.<br><br>
Note: the normalization, resolution, and chromosome/regions must already exist in the .hic to be read 
(i.e. they are not calculated by straw, only read from the file if available)<br>


## Compile on Linux

```bash
g++ -std=c++0x -o straw main.cpp straw.cpp -lcurl -lz
```

Please see [the wiki](https://github.com/theaidenlab/straw/wiki) for more documentation.

For questions, please use
[the Google Group](https://groups.google.com/forum/#!forum/3d-genomics).

Ongoing development work is carried out by <a href="http://mshamim.com">Muhammad S. Shamim</a>.

If you use this tool in your work, please cite 

**Neva C. Durand, James T. Robinson, Muhammad S. Shamim, Ido Machol, Jill P. Mesirov, Eric S. Lander, and Erez Lieberman Aiden. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell Systems 3(1), 2016.**

