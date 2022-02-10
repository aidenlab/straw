## Quick Start Python

Straw is library which allows rapid streaming of contact data from .hic files. 
To learn more about Hi-C data and 3D genomics, visit https://aidenlab.gitbook.io/juicebox/

Once you've installed the library with `pip install hic-straw`, you can import your code with `import hicstraw`. 

## Legacy usage to fetch list of contacts

The new usage for straw allows you to create objects and retain intermediate variables.
This can speed up your code significantly when querying hundreds or thousands of regions
for a given chromosome/resolution/normalization.

First we import `numpy` and `hicstraw`.
```python
import numpy as np
import hicstraw
```

We then create a Hi-C file object. 
From this object, we can query genomeID, chromosomes, and resolutions.
```python
hic = hicstraw.HiCFile("HIC001.hic")
print(hic.getChromosomes())
print(hic.getGenomeID())
print(hic.getResolutions())
```

We can also collect a matrix zoom data object, which is specific to 
- specific matrix-type: `observed` (count) or `oe` (observed/expected ratio)
- chromosome-chromosome pair
- resolution
- normalization

This object retains information for fast future queries. 
Here's an example that pick the counts from the intrachromosomal region for chr4 
with KR normalization at 5kB resolution.
```python
mzd = hic.getMatrixZoomData('4', '4', "observed", "KR", "BP", 5000)
```

We can get numpy matrices for specific genomic windows by calling:
```python
numpy_matrix = mzd.getRecordsAsMatrix(10000000, 12000000, 10000000, 12000000)
```

### Usage
```
hic = hicstraw.HiCFile(filepath)
hic.getChromosomes()
hic.getGenomeID()
hic.getResolutions()

mzd = hic.getMatrixZoomData(chrom1, chrom2, data_type, normalization, "BP", resolution)

numpy_matrix = mzd.getRecordsAsMatrix(gr1, gr2, gc1, gc2)
records_list = mzd.getRecords(gr1, gr2, gc1, gc2)
```

`filepath`: path to file (local or URL)<br>
`data_type`: `'observed'` (previous default / "main" data) or `'oe'` (observed/expected)<br>
`normalization`: `NONE`, `VC`, `VC_SQRT`, `KR`, `SCALE`, etc.<br>
`resolution`: typically `2500000`, `1000000`, `500000`, `100000`, `50000`, `25000`, `10000`, `5000`, etc.<br><br>
Note: the normalization, resolution, and chromosome/regions must already exist in the .hic to be read 
(i.e. they are not calculated by straw, only read from the file if available)<br>
`gr1`: start genomic position along rows<br>
`gr2`: end genomic position along rows<br>
`gc1`: start genomic position along columns<br>
`gc2`: end genomic position along columns<br>


## Legacy usage to fetch list of contacts

For example, to fetch a list of all the raw contacts on chrX at 100Kb resolution:

```python
import hicstraw
result = hicstraw.straw('observed', 'NONE', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
for i in range(len(result)):
    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
```

To fetch a list of KR normalized contacts for the same region:
```python
import hicstraw
result = hicstraw.straw('observed', 'KR', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
for i in range(len(result)):
    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
```

To query observed/expected KR normalized data:
```python
import hicstraw
result = hicstraw.straw('oe', 'KR', 'HIC001.hic', 'X', 'X', 'BP', 1000000)
for i in range(len(result)):
    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))
```

### Usage
```
hicstraw.straw(data_type, normalization, file, region_x, region_y, 'BP', resolution)
```

`data_type`: `'observed'` (previous default / "main" data) or `'oe'` (observed/expected)<br>
`normalization`: `NONE`, `VC`, `VC_SQRT`, `KR`, `SCALE`, etc.<br>
`file`: filepath (local or URL)<br>
`region_x/y`: provide the `chromosome` or utilize the syntax `chromosome:start_position:end_position` if using a smaller window within the chromosome<br>
`resolution`: typically `2500000`, `1000000`, `500000`, `100000`, `50000`, `25000`, `10000`, `5000`, etc.<br><br>
Note: the normalization, resolution, and chromosome/regions must already exist in the .hic to be read 
(i.e. they are not calculated by straw, only read from the file if available)<br>

