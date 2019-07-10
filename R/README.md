# strawr
Straw is a fast implementation of reading/dump for .hic files

## Installation
```R
remotes::install_github("aidenlab/straw/R")
```

## Usage
```R
# <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
hic.data.frame <- strawr::straw_R("KR /path/to/file.hic 11 11 BP 10000")
```
