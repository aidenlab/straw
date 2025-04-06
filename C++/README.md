# Straw C++

This is a C++ implementation of the Straw tool for reading Hi-C data and converting it to a simpler slice format.

## Description:
The tool provides functionality to read .hic files and extract contact matrices at specific resolutions. It can output the data in a simplified binary slice format that is more efficient for downstream processing at single high-resolutions.

## Installation:
1. Requires CMake 3.13 or higher
2. Requires libcurl and zlib development libraries
3. Clone the repository
4. Create a build directory: `mkdir build`
5. Enter build directory: `cd build`
6. Run cmake: `cmake ..`
7. Build: `make`

## Usage:
The main executable 'straw' supports two modes:
1. Standard mode:
`straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG/MATRIX> <binsize>`
2. Dump mode (creates slice file):
`straw dump <observed/oe/expected> <NONE/VC/VC_SQRT/KR> <hicFile> <BP/FRAG> <binsize> <outputFile>`

## Examples:
1. Extract specific region:
`straw observed NONE input.hic chr1:0:1000000 chr2:0:1000000 BP 10000`
2. Create slice file at 10kb resolution:
`straw dump observed NONE input.hic BP 10000 output.slc`

## Slice Format:
The slice format (.slc) is a binary format that contains:
1. Magic string "HICSLICE"
2. Resolution (int32)
3. Number of chromosomes (int32)
4. Chromosome mapping (name lengths, names, and keys)
5. Contact records (chr1Key, binX, chr2Key, binY, value)

## Reading Slice Files:
A C++ reader is provided in the slice_reader directory. It provides methods to:
1. Read basic file information (resolution, chromosomes)
2. Read all contact records
3. Read records for specific chromosome pairs
The reader automatically handles the chromosome ordering convention.

## Notes:
The simplified slice format and reader is only intended for repeated analysis on a high resolution slice of the matrix. Otherwise, the original hic file format is more efficient.

## Bug Reports or Feature Requests:
For bug reports or feature requests, please open an issue on the repository.