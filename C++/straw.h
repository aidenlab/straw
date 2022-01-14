/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#ifndef STRAW_H
#define STRAW_H

#include <fstream>
#include <set>
#include <vector>
#include <map>

// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
    int64_t size;
    int64_t position;
};

// sparse matrixType entry
struct contactRecord {
  int32_t binX;
  int32_t binY;
  float counts;
};

// chromosome
struct chromosome {
    std::string name;
    int32_t index;
    int64_t length;
};

std::map<int32_t, indexEntry>
readMatrixZoomData(std::istream &fin, const std::string &myunit, int32_t mybinsize, float &mySumCounts,
                   int32_t &myBlockBinCount,
                   int32_t &myBlockColumnCount, bool &found);

std::map<int32_t, indexEntry>
readMatrix(std::istream &fin, int32_t myFilePosition, std::string unit, int32_t resolution, float &mySumCounts,
           int32_t &myBlockBinCount, int32_t &myBlockColumnCount);

std::vector<double> readNormalizationVector(std::istream &fin, indexEntry entry);

std::vector<contactRecord>
straw(const std::string& matrixType, const std::string& norm, const std::string& fname, const std::string& chr1loc, const std::string& chr2loc,
      const std::string &unit, int32_t binsize);

#endif
