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

// this is for creating a stream from a byte array for ease of use
// see https://stackoverflow.com/questions/41141175/how-to-implement-seekg-seekpos-on-an-in-memory-buffer
struct membuf : std::streambuf {
    membuf(char *begin, int32_t l) {
        setg(begin, begin, begin + l);
    }
};

struct memstream : virtual membuf, std::istream {
    memstream(char *begin, int32_t l) :
            membuf(begin, l),
            std::istream(static_cast<std::streambuf*>(this)) {
    }

    std::istream::pos_type seekpos(std::istream::pos_type sp, std::ios_base::openmode which) override {
        return seekoff(sp - std::istream::pos_type(std::istream::off_type(0)), std::ios_base::beg, which);
    }

    std::istream::pos_type seekoff(std::istream::off_type off,
                                    std::ios_base::seekdir dir,
                                    std::ios_base::openmode which = std::ios_base::in) override {
        if (dir == std::ios_base::cur)
            gbump(off);
        else if (dir == std::ios_base::end)
            setg(eback(), egptr() + off, egptr());
        else if (dir == std::ios_base::beg)
            setg(eback(), eback() + off, egptr());
        return gptr() - eback();
    }
};

// for holding data from URL call
struct MemoryStruct {
    char *memory;
    size_t size;
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
straw(const std::string& matrixType, const std::string& norm, const std::string& fname, const std::string& chr1loc,
      const std::string& chr2loc, const std::string &unit, int32_t binsize);

std::vector<std::vector<float>>
strawAsMatrix(const std::string &matrixType, const std::string &norm, const std::string &fileName,
              const std::string &chr1loc, const std::string &chr2loc, const std::string &unit, int32_t binsize);

int64_t getNumRecordsForFile(const std::string& filename, int32_t binsize, bool interOnly);

int64_t getNumRecordsForChromosomes(const std::string& filename, int32_t binsize, bool interOnly);

class HiCFile {
public:  // Make sure this is public
    string prefix = "http"; // HTTP code
    int64_t master = 0LL;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int32_t numChromosomes = 0;
    int32_t version = 0;
    int64_t nviPosition = 0LL;
    int64_t nviLength = 0LL;
    vector<int32_t> resolutions;
    static int64_t totalFileSize;
    string fileName;

    explicit HiCFile(const string &fileName);
    string getGenomeID() const;
    vector<int32_t> getResolutions() const;
    vector<chromosome> getChromosomes();
    MatrixZoomData* getMatrixZoomData(const string &chr1, const string &chr2, const string &matrixType,
                                    const string &norm, const string &unit, int32_t resolution);

    static size_t hdf(char *b, size_t size, size_t nitems, void *userdata);
    static CURL *oneTimeInitCURL(const char *url);
};

class MatrixZoomData {
public:  // Make sure this is public
    bool isIntra;
    string fileName;
    int64_t myFilePos = 0LL;
    vector<double> expectedValues;
    bool foundFooter = false;
    vector<double> c1Norm;
    vector<double> c2Norm;
    int32_t c1 = 0;
    int32_t c2 = 0;
    string matrixType;
    string norm;
    int32_t version = 0;
    int32_t resolution = 0;
    int32_t numBins1 = 0;
    int32_t numBins2 = 0;
    float sumCounts;
    int32_t blockBinCount, blockColumnCount;
    map<int32_t, indexEntry> blockMap;
    double avgCount;

    MatrixZoomData(const chromosome &chrom1, const chromosome &chrom2, const string &matrixType,
                   const string &norm, const string &unit, int32_t resolution,
                   int32_t &version, int64_t &master, int64_t &totalFileSize,
                   const string &fileName);

    vector<contactRecord> getRecords(int64_t gx0, int64_t gx1, int64_t gy0, int64_t gy1);
    // ... other public methods ...
};

#endif
