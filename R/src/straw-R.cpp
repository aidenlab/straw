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
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <streambuf>
#include "zlib.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*
  Quick dump: fast C++ implementation of dump. Not as fully featured as the
  Java version. Reads the .hic file, finds the appropriate matrix and slice
  of data, and outputs as text in sparse upper triangular format.

  Currently only supporting matrices.

  Usage: juicebox-quick-dump <observed/oe> <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
 */
// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
  int32_t size;
  int64_t position;
};

// sparse matrix entry
struct contactRecord {
  int32_t binX;
  int32_t binY;
  float counts;
};

// this is for creating a stream from a byte array for ease of use
struct membuf : std::streambuf
{
  membuf(char* begin, char* end) {
    this->setg(begin, begin, end);
  }
};

// version number
int32_t version;

// map of block numbers to pointers
map <int32_t, indexEntry> blockMap;

// returns whether or not this is valid HiC file
bool readMagicString(ifstream& fin) {
  string str;
  getline(fin, str, '\0' );
  return str[0]=='H' && str[1]=='I' && str[2]=='C';
}

// reads the header, storing the positions of the normalization vectors and returning the master pointer
int64_t readHeader(ifstream& fin, string chr1, string chr2, int32_t &c1pos1, int32_t &c1pos2, int32_t &c2pos1, int32_t &c2pos2, int32_t &chr1ind, int32_t &chr2ind, int32_t &chr1len, int32_t &chr2len) {
  if (!readMagicString(fin)) {
    stop("Hi-C magic string is missing, does not appear to be a hic file.");
    // exit(1);
  }

  fin.read((char*)&version, sizeof(int32_t));
  if (version < 6) {
    stop("Version %d no longer supported.", version);
    // exit(1);
  }
  int64_t master;
  fin.read((char*)&master, sizeof(int64_t));
  string genome;
  getline(fin, genome, '\0' );
  int32_t nattributes;
  fin.read((char*)&nattributes, sizeof(int32_t));
  // reading and ignoring attribute-value dictionary
  for (int i=0; i<nattributes; i++) {
    string key, value;
    getline(fin, key, '\0');
    getline(fin, value, '\0');
  }
  int32_t nChrs;
  fin.read((char*)&nChrs, sizeof(int32_t));
  // chromosome map for finding matrix
  bool found1 = false;
  bool found2 = false;
  for (int i=0; i<nChrs; i++) {
    string name;
    int32_t length;
    getline(fin, name, '\0');
    fin.read((char*)&length, sizeof(int32_t));
    if (name==chr1) {
      found1=true;
      chr1ind = i;
      chr1len = length;
      if (c1pos1 == -100) {
	c1pos1 = 0;
	c1pos2 = length;
      }
    }
    if (name==chr2) {
      found2=true;
      chr2ind = i;
      chr2len = length;
      if (c2pos1 == -100) {
	c2pos1 = 0;
	c2pos2 = length;
      }
    }
  }
  if (!found1 || !found2) {
    stop("One of the chromosomes wasn't found in the file. Check that the chromosome name matches the genome.");
    // exit(1);
  }
  return master;
}

// reads the footer from the master pointer location. takes in the chromosomes, norm, unit (BP or FRAG) and resolution or
// binsize, and sets the file position of the matrix and the normalization vectors for those chromosomes at the given
// normalization and resolution
void readFooter(ifstream& fin, int64_t master, int32_t c1, int32_t c2, string matrix, string norm, string unit, int32_t resolution, int64_t &myFilePos, indexEntry &c1NormEntry, indexEntry &c2NormEntry, vector<double> &expectedValues) {
  fin.seekg(master, ios::beg);
  int32_t nBytes;
  fin.read((char*)&nBytes, sizeof(int32_t));

  stringstream ss;
  ss << c1 << "_" << c2;
  string key = ss.str();

  int32_t nEntries;
  fin.read((char*)&nEntries, sizeof(int32_t));
  bool found = false;
  for (int i=0; i<nEntries; i++) {
    string str;
    getline(fin, str, '\0');
    int64_t fpos;
    fin.read((char*)&fpos, sizeof(int64_t));
    int32_t sizeinbytes;
    fin.read((char*)&sizeinbytes, sizeof(int32_t));
    if (str == key) {
      myFilePos = fpos;
      found=true;
    }
  }
  if (!found) {
    stop("File doesn't have the given chr_chr map.");
    // exit(1);
  }

  if ((matrix=="observed" && norm=="NONE") || ((matrix=="oe" || matrix=="expected") && norm=="NONE" && c1!=c2)) return; // no need to read norm vector index

  // read in and ignore expected value maps; don't store; reading these to
  // get to norm vector index
  int32_t nExpectedValues;
  fin.read((char*)&nExpectedValues, sizeof(int32_t));
  for (int i=0; i<nExpectedValues; i++) {
    string unit0;
    getline(fin, unit0, '\0'); //unit
    int32_t binSize;
    fin.read((char*)&binSize, sizeof(int32_t));

    int32_t nValues;
    fin.read((char*)&nValues, sizeof(int32_t));
    bool store = c1 == c2 && matrix == "oe" && norm == "NONE" && unit0 == unit && binSize == resolution;
    for (int j=0; j<nValues; j++) {
      double v;
      fin.read((char*)&v, sizeof(double));
      if (store) {
        expectedValues.push_back(v);
      }
    }

    int32_t nNormalizationFactors;
    fin.read((char*)&nNormalizationFactors, sizeof(int32_t));
    for (int j=0; j<nNormalizationFactors; j++) {
      int32_t chrIdx;
      fin.read((char*)&chrIdx, sizeof(int32_t));
      double v;
      fin.read((char*)&v, sizeof(double));
      if (store && chrIdx == c1) {
        for (vector<double>::iterator it=expectedValues.begin(); it!=expectedValues.end(); ++it) {
          *it = *it / v;
        }
      }
    }
  }
  if (c1 == c2 && matrix == "oe" && norm == "NONE") {
    if (expectedValues.size() == 0) {
      stop("File did not contain expected values vectors at %d %s.", resolution, unit);
      // exit(1);
    }
    return;
  }
  fin.read((char*)&nExpectedValues, sizeof(int32_t));
  for (int i=0; i<nExpectedValues; i++) {
    string type;
    getline(fin, type, '\0'); //typeString
    string unit0;
    getline(fin, unit0, '\0'); //unit
    int32_t binSize;
    fin.read((char*)&binSize, sizeof(int32_t));

    int32_t nValues;
    fin.read((char*)&nValues, sizeof(int32_t));
    bool store = c1 == c2 && matrix == "oe" && type == norm && unit0 == unit && binSize == resolution;
    for (int j=0; j<nValues; j++) {
      double v;
      fin.read((char*)&v, sizeof(double));
      if (store) {
        expectedValues.push_back(v);
      }
    }
    int32_t nNormalizationFactors;
    fin.read((char*)&nNormalizationFactors, sizeof(int32_t));
    for (int j=0; j<nNormalizationFactors; j++) {
      int32_t chrIdx;
      fin.read((char*)&chrIdx, sizeof(int32_t));
      double v;
      fin.read((char*)&v, sizeof(double));
      if (store && chrIdx == c1) {
        for (vector<double>::iterator it=expectedValues.begin(); it!=expectedValues.end(); ++it) {
          *it = *it / v;
        }
      }
    }
  }
  if (c1 == c2 && matrix == "oe" && norm != "NONE") {
    if (expectedValues.size() == 0) {
      stop("File did not contain normalized expected values vectors at %d %s.", resolution, unit);
      // exit(1);
    }
  }
  // Index of normalization vectors
  fin.read((char*)&nEntries, sizeof(int32_t));
  bool found1 = false;
  bool found2 = false;
  for (int i = 0; i < nEntries; i++) {
    string normtype;
    getline(fin, normtype, '\0'); //normalization type
    int32_t chrIdx;
    fin.read((char*)&chrIdx, sizeof(int32_t));
    string unit1;
    getline(fin, unit1, '\0'); //unit
    int32_t resolution1;
    fin.read((char*)&resolution1, sizeof(int32_t));
    int64_t filePosition;
    fin.read((char*)&filePosition, sizeof(int64_t));
    int32_t sizeInBytes;
    fin.read((char*)&sizeInBytes, sizeof(int32_t));
    if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
      c1NormEntry.position=filePosition;
      c1NormEntry.size=sizeInBytes;
      found1 = true;
    }
    if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
      c2NormEntry.position=filePosition;
      c2NormEntry.size=sizeInBytes;
      found2 = true;
    }
  }
  if (!found1 || !found2) {
    stop("File did not contain %s normalization vectors for one or both chromosomes at %d %s.", norm, resolution, unit);
    // exit(1);
  }
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
bool readMatrixZoomData(ifstream& fin, string myunit, int32_t mybinsize, float &mySumCounts, int32_t &myBlockBinCount, int32_t &myBlockColumnCount) {
  string unit;
  getline(fin, unit, '\0' ); // unit
  int32_t tmp;
  fin.read((char*)&tmp, sizeof(int32_t)); // Old "zoom" index -- not used
  float sumCounts;
  fin.read((char*)&sumCounts, sizeof(float)); // sumCounts
  float tmp2;
  fin.read((char*)&tmp2, sizeof(float)); // occupiedCellCount
  fin.read((char*)&tmp2, sizeof(float)); // stdDev
  fin.read((char*)&tmp2, sizeof(float)); // percent95
  int32_t binSize;
  fin.read((char*)&binSize, sizeof(int32_t));
  int32_t blockBinCount;
  fin.read((char*)&blockBinCount, sizeof(int32_t));
  int32_t blockColumnCount;
  fin.read((char*)&blockColumnCount, sizeof(int32_t));

  bool storeBlockData = false;
  if (myunit==unit && mybinsize==binSize) {
    mySumCounts = sumCounts;
    myBlockBinCount = blockBinCount;
    myBlockColumnCount = blockColumnCount;
    storeBlockData = true;
  }

  int32_t nBlocks;
  fin.read((char*)&nBlocks, sizeof(int32_t));

  for (int b = 0; b < nBlocks; b++) {
    int32_t blockNumber;
    fin.read((char*)&blockNumber, sizeof(int32_t));
    int64_t filePosition;
    fin.read((char*)&filePosition, sizeof(int64_t));
    int32_t blockSizeInBytes;
    fin.read((char*)&blockSizeInBytes, sizeof(int32_t));
    indexEntry entry;
    entry.size = blockSizeInBytes;
    entry.position = filePosition;
    if (storeBlockData) blockMap[blockNumber] = entry;
  }
  return storeBlockData;
}

// goes to the specified file pointer and finds the raw contact matrix at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
void readMatrix(ifstream& fin, int64_t myFilePosition, string unit, int32_t resolution, float &mySumCounts, int32_t &myBlockBinCount, int32_t &myBlockColumnCount) {
  fin.seekg(myFilePosition, ios::beg);
  int32_t c1,c2;
  fin.read((char*)&c1, sizeof(int32_t)); //chr1
  fin.read((char*)&c2, sizeof(int32_t)); //chr2
  int32_t nRes;
  fin.read((char*)&nRes, sizeof(int32_t));
  int32_t i=0;
  bool found=false;
  while (i<nRes && !found) {
    found = readMatrixZoomData(fin, unit, resolution, mySumCounts, myBlockBinCount, myBlockColumnCount);
    i++;
  }
  if (!found) {
    stop("Error finding block data.");
    // exit(1);
  }
}
// gets the blocks that need to be read for this slice of the data.  needs blockbincount, blockcolumncount, and whether
// or not this is intrachromosomal.
set<int32_t> getBlockNumbersForRegionFromBinPosition(int32_t* regionIndices, int32_t blockBinCount, int32_t blockColumnCount, bool intra) {
   int32_t col1 = regionIndices[0] / blockBinCount;
   int32_t col2 = (regionIndices[1] + 1) / blockBinCount;
   int32_t row1 = regionIndices[2] / blockBinCount;
   int32_t row2 = (regionIndices[3] + 1) / blockBinCount;

   set<int32_t> blocksSet;
   // first check the upper triangular matrix
   for (int32_t r = row1; r <= row2; r++) {
     for (int32_t c = col1; c <= col2; c++) {
       int32_t blockNumber = r * blockColumnCount + c;
       blocksSet.insert(blockNumber);
     }
   }
   // check region part that overlaps with lower left triangle
   // but only if intrachromosomal
   if (intra) {
     for (int32_t r = col1; r <= col2; r++) {
       for (int32_t c = row1; c <= row2; c++) {
	 int32_t blockNumber = r * blockColumnCount + c;
	 blocksSet.insert(blockNumber);
       }
     }
   }

   return blocksSet;
}

// this is the meat of reading the data.  takes in the block number and returns the set of contact records corresponding to
// that block.  the block data is compressed and must be decompressed using the zlib library functions
vector<contactRecord> readBlock(ifstream& fin, int32_t blockNumber) {
  indexEntry idx = blockMap[blockNumber];
  if (idx.size == 0) {
    vector<contactRecord> v;
    return v;
  }
  char compressedBytes[idx.size];
  char* uncompressedBytes = new char[idx.size*10]; //biggest seen so far is 3
  fin.seekg(idx.position, ios::beg);
  fin.read(compressedBytes, idx.size);
  // Decompress the block
  // zlib struct
  z_stream infstream;
  infstream.zalloc = Z_NULL;
  infstream.zfree = Z_NULL;
  infstream.opaque = Z_NULL;
  infstream.avail_in = (uInt)(idx.size); // size of input
  infstream.next_in = (Bytef *)compressedBytes; // input char array
  infstream.avail_out = (uint32_t)idx.size*10; // size of output
  infstream.next_out = (Bytef *)uncompressedBytes; // output char array
  // the actual decompression work.
  inflateInit(&infstream);
  inflate(&infstream, Z_NO_FLUSH);
  inflateEnd(&infstream);
  int32_t uncompressedSize=infstream.total_out;

  // create stream from buffer for ease of use
  membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
  istream bufferin(&sbuf);
  int32_t nRecords;
  bufferin.read((char*)&nRecords, sizeof(int32_t));
  vector<contactRecord> v(nRecords);
  // different versions have different specific formats
  if (version < 7) {
    for (int32_t i = 0; i < nRecords; i++) {
      int32_t binX, binY;
      bufferin.read((char*)&binX, sizeof(int32_t));
      bufferin.read((char*)&binY, sizeof(int32_t));
      float counts;
      bufferin.read((char*)&counts, sizeof(float));
      contactRecord record;
      record.binX = binX;
      record.binY = binY;
      record.counts = counts;
      v[i] = record;
    }
  }
  else {
    int32_t binXOffset, binYOffset;
    bufferin.read((char*)&binXOffset, sizeof(int32_t));
    bufferin.read((char*)&binYOffset, sizeof(int32_t));
    char useShort;
    bufferin.read((char*)&useShort, sizeof(char));
    char type;
    bufferin.read((char*)&type, sizeof(char));
    int32_t index=0;
    if (type == 1) {
      // List-of-rows representation
      int16_t rowCount;
      bufferin.read((char*)&rowCount, sizeof(int16_t));
      for (int32_t i = 0; i < rowCount; i++) {
	int16_t y;
	bufferin.read((char*)&y, sizeof(int16_t));
	int32_t binY = y + binYOffset;
	int16_t colCount;
	bufferin.read((char*)&colCount, sizeof(int16_t));
	for (int32_t j = 0; j < colCount; j++) {
	  int16_t x;
	  bufferin.read((char*)&x, sizeof(int16_t));
	  int32_t binX = binXOffset + x;
	  float counts;
	  if (useShort == 0) { // yes this is opposite of usual
	    int16_t c;
	    bufferin.read((char*)&c, sizeof(int16_t));
	    counts = c;
	  }
	  else {
	    bufferin.read((char*)&counts, sizeof(float));
	  }
	  contactRecord record;
	  record.binX = binX;
	  record.binY = binY;
	  record.counts = counts;
	  v[index]=record;
	  index++;
	}
      }
    }
    else if (type == 2) { // have yet to find test file where this is true, possibly entirely deprecated
      int32_t nPts;
      bufferin.read((char*)&nPts, sizeof(int32_t));
      int16_t w;
      bufferin.read((char*)&w, sizeof(int16_t));

      for (int32_t i = 0; i < nPts; i++) {
	//int32_t idx = (p.y - binOffset2) * w + (p.x - binOffset1);
	int32_t row = i / w;
	int32_t col = i - row * w;
	int32_t bin1 = binXOffset + col;
	int32_t bin2 = binYOffset + row;

	float counts;
	if (useShort == 0) { // yes this is opposite of the usual
	  int16_t c;
	  bufferin.read((char*)&c, sizeof(int16_t));
	  if (c != -32768) {
	    contactRecord record;
	    record.binX = bin1;
	    record.binY = bin2;
	    record.counts = c;
	    v[index]=record;
	    index++;
	  }
	}
	else {
	  bufferin.read((char*)&counts, sizeof(float));
	  if (counts != 0x7fc00000) { // not sure this works
	    //	  if (!Float.isNaN(counts)) {
	    contactRecord record;
	    record.binX = bin1;
	    record.binY = bin2;
	    record.counts = counts;
	    v[index]=record;
	    index++;
	  }
	}
      }
    }
  }
  delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
  return v;
}

// reads the normalization vector from the file at the specified location
vector<double> readNormalizationVector(ifstream& fin, indexEntry entry) {
  char buffer[entry.size];
  fin.seekg(entry.position, ios::beg);
  fin.read(buffer, entry.size);
  membuf sbuf(buffer, buffer + entry.size);
  istream bufferin(&sbuf);
  int32_t nValues;
  bufferin.read((char*)&nValues, sizeof(int32_t));
  vector<double> values(nValues);
  //  bool allNaN = true;

  for (int i = 0; i < nValues; i++) {
    double d;
    bufferin.read((char*)&d, sizeof(double));
    values[i] = d;
    /* if (!Double.isNaN(values[i])) {
      allNaN = false;
      }*/
  }
  //  if (allNaN) return null;
  return values;
}

//' Straw Quick Dump
//'
//' fast C++ implementation of dump. Not as fully featured as the
//' Java version. Reads the .hic file, finds the appropriate matrix and slice
//' of data, and outputs as data.frame in sparse upper triangular format.
//' Currently only supporting matrices.
//'
//' Usage: straw <observed/oe> <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
//'
//' @param norm Normalization to apply. Must be one of NONE/VC/VC_SQRT/KR.
//'     VC is vanilla coverage, VC_SQRT is square root of vanilla coverage, and KR is Knight-Ruiz or
//'     Balanced normalization.
//' @param fname path to .hic file
//' @param chr1loc first chromosome location
//' @param chr2loc second chromosome location
//' @param unit BP (BasePair) or FRAG (FRAGment)
//' @param binsize The bin size. By default, for BP, this is one of <2500000, 1000000, 500000,
//'     250000, 100000, 50000, 25000, 10000, 5000> and for FRAG this is one of <500, 200,
//'     100, 50, 20, 5, 2, 1>.
//' @param matrix Type of matrix to output. Must be one of observed/oe/expected.
//'     observed is observed counts, oe is observed/expected counts, expected is expected counts.
//' @return Data.frame of a sparse matrix of data from hic file. x,y,counts
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame straw(std::string norm, std::string fname, std::string chr1loc, std::string chr2loc, std::string unit, int32_t binsize, std::string matrix = "observed")
{
  blockMap.clear();
  if (!(unit=="BP"||unit=="FRAG")) {
    stop("Norm specified incorrectly, must be one of <BP/FRAG>.\nUsage: juicebox-quick-dump <observed/oe> <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>");
  }

  ifstream fin(fname, fstream::in | fstream::binary);
  if (!fin) {
    stop("File %s cannot be opened for reading.", fname);
  }
  stringstream ss(chr1loc);
  string chr1, chr2, x, y;
  int32_t c1pos1=-100, c1pos2=-100, c2pos1=-100, c2pos2=-100;
  int32_t chr1len, chr2len;
  getline(ss, chr1, ':');
  if (getline(ss, x, ':') && getline(ss, y, ':')) {
    c1pos1 = stoi(x);
    c1pos2 = stoi(y);
  }
  stringstream ss1(chr2loc);
  getline(ss1, chr2, ':');
  if (getline(ss1, x, ':') && getline(ss1, y, ':')) {
    c2pos1 = stoi(x);
    c2pos2 = stoi(y);
  }
  int32_t chr1ind, chr2ind;
  int64_t master = readHeader(fin, chr1, chr2, c1pos1, c1pos2, c2pos1, c2pos2, chr1ind, chr2ind, chr1len, chr2len);

  int32_t c1=min(chr1ind,chr2ind);
  int32_t c2=max(chr1ind,chr2ind);
  int32_t origRegionIndices[4]; // as given by user
  int32_t regionIndices[4]; // used to find the blocks we need to access
  // reverse order if necessary
  if (chr1ind > chr2ind) {
    origRegionIndices[0] = c2pos1;
    origRegionIndices[1] = c2pos2;
    origRegionIndices[2] = c1pos1;
    origRegionIndices[3] = c1pos2;
    regionIndices[0] = c2pos1 / binsize;
    regionIndices[1] = c2pos2 / binsize;
    regionIndices[2] = c1pos1 / binsize;
    regionIndices[3] = c1pos2 / binsize;
  }
  else {
    origRegionIndices[0] = c1pos1;
    origRegionIndices[1] = c1pos2;
    origRegionIndices[2] = c2pos1;
    origRegionIndices[3] = c2pos2;
    regionIndices[0] = c1pos1 / binsize;
    regionIndices[1] = c1pos2 / binsize;
    regionIndices[2] = c2pos1 / binsize;
    regionIndices[3] = c2pos2 / binsize;
  }

  indexEntry c1NormEntry, c2NormEntry;
  int64_t myFilePos;
  vector<double> expectedValues;

  // readFooter will assign the above variables
  readFooter(fin, master, c1, c2, matrix, norm, unit, binsize, myFilePos, c1NormEntry, c2NormEntry, expectedValues);

  vector<double> c1Norm;
  vector<double> c2Norm;

  if (norm != "NONE") {
    c1Norm = readNormalizationVector(fin, c1NormEntry);
    c2Norm = readNormalizationVector(fin, c2NormEntry);
  }
  float sumCounts;
  int32_t blockBinCount, blockColumnCount;
  // readMatrix will assign sumCounts, blockBinCount, and blockColumnCount
  readMatrix(fin, myFilePos, unit, binsize, sumCounts, blockBinCount, blockColumnCount);
  double avgCount;
  if (c1 != c2) {
    int64_t nBins1 = chr1len / binsize;
    int64_t nBins2 = chr2len / binsize;
    avgCount = (sumCounts / nBins1) / nBins2;   // <= trying to avoid overflows
  }

  set<int32_t> blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount, c1==c2);

  // getBlockIndices
  vector<contactRecord> records;
  vector<int32_t> xActual_vec, yActual_vec;
  vector<float> counts_vec;
  for (set<int32_t>::iterator it=blockNumbers.begin(); it!=blockNumbers.end(); ++it) {
    // get contacts in this block
    records = readBlock(fin, *it);
    for (vector<contactRecord>::iterator it2=records.begin(); it2!=records.end(); ++it2) {
      contactRecord rec = *it2;

      int32_t xActual = rec.binX * binsize;
      int32_t yActual = rec.binY * binsize;
      float counts = rec.counts;
      if (norm != "NONE") {
        counts = counts / (c1Norm[rec.binX] * c2Norm[rec.binY]);
      }
      if (matrix == "oe") {
        if (c1 == c2) {
          counts = counts / expectedValues[min(expectedValues.size() - 1, (size_t)floor(abs(yActual - xActual) / binsize))];
        }
        else {
          counts = counts / avgCount;
        }
      }
      else if (matrix == "expected") {
        if (c1 == c2) {
          counts = expectedValues[min(expectedValues.size() - 1, (size_t)floor(abs(yActual - xActual) / binsize))];
        }
        else {
          counts = avgCount;
        }
      }
//      cout << xActual << " " << yActual << " " << counts << endl;
      if ((xActual >= origRegionIndices[0] && xActual <= origRegionIndices[1] &&
	   yActual >= origRegionIndices[2] && yActual <= origRegionIndices[3]) ||
	  // or check regions that overlap with lower left
	  ((c1==c2) && yActual >= origRegionIndices[0] && yActual <= origRegionIndices[1] && xActual >= origRegionIndices[2] && xActual <= origRegionIndices[3])) {
	  //printf("%d\t%d\t%.14g\n", xActual, yActual, counts);
    xActual_vec.push_back(xActual);
    yActual_vec.push_back(yActual);
    counts_vec.push_back(counts);
      }
    }
  }
  return Rcpp::DataFrame::create(Rcpp::Named("x") = xActual_vec, Rcpp::Named("y") = yActual_vec, Rcpp::Named("counts") = counts_vec);
}

//' Function for reading basepair resolutions from .hic file
//'
//' @param fname path to .hic file
//' @return Vector of basepair resolutions
//' @export
// [[Rcpp::export]]
NumericVector readHicBpResolutions(std::string fname)
{
  ifstream fin(fname, ios::in | ios::binary);
  if (!fin) {
    stop("File %s cannot be opened for reading.", fname);
  }

  if (!readMagicString(fin)) {
    fin.close();
    stop("Hi-C magic string is missing, does not appear to be a hic file.");
  }

  int version;
  fin.read((char*)&version, sizeof(int));
  if (version < 6) {
    fin.close();
    stop("Version %d no longer supported.", version);
  }
  long master;
  fin.read((char*)&master, sizeof(long));
  string genome;
  getline(fin, genome, '\0' );
  int nattributes;
  fin.read((char*)&nattributes, sizeof(int));
  // reading and ignoring attribute-value dictionary
  for (int i=0; i<nattributes; i++) {
    string key, value;
    getline(fin, key, '\0');
    getline(fin, value, '\0');
  }
  int nChrs;
  fin.read((char*)&nChrs, sizeof(int));
  // chromosome map for finding matrix
  for (int i=0; i<nChrs; i++) {
    string name;
    int length;
    getline(fin, name, '\0');
    fin.read((char*)&length, sizeof(int));
  }
  int nBpResolutions;
  fin.read((char*)&nBpResolutions, sizeof(int));
  NumericVector bpResolutions(nBpResolutions);
  for (int i=0; i<nBpResolutions; i++) {
    int resBP;
    fin.read((char*)&resBP, sizeof(int));
    bpResolutions[i] = resBP;
  }

  fin.close();

  return bpResolutions;
}

//' Function for reading chromosomes from .hic file
//'
//' @param fname path to .hic file
//' @return Data frame of chromosome names and lengths
//' @export
// [[Rcpp::export]]
DataFrame readHicChroms(std::string fname)
{
  ifstream fin(fname, ios::in | ios::binary);
  if (!fin) {
    stop("File %s cannot be opened for reading.", fname);
  }

  if (!readMagicString(fin)) {
    fin.close();
    stop("Hi-C magic string is missing, does not appear to be a hic file.");
  }

  int version;
  fin.read((char*)&version, sizeof(int));
  if (version < 6) {
    fin.close();
    stop("Version %d no longer supported.", version);
  }
  long master;
  fin.read((char*)&master, sizeof(long));
  string genome;
  getline(fin, genome, '\0' );
  int nattributes;
  fin.read((char*)&nattributes, sizeof(int));
  // reading and ignoring attribute-value dictionary
  for (int i=0; i<nattributes; i++) {
    string key, value;
    getline(fin, key, '\0');
    getline(fin, value, '\0');
  }
  int nChrs;
  fin.read((char*)&nChrs, sizeof(int));
  StringVector chrom_names(nChrs);
  NumericVector chrom_lengths(nChrs);
  for (int i=0; i<nChrs; i++) {
    string name;
    int length;
    getline(fin, name, '\0');
    fin.read((char*)&length, sizeof(int));
    chrom_names[i] = name;
    chrom_lengths[i] = length;
  }
  fin.close();

  return DataFrame::create(Named("name") = chrom_names , Named("length") = chrom_lengths);
}
