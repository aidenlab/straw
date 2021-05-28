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
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <set>
#include <vector>
#include <streambuf>
#include <curl/curl.h>
#include "zlib.h"
#include "straw.h"
using namespace std;

/*
  Straw: fast C++ implementation of dump. Not as fully featured as the
  Java version. Reads the .hic file, finds the appropriate matrix and slice
  of data, and outputs as text in sparse upper triangular format.

  Currently only supporting matrices.

  Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>
 */
// this is for creating a stream from a byte array for ease of use
struct membuf : std::streambuf {
    membuf(char *begin, char *end) {
        this->setg(begin, begin, end);
    }
};

// for holding data from URL call
struct MemoryStruct {
    char *memory;
    size_t size;
};

// callback for libcurl. data written to this buffer
static size_t
WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
    size_t realsize = size * nmemb;
    struct MemoryStruct *mem = (struct MemoryStruct *) userp;

    mem->memory = static_cast<char *>(realloc(mem->memory, mem->size + realsize + 1));
    if (mem->memory == nullptr) {
        /* out of memory! */
        printf("not enough memory (realloc returned NULL)\n");
        return 0;
    }

    std::memcpy(&(mem->memory[mem->size]), contents, realsize);
    mem->size += realsize;
    mem->memory[mem->size] = 0;

    return realsize;
}

// get a buffer that can be used as an input stream from the URL
char *getData(CURL *curl, long long position, long long chunksize) {
    std::ostringstream oss;
    struct MemoryStruct chunk;

    chunk.memory = static_cast<char *>(malloc(1));
    chunk.size = 0;    /* no data at this point */
    oss << position << "-" << position + chunksize;
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *) &chunk);
    curl_easy_setopt(curl, CURLOPT_RANGE, oss.str().c_str());
    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        fprintf(stderr, "curl_easy_perform() failed: %s\n",
                curl_easy_strerror(res));
    }
    //  printf("%lu bytes retrieved\n", (long)chunk.size);

    return chunk.memory;
}

// returns whether or not this is valid HiC file
bool readMagicString(istream &fin) {
    string str;
    getline(fin, str, '\0');
    return str[0] == 'H' && str[1] == 'I' && str[2] == 'C';
}

char readCharFromFile(istream &fin) {
    char tempChar;
    fin.read(&tempChar, sizeof(char));
    return tempChar;
}

short readShortFromFile(istream &fin) {
    short tempShort;
    fin.read((char *) &tempShort, sizeof(short));
    return tempShort;
}

int readIntFromFile(istream &fin) {
    int tempInt;
    fin.read((char *) &tempInt, sizeof(int));
    return tempInt;
}

long long readLongLongFromFile(istream &fin) {
    long long tempLong;
    fin.read((char *) &tempLong, sizeof(long));
    return tempLong;
}

float readFloatFromFile(istream &fin) {
    float tempFloat;
    fin.read((char *) &tempFloat, sizeof(float));
    return tempFloat;
}

double readDoubleFromFile(istream &fin) {
    double tempDouble;
    fin.read((char *) &tempDouble, sizeof(double));
    return tempDouble;
}

// reads the header, storing the positions of the normalization vectors and returning the masterIndexPosition pointer
map<string, chromosome> readHeader(istream &fin, long long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long long &nviPosition, long long &nviLength) {
    map<string, chromosome> chromosomeMap;
    if (!readMagicString(fin)) {
        cerr << "Hi-C magic string is missing, does not appear to be a hic file" << endl;
        masterIndexPosition = -1;
        return chromosomeMap;
    }

    fin.read((char *) &version, sizeof(int));
    if (version < 6) {
        cerr << "Version " << version << " no longer supported" << endl;
        masterIndexPosition = -1;
        return chromosomeMap;
    }
    fin.read((char *) &masterIndexPosition, sizeof(long));
    getline(fin, genomeID, '\0');

    if (version > 8) {
        nviPosition = readLongLongFromFile(fin);
        nviLength = readLongLongFromFile(fin);
    }

    int nattributes = readIntFromFile(fin);

    // reading and ignoring attribute-value dictionary
    for (int i = 0; i < nattributes; i++) {
        string key, value;
        getline(fin, key, '\0');
        getline(fin, value, '\0');
    }

    numChromosomes = readIntFromFile(fin);
    // chromosome map for finding matrixType
    for (int i = 0; i < numChromosomes; i++) {
        string name;
        long long length;
        getline(fin, name, '\0');
        if (version > 8) {
            fin.read((char *) &length, sizeof(long));
        } else {
            length = (long) readIntFromFile(fin);
        }

        chromosome chr;
        chr.index = i;
        chr.name = name;
        chr.length = length;
        chromosomeMap[name] = chr;
    }
    return chromosomeMap;
}

// reads the footer from the master pointer location. takes in the chromosomes,
// norm, unit (BP or FRAG) and resolution or binsize, and sets the file
// position of the matrix and the normalization vectors for those chromosomes
// at the given normalization and resolution
bool readFooter(istream &fin, long long master, int version, int c1, int c2, const string &matrixType, const string &norm,
                const string &unit, int resolution, long long &myFilePos, indexEntry &c1NormEntry, indexEntry &c2NormEntry,
                vector<double> &expectedValues) {
    if (version > 8) {
        long long nBytes = readLongLongFromFile(fin);
    } else {
        int nBytes = readIntFromFile(fin);
    }

    stringstream ss;
    ss << c1 << "_" << c2;
    string key = ss.str();

    int nEntries = readIntFromFile(fin);
    bool found = false;
    for (int i = 0; i < nEntries; i++) {
        string str;
        getline(fin, str, '\0');
        long long fpos = readLongLongFromFile(fin);
        int sizeinbytes = readIntFromFile(fin);
        if (str == key) {
            myFilePos = fpos;
            found = true;
        }
    }
    if (!found) {
        cerr << "File doesn't have the given chr_chr map " << key << endl;
        return false;
    }

    if ((matrixType == "observed" && norm == "NONE") || (matrixType == "oe" && norm == "NONE" && c1 != c2))
        return true; // no need to read norm vector index

    // read in and ignore expected value maps; don't store; reading these to
    // get to norm vector index
    int nExpectedValues = readIntFromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string unit0;
        getline(fin, unit0, '\0'); //unit
        int binSize = readIntFromFile(fin);

        long long nValues;
        if (version > 8) {
            fin.read((char *) &nValues, sizeof(long));
        } else {
            nValues = (long) readIntFromFile(fin);
        }

        bool store = c1 == c2 && matrixType == "oe" && norm == "NONE" && unit0 == unit && binSize == resolution;

        if (version > 8) {
            for (long long j = 0; j < nValues; j++) {
                double v = readFloatFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }
        } else {
            for (int j = 0; j < nValues; j++) {
                double v = readDoubleFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }
        }

        int nNormalizationFactors = readIntFromFile(fin);
        for (int j = 0; j < nNormalizationFactors; j++) {
            int chrIdx = readIntFromFile(fin);
            double v;
            if (version > 8) {
                v = readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (store && chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
    }

    if (c1 == c2 && matrixType == "oe" && norm == "NONE") {
        if (expectedValues.empty()) {
            cerr << "File did not contain expected values vectors at " << resolution << " " << unit << endl;
            return false;
        }
        return true;
    }

    nExpectedValues = readIntFromFile(fin);
    for (int i = 0; i < nExpectedValues; i++) {
        string type, unit0;
        getline(fin, type, '\0'); //typeString
        getline(fin, unit0, '\0'); //unit
        int binSize = readIntFromFile(fin);

        long long nValues;
        if (version > 8) {
            fin.read((char *) &nValues, sizeof(long));
        } else {
            nValues = (long) readIntFromFile(fin);
        }
        bool store = c1 == c2 && matrixType == "oe" && type == norm && unit0 == unit && binSize == resolution;

        if (version > 8) {
            for (long long j = 0; j < nValues; j++) {
                double v = readFloatFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }
        } else {
            for (int j = 0; j < nValues; j++) {
                double v = readDoubleFromFile(fin);
                if (store) {
                    expectedValues.push_back(v);
                }
            }

        }

        int nNormalizationFactors = readIntFromFile(fin);
        for (int j = 0; j < nNormalizationFactors; j++) {
            int chrIdx = readIntFromFile(fin);
            double v;
            if (version > 8) {
                v = (double) readFloatFromFile(fin);
            } else {
                v = readDoubleFromFile(fin);
            }
            if (store && chrIdx == c1) {
                for (double &expectedValue : expectedValues) {
                    expectedValue = expectedValue / v;
                }
            }
        }
    }

    if (c1 == c2 && matrixType == "oe" && norm != "NONE") {
        if (expectedValues.empty()) {
            cerr << "File did not contain normalized expected values vectors at " << resolution << " " << unit << endl;
            return false;
        }
    }

    // Index of normalization vectors
    nEntries = readIntFromFile(fin);
    bool found1 = false;
    bool found2 = false;
    for (int i = 0; i < nEntries; i++) {
        string normtype;
        getline(fin, normtype, '\0'); //normalization type
        int chrIdx = readIntFromFile(fin);
        string unit1;
        getline(fin, unit1, '\0'); //unit
        int resolution1 = readIntFromFile(fin);
        long long filePosition = readLongLongFromFile(fin);
        long long sizeInBytes;
        if (version > 8) {
            fin.read((char *) &sizeInBytes, sizeof(long));
        } else {
            sizeInBytes = (long) readIntFromFile(fin);
        }

        if (chrIdx == c1 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c1NormEntry.position = filePosition;
            c1NormEntry.size = sizeInBytes;
            found1 = true;
        }
        if (chrIdx == c2 && normtype == norm && unit1 == unit && resolution1 == resolution) {
            c2NormEntry.position = filePosition;
            c2NormEntry.size = sizeInBytes;
            found2 = true;
        }
    }
    if (!found1 || !found2) {
        cerr << "File did not contain " << norm << " normalization vectors for one or both chromosomes at "
             << resolution << " " << unit << endl;
    }
    return true;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map<int, indexEntry> readMatrixZoomData(istream &fin, const string &myunit, int mybinsize, float &mySumCounts,
                                        int &myBlockBinCount, int &myBlockColumnCount, bool &found) {

    map<int, indexEntry> blockMap;
    string unit;
    getline(fin, unit, '\0'); // unit
    readIntFromFile(fin); // Old "zoom" index -- not used
    float sumCounts = readFloatFromFile(fin); // sumCounts
    readFloatFromFile(fin); // occupiedCellCount
    readFloatFromFile(fin); // stdDev
    readFloatFromFile(fin); // percent95
    int binSize = readIntFromFile(fin);
    int blockBinCount = readIntFromFile(fin);
    int blockColumnCount = readIntFromFile(fin);

    found = false;
    if (myunit == unit && mybinsize == binSize) {
        mySumCounts = sumCounts;
        myBlockBinCount = blockBinCount;
        myBlockColumnCount = blockColumnCount;
        found = true;
    }

    int nBlocks = readIntFromFile(fin);

    for (int b = 0; b < nBlocks; b++) {
        int blockNumber = readIntFromFile(fin);
        long long filePosition = readLongLongFromFile(fin);
        int blockSizeInBytes = readIntFromFile(fin);
        indexEntry entry = indexEntry();
        entry.size = (long) blockSizeInBytes;
        entry.position = filePosition;
        if (found) blockMap[blockNumber] = entry;
    }
    return blockMap;
}

// reads the raw binned contact matrix at specified resolution, setting the block bin count and block column count
map<int, indexEntry> readMatrixZoomDataHttp(CURL *curl, long long &myFilePosition, const string &myunit, int mybinsize,
                                            float &mySumCounts, int &myBlockBinCount, int &myBlockColumnCount,
                                            bool &found) {

    map<int, indexEntry> blockMap;
    char *buffer;
    int header_size = 5 * sizeof(int) + 4 * sizeof(float);
    char *first;
    first = getData(curl, myFilePosition, 1);
    if (first[0] == 'B') {
        header_size += 3;
    } else if (first[0] == 'F') {
        header_size += 5;
    } else {
        cerr << "Unit not understood" << endl;
        return blockMap;
    }
    buffer = getData(curl, myFilePosition, header_size);
    membuf sbuf(buffer, buffer + header_size);
    istream fin(&sbuf);

    string unit;
    getline(fin, unit, '\0'); // unit
    readIntFromFile(fin); // Old "zoom" index -- not used
    float sumCounts = readFloatFromFile(fin); // sumCounts
    readFloatFromFile(fin); // occupiedCellCount
    readFloatFromFile(fin); // stdDev
    readFloatFromFile(fin); // percent95
    int binSize = readIntFromFile(fin);
    int blockBinCount = readIntFromFile(fin);
    int blockColumnCount = readIntFromFile(fin);

    found = false;
    if (myunit == unit && mybinsize == binSize) {
        mySumCounts = sumCounts;
        myBlockBinCount = blockBinCount;
        myBlockColumnCount = blockColumnCount;
        found = true;
    }

    int nBlocks = readIntFromFile(fin);

    if (found) {
        int chunkSize = nBlocks * (sizeof(int) + sizeof(long) + sizeof(int));
        buffer = getData(curl, myFilePosition + header_size, chunkSize);
        membuf sbuf2(buffer, buffer + chunkSize);
        istream fin2(&sbuf2);
        for (int b = 0; b < nBlocks; b++) {
            int blockNumber = readIntFromFile(fin2);
            long long filePosition = readLongLongFromFile(fin2);
            int blockSizeInBytes = readIntFromFile(fin2);
            indexEntry entry = indexEntry();
            entry.size = (long) blockSizeInBytes;
            entry.position = filePosition;
            blockMap[blockNumber] = entry;
        }
    } else {
        myFilePosition = myFilePosition + header_size + (nBlocks * (sizeof(int) + sizeof(long) + sizeof(int)));
    }
    delete buffer;
    return blockMap;
}

// goes to the specified file pointer in http and finds the raw contact matrixType at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map<int, indexEntry> readMatrixHttp(CURL *curl, long long myFilePosition, const string &unit, int resolution,
                                    float &mySumCounts, int &myBlockBinCount, int &myBlockColumnCount) {
    char *buffer;
    int size = sizeof(int) * 3;
    buffer = getData(curl, myFilePosition, size);
    membuf sbuf(buffer, buffer + size);
    istream bufin(&sbuf);

    int c1 = readIntFromFile(bufin);
    int c2 = readIntFromFile(bufin);
    int nRes = readIntFromFile(bufin);
    int i = 0;
    bool found = false;
    myFilePosition = myFilePosition + size;
    delete buffer;
    map<int, indexEntry> blockMap;

    while (i < nRes && !found) {
        // myFilePosition gets updated within call
        blockMap = readMatrixZoomDataHttp(curl, myFilePosition, unit, resolution, mySumCounts, myBlockBinCount, myBlockColumnCount,
                                          found);
        i++;
    }
    if (!found) {
        cerr << "Error finding block data" << endl;
    }
    return blockMap;
}

// goes to the specified file pointer and finds the raw contact matrixType at specified resolution, calling readMatrixZoomData.
// sets blockbincount and blockcolumncount
map<int, indexEntry> readMatrix(istream &fin, long long myFilePosition, const string &unit, int resolution,
                                float &mySumCounts, int &myBlockBinCount, int &myBlockColumnCount) {
    map<int, indexEntry> blockMap;

    fin.seekg(myFilePosition, ios::beg);
    int c1 = readIntFromFile(fin);
    int c2 = readIntFromFile(fin);
    int nRes = readIntFromFile(fin);
    int i = 0;
    bool found = false;
    while (i < nRes && !found) {
        blockMap = readMatrixZoomData(fin, unit, resolution, mySumCounts, myBlockBinCount, myBlockColumnCount, found);
        i++;
    }
    if (!found) {
        cerr << "Error finding block data" << endl;
    }
    return blockMap;
}

// gets the blocks that need to be read for this slice of the data.  needs blockbincount, blockcolumncount, and whether
// or not this is intrachromosomal.
set<int> getBlockNumbersForRegionFromBinPosition(const long long *regionIndices, int blockBinCount, int blockColumnCount,
                                                 bool intra) {
    int col1 = static_cast<int>(regionIndices[0] / blockBinCount);
    int col2 = static_cast<int>((regionIndices[1] + 1) / blockBinCount);
    int row1 = static_cast<int>(regionIndices[2] / blockBinCount);
    int row2 = static_cast<int>((regionIndices[3] + 1) / blockBinCount);

    set<int> blocksSet;
    // first check the upper triangular matrixType
    for (int r = row1; r <= row2; r++) {
        for (int c = col1; c <= col2; c++) {
            int blockNumber = r * blockColumnCount + c;
            blocksSet.insert(blockNumber);
        }
    }
    // check region part that overlaps with lower left triangle but only if intrachromosomal
    if (intra) {
        for (int r = col1; r <= col2; r++) {
            for (int c = row1; c <= row2; c++) {
                int blockNumber = r * blockColumnCount + c;
                blocksSet.insert(blockNumber);
            }
        }
    }
    return blocksSet;
}

set<int> getBlockNumbersForRegionFromBinPositionV9Intra(long long *regionIndices, int blockBinCount, int blockColumnCount) {
    // regionIndices is binX1 binX2 binY1 binY2
    set<int> blocksSet;
    int translatedLowerPAD = static_cast<int>((regionIndices[0] + regionIndices[2]) / 2 / blockBinCount);
    int translatedHigherPAD = static_cast<int>((regionIndices[1] + regionIndices[3]) / 2 / blockBinCount + 1);
    int translatedNearerDepth = static_cast<int>(log2(
            1 + abs(regionIndices[0] - regionIndices[3]) / sqrt(2) / blockBinCount));
    int translatedFurtherDepth = static_cast<int>(log2(
            1 + abs(regionIndices[1] - regionIndices[2]) / sqrt(2) / blockBinCount));

    // because code above assume above diagonal; but we could be below diagonal
    int nearerDepth = min(translatedNearerDepth, translatedFurtherDepth);
    if ((regionIndices[0] > regionIndices[3] && regionIndices[1] < regionIndices[2]) ||
        (regionIndices[1] > regionIndices[2] && regionIndices[0] < regionIndices[3])) {
        nearerDepth = 0;
    }
    int furtherDepth = max(translatedNearerDepth, translatedFurtherDepth) + 1; // +1; integer divide rounds down

    for (int depth = nearerDepth; depth <= furtherDepth; depth++) {
        for (int pad = translatedLowerPAD; pad <= translatedHigherPAD; pad++) {
            int blockNumber = depth * blockColumnCount + pad;
            blocksSet.insert(blockNumber);
        }
    }

    return blocksSet;
}

void appendRecord(vector<contactRecord> &vector, int index, int binX, int binY, float counts) {
    contactRecord record = contactRecord();
    record.binX = binX;
    record.binY = binY;
    record.counts = counts;
    vector[index] = record;
}

// this is the meat of reading the data.  takes in the block number and returns the set of contact records corresponding to
// that block.  the block data is compressed and must be decompressed using the zlib library functions
vector<contactRecord> readBlock(istream &fin, CURL *curl, bool isHttp, indexEntry idx, int version) {
    if (idx.size <= 0) {
        vector<contactRecord> v;
        return v;
    }
    char *compressedBytes = new char[idx.size];
    char *uncompressedBytes = new char[idx.size * 10]; //biggest seen so far is 3

    if (isHttp) {
        compressedBytes = getData(curl, idx.position, idx.size);
    } else {
        fin.seekg(idx.position, ios::beg);
        fin.read(compressedBytes, idx.size);
    }
    // Decompress the block
    // zlib struct
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = static_cast<uInt>(idx.size); // size of input
    infstream.next_in = (Bytef *) compressedBytes; // input char array
    infstream.avail_out = static_cast<uInt>(idx.size * 10); // size of output
    infstream.next_out = (Bytef *) uncompressedBytes; // output char array
    // the actual decompression work.
    inflateInit(&infstream);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    int uncompressedSize = static_cast<int>(infstream.total_out);

    // create stream from buffer for ease of use
    membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
    istream bufferin(&sbuf);
    unsigned long long nRecords = static_cast<unsigned long long>(readIntFromFile(bufferin));
    vector<contactRecord> v(nRecords);
    // different versions have different specific formats
    if (version < 7) {
        for (int i = 0; i < nRecords; i++) {
            int binX = readIntFromFile(bufferin);
            int binY = readIntFromFile(bufferin);
            float counts = readFloatFromFile(bufferin);
            appendRecord(v, i, binX, binY, counts);
        }
    } else {
        int binXOffset = readIntFromFile(bufferin);
        int binYOffset = readIntFromFile(bufferin);
        bool useShort = readCharFromFile(bufferin) == 0; // yes this is opposite of usual

        bool useShortBinX = true;
        bool useShortBinY = true;
        if (version > 8) {
            useShortBinX = readCharFromFile(bufferin) == 0;
            useShortBinY = readCharFromFile(bufferin) == 0;
        }

        char type = readCharFromFile(bufferin);
        int index = 0;
        if (type == 1) {
            if (useShortBinX && useShortBinY) {
                short rowCount = readShortFromFile(bufferin);
                for (short i = 0; i < rowCount; i++) {
                    int binY = binYOffset + readShortFromFile(bufferin);
                    short colCount = readShortFromFile(bufferin);
                    for (short j = 0; j < colCount; j++) {
                        int binX = binXOffset + readShortFromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readShortFromFile(bufferin);
                        } else {
                            bufferin.read((char *) &counts, sizeof(float));
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (useShortBinX && !useShortBinY) {
                int rowCount = readIntFromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int binY = binYOffset + readIntFromFile(bufferin);
                    short colCount = readShortFromFile(bufferin);
                    for (short j = 0; j < colCount; j++) {
                        int binX = binXOffset + readShortFromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readShortFromFile(bufferin);
                        } else {
                            bufferin.read((char *) &counts, sizeof(float));
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else if (!useShortBinX && useShortBinY) {
                short rowCount = readShortFromFile(bufferin);
                for (short i = 0; i < rowCount; i++) {
                    int binY = binYOffset + readShortFromFile(bufferin);
                    int colCount = readIntFromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int binX = binXOffset + readIntFromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readShortFromFile(bufferin);
                        } else {
                            bufferin.read((char *) &counts, sizeof(float));
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            } else {
                int rowCount = readIntFromFile(bufferin);
                for (int i = 0; i < rowCount; i++) {
                    int binY = binYOffset + readIntFromFile(bufferin);
                    int colCount = readIntFromFile(bufferin);
                    for (int j = 0; j < colCount; j++) {
                        int binX = binXOffset + readIntFromFile(bufferin);
                        float counts;
                        if (useShort) {
                            counts = readShortFromFile(bufferin);
                        } else {
                            bufferin.read((char *) &counts, sizeof(float));
                        }
                        appendRecord(v, index++, binX, binY, counts);
                    }
                }
            }
        } else if (type == 2) {
            int nPts = readIntFromFile(bufferin);
            short w = readShortFromFile(bufferin);

            for (int i = 0; i < nPts; i++) {
                //int idx = (p.y - binOffset2) * w + (p.x - binOffset1);
                int row = i / w;
                int col = i - row * w;
                int bin1 = binXOffset + col;
                int bin2 = binYOffset + row;

                float counts;
                if (useShort) {
                    short c = readShortFromFile(bufferin);
                    if (c != -32768) {
                        appendRecord(v, index++, bin1, bin2, c);
                    }
                } else {
                    bufferin.read((char *) &counts, sizeof(float));
                    if (!isnan(counts)) {
                        appendRecord(v, index++, bin1, bin2, counts);
                    }
                }
            }
        }
    }
    delete[] compressedBytes;
    delete[] uncompressedBytes; // don't forget to delete your heap arrays in C++!
    return v;
}

// reads the normalization vector from the file at the specified location
vector<double> readNormalizationVector(istream &bufferin, int version) {
    long long nValues;
    if (version > 8) {
        bufferin.read((char *) &nValues, sizeof(long));
    } else {
        nValues = (long) readIntFromFile(bufferin);
    }

    unsigned long long numValues = static_cast<unsigned long long>(nValues);
    vector<double> values(numValues);

    if (version > 8) {
        for (long long i = 0; i < nValues; i++) {
            values[i] = (double) readFloatFromFile(bufferin);
        }
    } else {
        for (int i = 0; i < nValues; i++) {
            values[i] = readDoubleFromFile(bufferin);
        }
    }

    return values;
}

class FileReader {
public:
    string prefix = "http"; // HTTP code
    ifstream fin;
    CURL *curl;
    bool isHttp = false;

    static CURL *initCURL(const char *url) {
        CURL *curl = curl_easy_init();
        if (curl) {
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
            curl_easy_setopt(curl, CURLOPT_URL, url);
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw");
        }
        return curl;
    }

    explicit FileReader(const string &fname) {

        // read header into buffer; 100K should be sufficient
        if (std::strncmp(fname.c_str(), prefix.c_str(), prefix.size()) == 0) {
            isHttp = true;
            curl = initCURL(fname.c_str());
            if (!curl) {
                cerr << "URL " << fname << " cannot be opened for reading" << endl;
                exit(1);
            }
        } else {
            fin.open(fname, fstream::in | fstream::binary);
            if (!fin) {
                cerr << "File " << fname << " cannot be opened for reading" << endl;
                exit(2);
            }
        }
    }

    void close(){
        if(isHttp){
            curl_easy_cleanup(curl);
        } else {
            fin.close();
        }
    }
};

class HiCFile {
public:
    string prefix = "http"; // HTTP code
    bool isHttp = false;
    ifstream fin;
    CURL *curl;
    long long master = 0LL;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int numChromosomes = 0;
    int version = 0;
    long long nviPosition = 0;
    long long nviLength = 0;
    static long long totalFileSize;

    static size_t hdf(char *b, size_t size, size_t nitems, void *userdata) {
        size_t numbytes = size * nitems;
        b[numbytes + 1] = '\0';
        string s(b);
        int found = static_cast<int>(s.find("Content-Range"));
        if (found != string::npos) {
            int found2 = static_cast<int>(s.find("/"));
            //Content-Range: bytes 0-100000/891471462
            if (found2 != string::npos) {
                string total = s.substr(found2 + 1);
                totalFileSize = stol(total);
            }
        }

        return numbytes;
    }

    static CURL *initCURL(const char *url) {
        CURL *curl = curl_easy_init();
        if (curl) {
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteMemoryCallback);
            curl_easy_setopt(curl, CURLOPT_URL, url);
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, hdf);
            curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw");
        }
        return curl;
    }

    explicit HiCFile(const string &fname) {

        // read header into buffer; 100K should be sufficient
        if (std::strncmp(fname.c_str(), prefix.c_str(), prefix.size()) == 0) {
            isHttp = true;
            char *buffer;
            curl = initCURL(fname.c_str());
            if (curl) {
                buffer = getData(curl, 0, 100000);
            } else {
                cerr << "URL " << fname << " cannot be opened for reading" << endl;
                exit(1);
            }
            membuf sbuf(buffer, buffer + 100000);
            istream bufin(&sbuf);
            chromosomeMap = readHeader(bufin, master, genomeID, numChromosomes,
                                       version, nviPosition, nviLength);
            delete buffer;
        } else {
            fin.open(fname, fstream::in | fstream::binary);
            if (!fin) {
                cerr << "File " << fname << " cannot be opened for reading" << endl;
                exit(2);
            }
            chromosomeMap = readHeader(fin, master, genomeID, numChromosomes,
                                       version, nviPosition, nviLength);
        }
    }

    void close(){
        if(isHttp){
            curl_easy_cleanup(curl);
        } else {
            fin.close();
        }
    }

    vector<double> readNormalizationVectorFromFooter(indexEntry cNormEntry) {
        char *buffer;
        if (isHttp) {
            buffer = getData(curl, cNormEntry.position, cNormEntry.size);
        } else {
            buffer = new char[cNormEntry.size];
            fin.seekg(cNormEntry.position, ios::beg);
            fin.read(buffer, cNormEntry.size);
        }
        membuf sbuf3(buffer, buffer + cNormEntry.size);
        istream bufferin(&sbuf3);
        vector<double> cNorm = readNormalizationVector(bufferin, version);
        delete buffer;
        return cNorm;
    }
};

long long HiCFile::totalFileSize = 0LL;

class MatrixZoomData {
public:
    indexEntry c1NormEntry, c2NormEntry;
    long long myFilePos = 0LL;
    vector<double> expectedValues;
    bool foundFooter = false;
    vector<double> c1Norm;
    vector<double> c2Norm;
    int c1 = 0;
    int c2 = 0;
    string matrixType;
    string norm;
    string unit;
    int resolution = 0;
    int numBins1 = 0;
    int numBins2 = 0;

    MatrixZoomData(HiCFile *hiCFile, const chromosome &chrom1, const chromosome &chrom2, const string &matrixType,
                   const string &norm, const string &unit, int resolution) {

        int c01 = chrom1.index;
        int c02 = chrom2.index;
        if (c01 <= c02) { // default is ok
            this->c1 = c01;
            this->c2 = c02;
            this->numBins1 = static_cast<int>(chrom1.length / resolution);
            this->numBins2 = static_cast<int>(chrom2.length / resolution);
        } else { // flip
            this->c1 = c02;
            this->c2 = c01;
            this->numBins1 = static_cast<int>(chrom2.length / resolution);
            this->numBins2 = static_cast<int>(chrom1.length / resolution);
        }

        this->matrixType = matrixType;
        this->norm = norm;
        this->unit = unit;
        this->resolution = resolution;

        if (hiCFile->isHttp) {
            char *buffer2;
            long long bytes_to_read = hiCFile->totalFileSize - hiCFile->master;
            buffer2 = getData(hiCFile->curl, hiCFile->master, bytes_to_read);
            membuf sbuf2(buffer2, buffer2 + bytes_to_read);
            istream bufin2(&sbuf2);
            foundFooter = readFooter(bufin2, hiCFile->master, hiCFile->version, c1, c2, matrixType, norm, unit,
                                     resolution,
                                     myFilePos,
                                     c1NormEntry, c2NormEntry, expectedValues);
            delete buffer2;
        } else {
            hiCFile->fin.seekg(hiCFile->master, ios::beg);
            foundFooter = readFooter(hiCFile->fin, hiCFile->master, hiCFile->version, c1, c2, matrixType, norm,
                                     unit,
                                     resolution, myFilePos,
                                     c1NormEntry, c2NormEntry, expectedValues);
        }

        if (!foundFooter) {
            return;
        }

        if (norm != "NONE") {
            c1Norm = hiCFile->readNormalizationVectorFromFooter(c1NormEntry);
            if (c1 == c2) {
                c2Norm = c1Norm;
            } else {
                c2Norm = hiCFile->readNormalizationVectorFromFooter(c2NormEntry);
            }
        }
    }
};

MatrixZoomData *
getMatrixZoomData(HiCFile *hiCFile, const string &chr1, const string &chr2, string matrixType, string norm,
                  string unit, int resolution) {

    chromosome chrom1 = hiCFile->chromosomeMap[chr1];
    chromosome chrom2 = hiCFile->chromosomeMap[chr2];
    return new MatrixZoomData(hiCFile, chrom1, chrom2, std::move(matrixType), std::move(norm), std::move(unit),
                              resolution);
}

void parsePositions(const string &chrLoc, string &chrom, long long &pos1, long long &pos2, map<string, chromosome> map) {
    string x, y;
    stringstream ss(chrLoc);
    getline(ss, chrom, ':');
    if (map.count(chrom) == 0) {
        cerr << chrom << " not found in the file." << endl;
        exit(6);
    }

    if (getline(ss, x, ':') && getline(ss, y, ':')) {
        pos1 = stol(x);
        pos2 = stol(y);
    } else {
        pos1 = 0LL;
        pos2 = map[chrom].length;
    }
}

class BlocksRecords {
public:
    float sumCounts;
    int blockBinCount, blockColumnCount;
    map<int, indexEntry> blockMap;
    double avgCount;
    bool isIntra;

    BlocksRecords(FileReader *fileReader, const footerInfo &footer) {

        isIntra = footer.c1 == footer.c2;

        if (fileReader->isHttp) {
            // readMatrix will assign blockBinCount and blockColumnCount
            blockMap = readMatrixHttp(fileReader->curl, footer.myFilePos, footer.unit, footer.resolution, sumCounts,
                                      blockBinCount,
                                      blockColumnCount);
        } else {
            // readMatrix will assign blockBinCount and blockColumnCount
            blockMap = readMatrix(fileReader->fin, footer.myFilePos, footer.unit, footer.resolution, sumCounts,
                                  blockBinCount,
                                  blockColumnCount);
        }

        if (!isIntra) {
            avgCount = (sumCounts / footer.numBins1) / footer.numBins2;   // <= trying to avoid overflows
        }
    }

    set<int> getBlockNumbers(int version, bool isIntra, long long *regionIndices, int blockBinCount, int blockColumnCount) {
        set<int> blockNumbers;
        if (version > 8 && isIntra) {
            blockNumbers = getBlockNumbersForRegionFromBinPositionV9Intra(regionIndices, blockBinCount,
                                                                          blockColumnCount);
        } else {
            blockNumbers = getBlockNumbersForRegionFromBinPosition(regionIndices, blockBinCount, blockColumnCount,
                                                                   isIntra);
        }
        return blockNumbers;
    }

    vector<contactRecord>
    getRecords(FileReader *fileReader, long long regionIndices[4],
               const long long origRegionIndices[4], const footerInfo &footer) {

        set<int> blockNumbers = getBlockNumbers(footer.version, isIntra, regionIndices, blockBinCount,
                                                blockColumnCount);

        vector<contactRecord> records;
        for (int blockNumber : blockNumbers) {
            // get contacts in this block
            //cout << *it << " -- " << blockMap.size() << endl;
            //cout << blockMap[*it].size << " " <<  blockMap[*it].position << endl;
            vector<contactRecord> tmp_records = readBlock(fileReader->fin, fileReader->curl, fileReader->isHttp,
                                                          blockMap[blockNumber], footer.version);
            for (contactRecord rec : tmp_records) {
                long long x = rec.binX * footer.resolution;
                long long y = rec.binY * footer.resolution;

                if ((x >= origRegionIndices[0] && x <= origRegionIndices[1] &&
                     y >= origRegionIndices[2] && y <= origRegionIndices[3]) ||
                    // or check regions that overlap with lower left
                    (isIntra && y >= origRegionIndices[0] && y <= origRegionIndices[1] && x >= origRegionIndices[2] &&
                     x <= origRegionIndices[3])) {

                    float c = rec.counts;
                    if (footer.norm != "NONE") {
                        c = static_cast<float>(c / (footer.c1Norm[rec.binX] * footer.c2Norm[rec.binY]));
                    }
                    if (footer.matrixType == "oe") {
                        if (isIntra) {
                            c = static_cast<float>(c / footer.expectedValues[min(footer.expectedValues.size() - 1,
                                                                                 (size_t) floor(abs(y - x) /
                                                                                                footer.resolution))]);
                        } else {
                            c = static_cast<float>(c / avgCount);
                        }
                    }

                    contactRecord record = contactRecord();
                    record.binX = static_cast<int>(x);
                    record.binY = static_cast<int>(y);
                    record.counts = c;
                    records.push_back(record);
                }
            }
        }
        return records;
    }
};

vector<contactRecord> getBlockRecords(FileReader *fileReader, long long origRegionIndices[4], const footerInfo &footer) {
    if (!footer.foundFooter) {
        vector<contactRecord> v;
        return v;
    }

    long long regionIndices[4]; // used to find the blocks we need to access
    regionIndices[0] = origRegionIndices[0] / footer.resolution;
    regionIndices[1] = origRegionIndices[1] / footer.resolution;
    regionIndices[2] = origRegionIndices[2] / footer.resolution;
    regionIndices[3] = origRegionIndices[3] / footer.resolution;

    BlocksRecords *blocksRecords = new BlocksRecords(fileReader, footer);
    return blocksRecords->getRecords(fileReader, regionIndices, origRegionIndices, footer);
}

footerInfo getNormalizationInfoForRegion(string fname, string chr1, string chr2,
                                         const string &matrixType, const string &norm,
                                         const string &unit, int binsize) {

    HiCFile *hiCFile = new HiCFile(std::move(fname));
    MatrixZoomData *mzd = getMatrixZoomData(hiCFile, chr1, chr2, std::move(matrixType), std::move(norm), unit,
                                            binsize);
    footerInfo footer = footerInfo();
    footer.resolution = mzd->resolution;
    footer.foundFooter = mzd->foundFooter;
    footer.version = hiCFile->version;
    footer.c1 = mzd->c1;
    footer.c2 = mzd->c2;
    footer.numBins1 = mzd->numBins1;
    footer.numBins2 = mzd->numBins2;
    footer.myFilePos = mzd->myFilePos;
    footer.unit = mzd->unit;
    footer.norm = mzd->norm;
    footer.matrixType = mzd->matrixType;
    footer.c1Norm = mzd->c1Norm;
    footer.c2Norm = mzd->c2Norm;
    footer.expectedValues = mzd->expectedValues;
    hiCFile->close();
    return footer;
}

vector<contactRecord>
getBlockRecordsWithNormalization(string fname,
                                 long long c1pos1, long long c1pos2, long long c2pos1, long long c2pos2,
                                 int resolution, bool foundFooter, int version, int c1, int c2,
                                 int numBins1, int numBins2, long long myFilePos, string unit, string norm,
                                 string matrixType,
                                 vector<double> c1Norm, vector<double> c2Norm, vector<double> expectedValues) {
    long long origRegionIndices[4]; // as given by user
    origRegionIndices[0] = c1pos1;
    origRegionIndices[1] = c1pos2;
    origRegionIndices[2] = c2pos1;
    origRegionIndices[3] = c2pos2;

    FileReader *fileReader = new FileReader(std::move(fname));
    footerInfo footer = footerInfo();
    footer.resolution = resolution;
    footer.foundFooter = foundFooter;
    footer.version = version;
    footer.c1 = c1;
    footer.c2 = c2;
    footer.numBins1 = numBins1;
    footer.numBins2 = numBins2;
    footer.myFilePos = myFilePos;
    footer.unit = unit;
    footer.norm = norm;
    footer.matrixType = matrixType;
    footer.c1Norm = c1Norm;
    footer.c2Norm = c2Norm;
    footer.expectedValues = expectedValues;
    vector<contactRecord> v = getBlockRecords(fileReader, origRegionIndices, footer);
    fileReader->close();
    return v;
}

vector<contactRecord>
straw(string matrixType, string norm, string fname, string chr1loc, string chr2loc, const string &unit, int binsize) {
    if (!(unit == "BP" || unit == "FRAG")) {
        cerr << "Norm specified incorrectly, must be one of <BP/FRAG>" << endl;
        cerr << "Usage: straw <NONE/VC/VC_SQRT/KR> <hicFile(s)> <chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>"
             << endl;
        vector<contactRecord> v;
        return v;
    }

    HiCFile *hiCFile = new HiCFile(std::move(fname));

    string chr1, chr2;
    long long c1pos1 = -100, c1pos2 = -100, c2pos1 = -100, c2pos2 = -100;
    parsePositions(std::move(chr1loc), chr1, c1pos1, c1pos2, hiCFile->chromosomeMap);
    parsePositions(std::move(chr2loc), chr2, c2pos1, c2pos2, hiCFile->chromosomeMap);

    // from header have size of chromosomes, set region to read

    long long origRegionIndices[4]; // as given by user
    // reverse order if necessary
    if (hiCFile->chromosomeMap[chr1].index > hiCFile->chromosomeMap[chr2].index) {
        origRegionIndices[0] = c2pos1;
        origRegionIndices[1] = c2pos2;
        origRegionIndices[2] = c1pos1;
        origRegionIndices[3] = c1pos2;
    } else {
        origRegionIndices[0] = c1pos1;
        origRegionIndices[1] = c1pos2;
        origRegionIndices[2] = c2pos1;
        origRegionIndices[3] = c2pos2;
    }
    hiCFile->close();

    footerInfo footer = getNormalizationInfoForRegion(fname, chr1, chr2, matrixType, norm, unit, binsize);

    return getBlockRecordsWithNormalization(fname,
                                            origRegionIndices[0], origRegionIndices[1],
                                            origRegionIndices[2], origRegionIndices[3],
                                            footer.resolution, footer.foundFooter, footer.version,
                                            footer.c1, footer.c2, footer.numBins1, footer.numBins2,
                                            footer.myFilePos, footer.unit, footer.norm, footer.matrixType,
                                            footer.c1Norm, footer.c2Norm, footer.expectedValues);
}

vector<chromosome> getChromosomes(string fname){
    HiCFile *hiCFile = new HiCFile(std::move(fname));
    vector<chromosome> chromosomes;
    std::map<std::string, chromosome>::iterator iter = hiCFile->chromosomeMap.begin();
    while (iter != hiCFile->chromosomeMap.end()) {
        chromosomes.push_back(static_cast<chromosome>(iter->second));
        iter++;
    }
    hiCFile->close();
    return chromosomes;
}
