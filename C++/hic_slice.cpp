#include "hic_slice.h"
#include "straw.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <zlib.h>

void writeCompressedBuffer(gzFile& file, const char* buffer, size_t size) {
    gzwrite(file, buffer, size);
}

void writeHeader(gzFile& file, const HicSliceHeader& header) {
    // Write magic string
    writeCompressedBuffer(file, HICSLICE_MAGIC.c_str(), HICSLICE_MAGIC.length());
    
    // Write resolution
    writeCompressedBuffer(file, (char*)&header.resolution, sizeof(int32_t));
    
    // Write number of chromosomes
    writeCompressedBuffer(file, (char*)&header.numChromosomes, sizeof(int32_t));
    
    // Write chromosome mapping
    for (const auto& chr : header.chromosomeKeys) {
        // Write chromosome name length
        int32_t nameLength = chr.first.length();
        writeCompressedBuffer(file, (char*)&nameLength, sizeof(int32_t));
        
        // Write chromosome name
        writeCompressedBuffer(file, chr.first.c_str(), nameLength);
        
        // Write chromosome key
        writeCompressedBuffer(file, (char*)&chr.second, sizeof(int16_t));
    }
}

void writeContactRecord(gzFile& file, const CompressedContactRecord& record) {
    writeCompressedBuffer(file, (char*)&record, sizeof(CompressedContactRecord));
}

void dumpGenomeWideDataAtResolution(const std::string& matrixType,
                                  const std::string& norm,
                                  const std::string& filePath,
                                  const std::string& unit,
                                  int32_t resolution,
                                  const std::string& outputPath) {
    // Open HiC file
    HiCFile* hicFile = new HiCFile(filePath);
    
    // Create header
    HicSliceHeader header;
    header.resolution = resolution;
    
    // Get chromosomes and create mapping
    std::vector<chromosome> chromosomes = hicFile->getChromosomes();
    int16_t chrKey = 0;
    for (const auto& chr : chromosomes) {
        if (chr.index > 0) {  // Skip chromosomes with index <= 0
            header.chromosomeKeys[chr.name] = chrKey++;
        }
    }
    header.numChromosomes = header.chromosomeKeys.size();
    
    // Open output file
    gzFile outFile = gzopen(outputPath.c_str(), "wb");
    if (!outFile) {
        std::cerr << "Error: Could not open output file " << outputPath << std::endl;
        return;
    }
    
    // Write header
    writeHeader(outFile, header);
    
    // Process each chromosome pair
    for (size_t i = 0; i < chromosomes.size(); i++) {
        if (chromosomes[i].index <= 0) continue;
        
        for (size_t j = i; j < chromosomes.size(); j++) {
            if (chromosomes[j].index <= 0) continue;
            
            // Get matrix data
            MatrixZoomData* mzd = hicFile->getMatrixZoomData(
                chromosomes[i].name, 
                chromosomes[j].name, 
                matrixType, 
                norm, 
                unit, 
                resolution
            );
            
            if (!mzd->foundFooter) continue;
            
            // Get records for entire chromosome pair
            std::vector<contactRecord> records = mzd->getRecords(
                0, chromosomes[i].length,
                0, chromosomes[j].length
            );
            
            // Write records
            for (const auto& record : records) {
                CompressedContactRecord compressedRecord;
                compressedRecord.chr1Key = header.chromosomeKeys[chromosomes[i].name];
                compressedRecord.binX = record.binX;
                compressedRecord.chr2Key = header.chromosomeKeys[chromosomes[j].name];
                compressedRecord.binY = record.binY;
                compressedRecord.value = record.counts;
                
                writeContactRecord(outFile, compressedRecord);
            }
            
            delete mzd;
        }
    }
    
    // Close files and cleanup
    gzclose(outFile);
    delete hicFile;
} 