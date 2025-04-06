#include "hic_slice_reader.h"
#include <zlib.h>
#include <stdexcept>
#include <iostream>
#include <cstring>
#include <cmath>  // for isnan, isinf

HicSliceReader::HicSliceReader(const std::string& filePath) : headerRead(false) {
    file = gzopen(filePath.c_str(), "rb");
    if (!file) {
        throw std::runtime_error("Could not open file: " + filePath);
    }
    readHeader();
}

HicSliceReader::~HicSliceReader() {
    if (file) {
        gzclose(file);
    }
}

void HicSliceReader::readHeader() {
    if (headerRead) return;

    // Read and verify magic string
    char magic[8];
    if (gzread(file, magic, 8) != 8 || strncmp(magic, "HICSLICE", 8) != 0) {
        throw std::runtime_error("Invalid file format: missing magic string");
    }

    // Read resolution
    if (gzread(file, &resolution, sizeof(int32_t)) != sizeof(int32_t)) {
        throw std::runtime_error("Failed to read resolution");
    }

    // Read number of chromosomes
    if (gzread(file, &numChromosomes, sizeof(int32_t)) != sizeof(int32_t)) {
        throw std::runtime_error("Failed to read chromosome count");
    }

    // Read chromosome mapping
    for (int i = 0; i < numChromosomes; i++) {
        int32_t nameLength;
        if (gzread(file, &nameLength, sizeof(int32_t)) != sizeof(int32_t)) {
            throw std::runtime_error("Failed to read chromosome name length");
        }

        std::vector<char> nameBuffer(nameLength + 1, 0);
        if (gzread(file, nameBuffer.data(), nameLength) != nameLength) {
            throw std::runtime_error("Failed to read chromosome name");
        }
        std::string chromosomeName(nameBuffer.data());

        int16_t key;
        if (gzread(file, &key, sizeof(int16_t)) != sizeof(int16_t)) {
            throw std::runtime_error("Failed to read chromosome key");
        }

        chromosomeKeyToName[key] = chromosomeName;
    }

    headerRead = true;
}

std::vector<std::string> HicSliceReader::getChromosomes() const {
    std::vector<std::string> chromosomes;
    chromosomes.reserve(chromosomeKeyToName.size());
    for (const auto& pair : chromosomeKeyToName) {
        chromosomes.push_back(pair.second);
    }
    return chromosomes;
}

std::string HicSliceReader::getChromosomeFromKey(int16_t key) const {
    auto it = chromosomeKeyToName.find(key);
    if (it == chromosomeKeyToName.end()) {
        throw std::runtime_error("Invalid chromosome key: " + std::to_string(key));
    }
    return it->second;
}

std::vector<SliceContactRecord> HicSliceReader::readAllRecords() {
    std::vector<SliceContactRecord> records;
    
    // Structure to read from file
    struct {
        int16_t chr1Key;
        int32_t binX;
        int16_t chr2Key;
        int32_t binY;
        float value;
    } compressedRecord;

    while (gzread(file, &compressedRecord, sizeof(compressedRecord)) == sizeof(compressedRecord)) {
        SliceContactRecord record;
        record.chr1 = getChromosomeFromKey(compressedRecord.chr1Key);
        record.binX = compressedRecord.binX;
        record.chr2 = getChromosomeFromKey(compressedRecord.chr2Key);
        record.binY = compressedRecord.binY;
        record.value = compressedRecord.value;
        records.push_back(record);
    }

    return records;
}

std::vector<SliceContactRecord> HicSliceReader::readRecordsForChromosomePair(
    const std::string& chr1, const std::string& chr2) {
    
    std::vector<SliceContactRecord> records;
    
    // Find chromosome keys
    int16_t chr1Key = -1, chr2Key = -1;
    for (const auto& pair : chromosomeKeyToName) {
        if (pair.second == chr1) chr1Key = pair.first;
        if (pair.second == chr2) chr2Key = pair.first;
    }
    
    if (chr1Key == -1 || chr2Key == -1) {
        throw std::runtime_error("Chromosome not found in file");
    }

    // Structure to read from file
    struct {
        int16_t chr1Key;
        int32_t binX;
        int16_t chr2Key;
        int32_t binY;
        float value;
    } compressedRecord;

    while (gzread(file, &compressedRecord, sizeof(compressedRecord)) == sizeof(compressedRecord)) {
        // Only look for records in the correct orientation (chr1 <= chr2)
        if (compressedRecord.chr1Key == chr1Key && compressedRecord.chr2Key == chr2Key &&
            compressedRecord.value > 0 && !std::isnan(compressedRecord.value) && !std::isinf(compressedRecord.value)) {
            
            SliceContactRecord record;
            record.chr1 = getChromosomeFromKey(compressedRecord.chr1Key);
            record.binX = compressedRecord.binX;
            record.chr2 = getChromosomeFromKey(compressedRecord.chr2Key);
            record.binY = compressedRecord.binY;
            record.value = compressedRecord.value;
            records.push_back(record);
        }
    }

    // Reset file position to start of records for next read
    gzrewind(file);
    readHeader(); // Skip header again

    return records;
} 