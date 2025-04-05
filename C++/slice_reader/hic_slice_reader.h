#ifndef HIC_SLICE_READER_H
#define HIC_SLICE_READER_H

#include <string>
#include <map>
#include <vector>

struct SliceContactRecord {
    std::string chr1;
    int32_t binX;
    std::string chr2;
    int32_t binY;
    float value;
};

class HicSliceReader {
public:
    explicit HicSliceReader(const std::string& filePath);
    ~HicSliceReader();

    // Get basic file information
    int32_t getResolution() const { return resolution; }
    std::vector<std::string> getChromosomes() const;
    
    // Read all records
    std::vector<SliceContactRecord> readAllRecords();
    
    // Read records for specific chromosome pair
    std::vector<SliceContactRecord> readRecordsForChromosomePair(const std::string& chr1, const std::string& chr2);

private:
    void readHeader();
    std::string getChromosomeFromKey(int16_t key) const;

    gzFile file;
    int32_t resolution;
    int32_t numChromosomes;
    std::map<int16_t, std::string> chromosomeKeyToName;
    bool headerRead;
};

#endif 