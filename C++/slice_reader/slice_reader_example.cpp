#include "hic_slice_reader.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <hicslice_file>" << std::endl;
        return 1;
    }

    try {
        HicSliceReader reader(argv[1]);
        
        // Print basic information
        std::cout << "Resolution: " << reader.getResolution() << std::endl;
        
        std::cout << "Chromosomes:" << std::endl;
        for (const auto& chr : reader.getChromosomes()) {
            std::cout << chr << std::endl;
        }

        // Example 1: Read all records
        std::cout << "\nReading all records:" << std::endl;
        auto records = reader.readAllRecords();
        std::cout << "Total records: " << records.size() << std::endl;
        
        // Print first few records
        for (size_t i = 0; i < std::min(size_t(5), records.size()); i++) {
            const auto& rec = records[i];
            std::cout << rec.chr1 << "\t" << rec.binX << "\t"
                     << rec.chr2 << "\t" << rec.binY << "\t"
                     << rec.value << std::endl;
        }

        // Example 2: Read specific chromosome pair
        std::cout << "\nReading chr1-chr2 records:" << std::endl;
        auto chr1chr2_records = reader.readRecordsForChromosomePair("chr1", "chr2");
        std::cout << "Chr1-Chr2 records: " << chr1chr2_records.size() << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
} 