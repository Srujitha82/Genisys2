// ============================================================================
//  MEMBER 2  -  Reading / writing reads (FASTA or one-read-per-line text)
// ============================================================================
#pragma once
#include <fstream>
#include <string>
#include <vector>
#include "seq_distance_2.hpp"

namespace genisys {

struct LoadResult {
    std::vector<std::string> reads;
    std::vector<std::string> names;
    int skipped = 0;   // invalid / empty reads that were ignored
    bool ok = false;   // false if the file could not be opened
};

// Supports multi-line FASTA ('>' headers) and plain files (one read per line).
inline LoadResult loadReads(const std::string& path) {
    LoadResult res;
    std::ifstream in(path);
    if (!in) return res;
    res.ok = true;

    std::string line, cur, name;
    bool inRecord = false;
    int plainCount = 0;

    auto flush = [&]() {
        if (!inRecord) return;
        std::string seq = normalizeRead(cur);
        if (isValidDNA(seq)) { res.reads.push_back(seq); res.names.push_back(name); }
        else ++res.skipped;
        cur.clear();
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == ';') continue;
        if (line[0] == '>') {
            flush();
            inRecord = true;
            name = line.substr(1);
        } else if (inRecord) {
            cur += line;
        } else {  // plain format
            std::string seq = normalizeRead(line);
            if (isValidDNA(seq)) {
                res.reads.push_back(seq);
                res.names.push_back("read_" + std::to_string(++plainCount));
            } else ++res.skipped;
        }
    }
    flush();
    return res;
}

inline bool writeFasta(const std::string& path, const std::vector<std::string>& reads) {
    std::ofstream out(path);
    if (!out) return false;
    for (size_t i = 0; i < reads.size(); ++i)
        out << ">read_" << (i + 1) << "\n" << reads[i] << "\n";
    return true;
}

}  // namespace genisys
