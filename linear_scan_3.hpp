// ============================================================================
//  MEMBER 3  -  Brute-force baseline (used for correctness checking)
// ============================================================================
#pragma once
#include <algorithm>
#include <string>
#include <vector>
#include "vantage_tree_1.hpp"

namespace genisys {

inline std::vector<Match> bruteRange(const std::vector<std::string>& data, const std::string& q,
                                     int tau, const DistFn& d) {
    std::vector<Match> out;
    for (size_t i = 0; i < data.size(); ++i) {
        int x = d(q, data[i]);
        if (x <= tau) out.push_back({static_cast<int>(i), x});
    }
    std::sort(out.begin(), out.end());
    return out;
}

inline std::vector<Match> bruteKnn(const std::vector<std::string>& data, const std::string& q,
                                   int k, const DistFn& d) {
    std::vector<Match> all;
    for (size_t i = 0; i < data.size(); ++i) all.push_back({static_cast<int>(i), d(q, data[i])});
    std::sort(all.begin(), all.end());
    if (k < 0) k = 0;
    if (all.size() > static_cast<size_t>(k)) all.resize(k);
    return all;
}

}  // namespace genisys
