
#pragma once
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>
#include <vector>

namespace genisys {

// Upper-case and strip whitespace.
inline std::string normalizeRead(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s)
        if (!std::isspace(c)) out.push_back(static_cast<char>(std::toupper(c)));
    return out;
}

// Valid read: non-empty, only A C G T N.
inline bool isValidDNA(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s)
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') return false;
    return true;
}

// Reverse complement (the opposite strand). N stays N.
inline std::string reverseComplement(const std::string& s) {
    std::string out(s.rbegin(), s.rend());
    for (char& c : out) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            default: break;
        }
    }
    return out;
}

// Classic Hamming distance - strings MUST have equal length.   O(L)
inline int hammingStrict(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) throw std::invalid_argument("hammingStrict: length mismatch");
    int d = 0;
    for (size_t i = 0; i < a.size(); ++i) d += (a[i] != b[i]);
    return d;
}

// Hamming for reads of different length: mismatches on the common prefix
// plus |len(a)-len(b)|. Equivalent to Hamming on strings padded with a symbol
// outside the alphabet, so it is still a metric.   O(max(L1,L2))
inline int hamming(const std::string& a, const std::string& b) {
    const size_t m = std::min(a.size(), b.size());
    int d = static_cast<int>(std::max(a.size(), b.size()) - m);
    for (size_t i = 0; i < m; ++i) d += (a[i] != b[i]);
    return d;
}

// Levenshtein edit distance (insert / delete / substitute, cost 1 each).
// Two-row dynamic programming: O(L1*L2) time, O(min(L1,L2)) space.
inline int levenshtein(const std::string& a, const std::string& b) {
    const std::string& s = (a.size() >= b.size()) ? a : b;  // longer
    const std::string& t = (a.size() >= b.size()) ? b : a;  // shorter
    std::vector<int> prev(t.size() + 1), cur(t.size() + 1);
    for (size_t j = 0; j <= t.size(); ++j) prev[j] = static_cast<int>(j);
    for (size_t i = 1; i <= s.size(); ++i) {
        cur[0] = static_cast<int>(i);
        for (size_t j = 1; j <= t.size(); ++j) {
            const int sub = prev[j - 1] + (s[i - 1] != t[j - 1]);
            cur[j] = std::min({sub, prev[j] + 1, cur[j - 1] + 1});
        }
        std::swap(prev, cur);
    }
    return prev[t.size()];
}

}  // namespace genisys
