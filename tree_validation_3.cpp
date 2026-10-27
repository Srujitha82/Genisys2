// MEMBER 3 - VP-tree vs brute-force correctness + edge cases.
#include <iostream>
#include <random>
#include "../src/vantage_tree_1.hpp"
#include "../src/seq_distance_2.hpp"
#include "../src/linear_scan_3.hpp"
using namespace genisys;
static int fails = 0;
#define CHECK(c) do { if (!(c)) { ++fails; std::cout << "FAIL line " << __LINE__ << ": " #c "\n"; } } while (0)

static bool sameDist(const std::vector<Match>& a, const std::vector<Match>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) if (a[i].dist != b[i].dist) return false;  // ties may pick other index
    return true;
}
static bool sameExact(const std::vector<Match>& a, const std::vector<Match>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) if (a[i].index != b[i].index || a[i].dist != b[i].dist) return false;
    return true;
}

int main() {
    // --- edge cases ---
    { VPTree t({}, levenshtein);
      CHECK(t.size() == 0 && t.height() == 0);
      CHECK(t.rangeSearch("ACGT", 5).empty());
      CHECK(t.kNearest("ACGT", 3).empty()); }
    { VPTree t({"ACGT"}, levenshtein);
      CHECK(t.height() == 1);
      auto r = t.rangeSearch("ACGT", 0);   CHECK(r.size() == 1 && r[0].dist == 0);
      CHECK(t.rangeSearch("TTTT", 1).empty());
      CHECK(t.kNearest("TTTT", 5).size() == 1); }
    { std::vector<std::string> same4(50, "ACGTACGT");   // all identical reads
      VPTree t(same4, levenshtein);
      CHECK(t.rangeSearch("ACGTACGT", 0).size() == 50);
      CHECK(t.kNearest("ACGTACGA", 10).size() == 10);
      CHECK(t.height() <= 8); }
    { std::vector<std::string> d = {"AAAA", "AAAT", "AATT", "ATTT", "TTTT"};
      VPTree t(d, levenshtein);
      CHECK(t.rangeSearch("AAAA", -1).empty());      // negative tau
      CHECK(t.kNearest("AAAA", 0).empty());          // k = 0
      CHECK(t.kNearest("AAAA", -3).empty());         // k < 0
      CHECK(t.kNearest("AAAA", 99).size() == 5);     // k > n
      CHECK(t.rangeSearch("AAAA", 0).size() == 1);
      CHECK(t.rangeSearch("AAAA", 100).size() == 5);
      CHECK(t.kNearest("AAAA", 2)[1].dist == 1); }
    { std::vector<std::string> d = {"A", "ACGTACGTAC", "AC", "", "GGGGGGGGGGGG"};  // mixed lengths
      VPTree t(d, levenshtein);
      CHECK(sameExact(t.rangeSearch("ACG", 3), bruteRange(d, "ACG", 3, levenshtein))); }

    // --- randomized comparison with brute force ---
    std::mt19937 rng(2024);
    auto rs = [&](int lo, int hi) { std::string s; int L = lo + rng() % (hi - lo + 1); while (L--) s += "ACGT"[rng() % 4]; return s; };
    for (int metricId = 0; metricId < 2; ++metricId) {
        DistFn m = metricId ? DistFn(hamming) : DistFn(levenshtein);
        for (int n : {2, 3, 10, 100, 500}) {
            std::vector<std::string> data;
            for (int i = 0; i < n; ++i) data.push_back(rs(8, 16));
            VPTree t(data, m);
            for (int qi = 0; qi < 60; ++qi) {
                std::string q = rs(8, 16);
                int tau = rng() % 8;
                int k = 1 + rng() % 12;
                CHECK(sameExact(t.rangeSearch(q, tau), bruteRange(data, q, tau, m)));
                CHECK(sameDist(t.kNearest(q, k), bruteKnn(data, q, k, m)));
            }
        }
    }

    // --- pruning actually helps ---
    { std::vector<std::string> data;
      for (int i = 0; i < 2000; ++i) data.push_back(rs(30, 30));
      VPTree t(data, hamming);
      t.rangeSearch(data[0], 2);
      std::cout << "info: n=2000, tau=2 -> " << t.lastDistanceCalls() << " distance calls (brute = 2000)\n";
      std::cout << "info: tree height = " << t.height() << " (log2(2000) ~ 11)\n"; }

    std::cout << (fails ? "SOME TESTS FAILED\n" : "m3 tests: ALL PASSED\n");
    return fails ? 1 : 0;
}
