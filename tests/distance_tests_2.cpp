// MEMBER 2 - unit tests for metrics and file I/O.
#include <cstdio>
#include <iostream>
#include <random>
#include "../src/seq_distance_2.hpp"
#include "../src/read_loader_2.hpp"
using namespace genisys;
static int fails = 0;
#define CHECK(c) do { if (!(c)) { ++fails; std::cout << "FAIL line " << __LINE__ << ": " #c "\n"; } } while (0)

int main() {
    CHECK(levenshtein("ACGT", "ACGT") == 0);
    CHECK(levenshtein("ACGT", "AGT") == 1);
    CHECK(levenshtein("", "ACG") == 3);
    CHECK(levenshtein("", "") == 0);
    CHECK(levenshtein("AAAA", "TTTT") == 4);
    CHECK(levenshtein("GATTACA", "GCATGCU") == levenshtein("GCATGCU", "GATTACA"));
    CHECK(hammingStrict("ACGT", "AGGT") == 1);
    bool threw = false;
    try { hammingStrict("AC", "ACG"); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);
    CHECK(hamming("ACGT", "ACGT") == 0);
    CHECK(hamming("ACGT", "ACG") == 1);
    CHECK(hamming("ACGT", "TCGA") == 2);
    CHECK(reverseComplement("AACG") == "CGTT");
    CHECK(reverseComplement(reverseComplement("ACGTN")) == "ACGTN");
    CHECK(isValidDNA("ACGTN") && !isValidDNA("") && !isValidDNA("ACGX") && !isValidDNA("acgt"));
    CHECK(normalizeRead(" ac gt\n") == "ACGT");

    // metric axioms on random strings (needed for VP-tree correctness)
    std::mt19937 rng(7);
    auto rs = [&]() { std::string s; int L = rng() % 12; while (L--) s += "ACGT"[rng() % 4]; return s; };
    for (int t = 0; t < 2000; ++t) {
        std::string a = rs(), b = rs(), c = rs();
        CHECK(levenshtein(a, c) <= levenshtein(a, b) + levenshtein(b, c));
        CHECK(hamming(a, c) <= hamming(a, b) + hamming(b, c));
        CHECK(levenshtein(a, b) == levenshtein(b, a));
        CHECK(hamming(a, b) == hamming(b, a));
    }

    // file I/O
    const char* path = "tmp_test_reads.fasta";
    CHECK(writeFasta(path, {"ACGT", "TTGA"}));
    auto r = loadReads(path);
    CHECK(r.ok && r.reads.size() == 2 && r.reads[1] == "TTGA");
    std::remove(path);
    { FILE* f = fopen(path, "w"); fputs(">s1\nACG\nTTA\n>bad\nAXZ\n>s3\nggg\n", f); fclose(f); }
    r = loadReads(path);
    CHECK(r.reads.size() == 2 && r.reads[0] == "ACGTTA" && r.reads[1] == "GGG" && r.skipped == 1);
    { FILE* f = fopen(path, "w"); fputs("acgt\n\nNNAC\nxx\n", f); fclose(f); }
    r = loadReads(path);
    CHECK(r.reads.size() == 2 && r.skipped == 1);
    std::remove(path);
    CHECK(!loadReads("does_not_exist.fa").ok);

    std::cout << (fails ? "SOME TESTS FAILED\n" : "m2 tests: ALL PASSED\n");
    return fails ? 1 : 0;
}
