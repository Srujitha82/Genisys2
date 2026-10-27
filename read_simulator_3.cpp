// MEMBER 3 - synthetic read generator.  Usage: ./gen_reads out.fasta N L [seed]
// Builds a random genome and samples N reads of length L with ~3% mutations.
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include "../src/read_loader_2.hpp"
int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: gen_reads out.fasta N L [seed]\n"; return 1; }
    int n = std::atoi(argv[2]), L = std::atoi(argv[3]);
    unsigned seed = argc > 4 ? std::atoi(argv[4]) : 42;
    if (n <= 0 || L <= 0) { std::cerr << "N and L must be positive\n"; return 1; }
    std::mt19937 rng(seed);
    const char* A = "ACGT";
    int G = std::max(L * 20, 2000);
    std::string genome(G, 'A');
    for (char& c : genome) c = A[rng() % 4];
    std::vector<std::string> reads;
    for (int i = 0; i < n; ++i) {
        std::string r = genome.substr(rng() % (G - L + 1), L);
        for (char& c : r) if (rng() % 100 < 3) c = A[rng() % 4];
        reads.push_back(r);
    }
    if (!genisys::writeFasta(argv[1], reads)) { std::cerr << "cannot write file\n"; return 1; }
    std::cout << "wrote " << n << " reads of length " << L << " to " << argv[1] << "\n";
}
