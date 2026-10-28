# Genisys2 — DNA Read Matching using a Vantage-Point (VP) Tree

Given a database of DNA reads (short sequences), find the reads closest to a query read:
* **Range search** – all reads within distance `tau` of the query
* **k-nearest search** – the `k` most similar reads
* Optional search on the **reverse-complement** strand (DNA is double-stranded)

Distances: **Levenshtein** (insert / delete / substitute) and **Hamming** (substitutions; length difference counted).
Both are true metrics, so the triangle inequality can be used to prune the VP-tree.

## Team and ownership
| Member | Files |
|---|---|
| Member 1 – VP-tree core | `src/vantage_tree_1.hpp`, `visual/dna_tree_visualizer_1.html`, `README.md`, `docs/design_analysis_1.md` |
| Member 2 – DNA metrics & I/O | `src/seq_distance_2.hpp`, `src/read_loader_2.hpp`, `tests/distance_tests_2.cpp` |
| Member 3 – App, tests, tooling | `src/genisys_app_3.cpp`, `src/linear_scan_3.hpp`, `tests/tree_validation_3.cpp`, `tools/read_simulator_3.cpp`, `data/sample_reads_3.fasta`, `Makefile` |

## Build & run
```bash
make all            # builds ./genisys and ./gen_reads
make test           # runs all unit / edge-case / brute-force comparison tests
./genisys data/sample_reads_3.fasta levenshtein
./gen_reads data/big.fasta 20000 100 3      # make a larger synthetic dataset
```
Requires `g++` with C++17 and `make`.



## Complexity (n reads, read length L, C = cost of one distance call)
| Operation | Time | Space |
|---|---|---|
| Build | O(n log n · C) | O(n) |
| Range / kNN query (typical, small tau) | O(log n · C) | O(log n) stack |
| Range / kNN query (worst case, tau ≥ data diameter) | O(n · C) | O(log n) |
| Brute force baseline | O(n · C) | O(1) |

C = O(L²) for Levenshtein, O(L) for Hamming. Details and viva answers: `docs/design_analysis_1.md`.
