# Genisys2 test cases (all results verified by running the code)

## A. Genisys2 Lab (`visual/genisys_lab_1.html`)

Open the page and press **F5** so it loads its default 11 reads (Levenshtein metric):

`TTCATT GATAGA ATATTT TACGTT GTTAGC ATTTTT CTCATT GTTATC ATATCT TTTATT GGTATC`

Type the query in the **CHECK A DNA READ** bar, choose the mode, press **Check it**.

| # | Case | Query | Mode | Expected result | Distance calls |
|---|------|-------|------|-----------------|----------------|
| 1 | Exact match | `GATAGA` | Range τ = 0 | GATAGA (d=0) | 3 of 11 |
| 2 | One substitution | `CTCATG` | Range τ = 1 | CTCATT (d=1) | 9 of 11 |
| 3 | Two substitutions | `GTTAGG` | Range τ = 2 | GTTAGC (d=1), GATAGA (d=2), GTTATC (d=2) | 11 of 11 |
| 4 | One insertion | `GATTAGA` | Range τ = 1 | GATAGA (d=1) | 3 of 11 |
| 5 | One deletion | `GATGA` | Range τ = 1 | GATAGA (d=1) | 3 of 11 |
| 6 | Unknown base N | `GANAGA` | Range τ = 1 | GATAGA (d=1) | 3 of 11 |
| 7 | No match | `AAAAAAAAAA` | Range τ = 1 | no results | 3 of 11 |
| 8 | Longest query (20 letters) | `ACGTACGTACGTACGTACGT` | Range τ = 2 | no results | 3 of 11 |
| 9 | Nearest read | `TTCATC` | k-nearest k = 1 | TTCATT (d=1) | 6 of 11 |
| 10 | 3 nearest | `ACGTAC` | k-nearest k = 3 | TACGTT (d=3), GGTATC (d=3), ATATTT (d=4) | 11 of 11 |
| 11 | 5 nearest | `GATAGA` | k-nearest k = 5 | GATAGA (0), GTTAGC (2), GTTATC (3), GGTATC (3), ATATTT (4) | 11 of 11 |

Every row ends with **"Brute-force check: identical distances"**.

Metric switch: change *Metric* to **Hamming** and run case 4 (`GATTAGA`, τ = 1): no results, because Hamming sees an
inserted letter as many differences. Levenshtein finds GATAGA. This shows why the metric matters.

Database editing:

| Action | Input | Expected |
|--------|-------|----------|
| Add one read | `TTGACCA` then **Add** | "Added 1 read. Tree rebuilt: 12 reads ..." |
| Then check it | query `TTGACCT`, k-nearest k = 1 | TTGACCA (d=1) |
| Add a duplicate | `GATAGA` | "Not added: GATAGA (already there)." |
| Add invalid text | `HELLO` | "Not added: HELLO. Use only the letters A, C, G, T (N = unknown)." |
| Add too long | 25 letters | "... (longer than 20)" |
| Delete | click **×** on any read | "Deleted ... A VP-tree is static, so the tree was rebuilt ..." |
| Delete all, then Check | any query | "The database is empty. Add some reads first ..." |
| Empty query box, then Check | (nothing) | "Type a DNA read first ..." |

## B. C++ console program (`genisys.exe`) on `data\sample_reads_3.fasta` (200 reads)

Start: `genisys.exe data\sample_reads_3.fasta levenshtein`
Expected first lines: `Loaded 200 reads (0 skipped)` and `VP-tree built: n=200, height=8, distance calls=1153`.

Query base read: `CCATCCGACTTAGTCTCCGGATAACGGAATCAAGCTGACG` (read_1).
One-letter change: `CCATCAGACTTAGTCTCCGGATAACGGAATCAAGCTGACG`.
Reverse complement of read_1: `CGTCAGCTTGATTCCGTTATCCGGAGACTAAGTCGGATGG`.

| # | Menu | Input after the menu number | Expected output |
|---|------|-----------------------------|-----------------|
| 1 | 3 | read_1, tau `0`, reverse complement `n` | `1 match(es), distance calls: 8` and `#0 read_1 dist=0` |
| 2 | 3 | one-letter change, tau `2`, `n` | `1 match(es), distance calls: 11` and `#0 read_1 dist=1` |
| 3 | 4 | one-letter change, k `3`, `n` | `#0 read_1 dist=1`, `#116 read_117 dist=7`, `#174 read_175 dist=17` |
| 4 | 3 | reverse complement, tau `0`, `y` | `1 match(es)` and `#0 read_1 dist=0` |
| 5 | 3 | reverse complement, tau `0`, `n` | `0 match(es)` and `(no matches)` |
| 6 | 5 | one-letter change, tau `2` | `VP-tree : 1 matches, 111 distance calls` / `Brute : 1 matches, 200 distance calls` / `Results identical: YES` |
| 7 | 6 | none | `n=200 height=8 build distance calls=1153 metric=levenshtein` |

Cases 4 and 5 show why the reverse-complement option exists: the read is found only when it is on.

Bad input (the program must not crash):

| Input | Expected message |
|-------|------------------|
| query `ACGX` | `Invalid DNA sequence.` |
| tau `-5` or `abc` | `Invalid number.` |
| k `abc` | `Invalid number.` |
| menu choice `hello` | `Please enter a number.` |
| menu choice `99` | `Unknown option.` |
| menu 1 with a missing file | `Cannot open file.` |
| start with a missing file, then menu 3 | `Cannot open ...` then `Load reads first (option 1).` |

## C. Automated tests

| Command | Expected last line |
|---------|--------------------|
| `test_distance.exe` | `m2 tests: ALL PASSED` |
| `test_tree.exe` | `m3 tests: ALL PASSED` |
| `run_demo_3.bat` | all five stages, then `ALL DONE - everything passed.` |
