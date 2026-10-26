# Design, justification and complexity (Member 1)

## 1. Why a VP-tree? (Data-structure justification)
A DNA read is **not a point in a coordinate space** – it is a string. The natural similarity
(edit distance / Hamming distance) is a **metric**, so we need an index that only uses `d(a,b)`.

| Alternative | Why it is worse for this problem |
|---|---|
| Linear scan / array | O(n) distance computations for every query (each is O(L²)) |
| KD-tree / R-tree | Needs numeric coordinates and axis-aligned splits; strings have none. Embedding to vectors loses the exact edit distance |
| Hash table | Only exact match; cannot answer "within 3 mismatches" |
| Trie | Exact / prefix match; approximate matching blows up (branching on every mismatch) |
| BK-tree | Works for integer metrics, but is **unbalanced** (depth depends on the data) and each node has many children; VP-tree guarantees a balanced binary tree by splitting at the median |
| k-mer inverted index | Good heuristic but not exact for arbitrary edit distance; needs tuned k |

**VP-tree** = balanced binary tree built only from distances. Each node stores a *vantage point* and the
*median radius* `mu`. Inside-ball reads go left, outside-ball reads go right.

## 2. Implementation strategy (array-based vs pointer-based)
Nodes are stored in one `std::vector<Node>`; children are **integer indices** (`-1` = none).
* ✔ contiguous memory → better cache behaviour, one allocation (`reserve(n)`), no per-node `new/delete`, no memory leaks / destructor recursion
* ✔ trivially copyable, easy to serialize
* ✘ tree is **static** – inserting/deleting needs a rebuild (acceptable: read databases are built once and queried many times)
* A pointer-based tree would allow cheap inserts but costs n allocations, pointer chasing and manual cleanup.

Other choices: reads are never copied around while building – only an `order_` array of indices is partitioned
(`std::nth_element`, O(n) per level). The vantage point is chosen **randomly with a fixed seed**
(reproducible; avoids adversarial worst cases). Recursion depth is O(log n) because of the median split, so no stack overflow.
kNN uses a **max-heap** of size k; the heap top is the current search radius `tau`.

## 3. Correctness (why pruning is valid)
Invariant: left ⇒ `d(vp,r) <= mu`, right ⇒ `d(vp,r) >= mu`. For query `q` and radius `tau`:
* A left read `r` has `d(q,r) >= d(q,vp) - d(vp,r) >= d(q,vp) - mu` (triangle inequality). If `d(q,vp) - tau > mu` no left read can be within `tau` → skip left.
* A right read has `d(q,r) >= d(vp,r) - d(q,vp) >= mu - d(q,vp)`. If `d(q,vp) + tau < mu` → skip right.
So we visit left iff `d(q,vp) - tau <= mu`, right iff `d(q,vp) + tau >= mu`. Ties at `mu` can legally sit on both sides, which the `<=` / `>=` conditions handle.

## 4. Time complexity
**Build.** At each node we compute `m-1` distances for a subtree of size m and run `nth_element` in O(m).
Median split ⇒ `T(n) = 2T(n/2) + Θ(n)` ⇒ `T(n) = Θ(n log n)` distance evaluations
(Master theorem, case 2). Total `O(n log n · C)`. Height = ⌈log₂(n+1)⌉ (verified: n=2000 → 11).

**Query.** One distance computation per visited node. Visiting a root-to-leaf path costs O(log n).
Pruning removes a subtree whenever the ball `B(q,tau)` does not intersect the corresponding region, so for small `tau`
(relative to the spread of the data) only O(log n) nodes are visited → **O(log n · C)** expected.
If `tau` is at least the diameter of the data, nothing can be pruned → **O(n · C)** worst case
(same as brute force; no metric index can beat this in the worst case).
Result sorting adds O(r log r) for r results; kNN heap operations add O(log k) per visited node.

**Space.** `n` nodes × 16 bytes + copy of reads = O(n) (plus O(log n) recursion stack).

## 5. Likely viva questions
* *Why median radius?* Guarantees both sides ≤ n/2 → height log n.
* *Why must the distance be a metric?* Pruning needs the triangle inequality; a non-metric distance would silently miss matches.
* *Is Hamming with different lengths a metric?* Our version = mismatches on the common prefix + |len difference| = Hamming on strings padded with an extra symbol → metric (checked on 2000 random triples in `distance_tests_2.cpp`).
* *Duplicate reads?* Distance 0, handled by the ≤/≥ tie rule; test with 50 identical reads.
* *Why not insert/delete?* Static tree by design (see §2); rebuild is O(n log n · C).
* *Curse of dimensionality?* On random strings of long length all distances concentrate near the mean, so pruning weakens; real reads with similarity structure prune well.
