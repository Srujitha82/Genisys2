
#pragma once
#include <algorithm>
#include <climits>
#include <cstddef>
#include <functional>
#include <queue>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace genisys {

using DistFn = std::function<int(const std::string&, const std::string&)>;

struct Match {
    int index;  // position of the read in the original dataset
    int dist;   // distance from the query
    bool operator<(const Match& o) const {
        return dist != o.dist ? dist < o.dist : index < o.index;
    }
};

class VPTree {
public:
    VPTree(const std::vector<std::string>& reads, DistFn fn)
        : data_(reads), dist_(std::move(fn)) { build(); }

    size_t size() const { return data_.size(); }
    int height() const { return height_; }
    long long lastDistanceCalls() const { return calls_; }
    long long buildDistanceCalls() const { return buildCalls_; }
    const std::string& read(int i) const { return data_[i]; }

    // All reads with distance <= tau, sorted by (dist, index).
    std::vector<Match> rangeSearch(const std::string& q, int tau) const {
        calls_ = 0;
        std::vector<Match> out;
        if (root_ >= 0 && tau >= 0) range(root_, q, tau, out);
        std::sort(out.begin(), out.end());
        return out;
    }

    // k closest reads, sorted by (dist, index).
    std::vector<Match> kNearest(const std::string& q, int k) const {
        calls_ = 0;
        std::vector<Match> out;
        if (root_ < 0 || k <= 0) return out;
        std::priority_queue<Match> heap;  // max-heap: worst candidate on top
        knn(root_, q, k, heap);
        while (!heap.empty()) { out.push_back(heap.top()); heap.pop(); }
        std::reverse(out.begin(), out.end());
        return out;
    }

private:
    // Array-based node: children are indices into nodes_, -1 means "none".
    struct Node { int point; int mu; int left; int right; };

    std::vector<std::string> data_;
    DistFn dist_;
    std::vector<Node> nodes_;
    std::vector<int> order_;
    int root_ = -1;
    int height_ = 0;
    mutable long long calls_ = 0;
    long long buildCalls_ = 0;

    int d(const std::string& a, const std::string& b) const {
        ++calls_;
        return dist_(a, b);
    }

    void build() {
        const int n = static_cast<int>(data_.size());
        order_.resize(n);
        for (int i = 0; i < n; ++i) order_[i] = i;
        nodes_.reserve(n);
        std::mt19937 rng(12345);  // fixed seed -> reproducible tree
        root_ = buildRange(0, n, 0, rng);
        buildCalls_ = calls_;
        calls_ = 0;
        order_.clear();
        order_.shrink_to_fit();
    }

    // Builds subtree over order_[lo, hi). Returns node index or -1.
    int buildRange(int lo, int hi, int depth, std::mt19937& rng) {
        if (lo >= hi) return -1;
        height_ = std::max(height_, depth + 1);

        std::uniform_int_distribution<int> pick(lo, hi - 1);
        std::swap(order_[lo], order_[pick(rng)]);  // random vantage point
        const int vp = order_[lo];

        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back({vp, 0, -1, -1});
        if (hi - lo == 1) return id;

        std::vector<std::pair<int, int>> tmp;  // (distance to vp, read index)
        tmp.reserve(hi - lo - 1);
        for (int i = lo + 1; i < hi; ++i)
            tmp.push_back({d(data_[vp], data_[order_[i]]), order_[i]});

        const size_t mid = (tmp.size() - 1) / 2;  // median position
        std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
        const int mu = tmp[mid].first;
        for (size_t i = 0; i < tmp.size(); ++i) order_[lo + 1 + i] = tmp[i].second;

        const int leftEnd = lo + 1 + static_cast<int>(mid) + 1;  // left = [lo+1, leftEnd)
        const int l = buildRange(lo + 1, leftEnd, depth + 1, rng);
        const int r = buildRange(leftEnd, hi, depth + 1, rng);
        nodes_[id].mu = mu;
        nodes_[id].left = l;
        nodes_[id].right = r;
        return id;
    }

    void range(int ni, const std::string& q, int tau, std::vector<Match>& out) const {
        const Node nd = nodes_[ni];
        const int dq = d(q, data_[nd.point]);
        if (dq <= tau) out.push_back({nd.point, dq});
        if (nd.left >= 0 && dq - tau <= nd.mu) range(nd.left, q, tau, out);
        if (nd.right >= 0 && dq + tau >= nd.mu) range(nd.right, q, tau, out);
    }

    void knn(int ni, const std::string& q, int k, std::priority_queue<Match>& heap) const {
        const Node nd = nodes_[ni];
        const int dq = d(q, data_[nd.point]);
        const Match cand{nd.point, dq};
        if (heap.size() < static_cast<size_t>(k)) heap.push(cand);
        else if (cand < heap.top()) { heap.pop(); heap.push(cand); }

        auto tau = [&]() -> long long {
            return heap.size() < static_cast<size_t>(k) ? (long long)INT_MAX : (long long)heap.top().dist;
        };
        auto goLeft  = [&]() { return nd.left  >= 0 && (long long)dq - tau() <= nd.mu; };
        auto goRight = [&]() { return nd.right >= 0 && (long long)dq + tau() >= nd.mu; };

        if (dq <= nd.mu) {                       // query is inside the ball: left first
            if (goLeft())  knn(nd.left, q, k, heap);
            if (goRight()) knn(nd.right, q, k, heap);
        } else {                                 // outside: right first
            if (goRight()) knn(nd.right, q, k, heap);
            if (goLeft())  knn(nd.left, q, k, heap);
        }
    }
};

}  // namespace genisys
