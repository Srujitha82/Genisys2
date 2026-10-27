// ============================================================================
//  MEMBER 3  -  Interactive console application (demo for the evaluation)
//  Usage:  ./genisys [reads_file] [hamming|levenshtein]
// ============================================================================
#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include "vantage_tree_1.hpp"
#include "seq_distance_2.hpp"
#include "read_loader_2.hpp"
#include "linear_scan_3.hpp"

using namespace genisys;

struct App {
    LoadResult data;
    std::string metricName = "levenshtein";
    DistFn metric = levenshtein;
    std::unique_ptr<VPTree> tree;
};

static void setMetric(App& a, const std::string& name) {
    if (name == "hamming") { a.metricName = name; a.metric = hamming; }
    else { a.metricName = "levenshtein"; a.metric = levenshtein; }
}

static void buildTree(App& a) {
    if (a.data.reads.empty()) { a.tree.reset(); std::cout << "No reads loaded.\n"; return; }
    auto t0 = std::chrono::steady_clock::now();
    a.tree = std::make_unique<VPTree>(a.data.reads, a.metric);
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::cout << "VP-tree built: n=" << a.tree->size() << ", height=" << a.tree->height()
              << ", distance calls=" << a.tree->buildDistanceCalls() << ", time=" << ms
              << " ms, metric=" << a.metricName << "\n";
}

static bool readQuery(std::string& q) {
    std::cout << "Enter query read (A/C/G/T/N): ";
    std::string raw;
    if (!(std::cin >> raw)) return false;
    q = normalizeRead(raw);
    if (!isValidDNA(q)) { std::cout << "Invalid DNA sequence.\n"; return false; }
    return true;
}

static bool readInt(const char* prompt, int& v, int minV) {
    std::cout << prompt;
    if (!(std::cin >> v) || v < minV) {
        std::cin.clear();
        std::cin.ignore(1000, '\n');
        std::cout << "Invalid number.\n";
        return false;
    }
    return true;
}

static bool askRevComp() {
    std::cout << "Also search reverse complement? (y/n): ";
    char yn = 'n';
    std::cin >> yn;
    return yn == 'y' || yn == 'Y';
}

// Merge forward + reverse-complement results, keeping the best distance per read.
static std::vector<Match> mergeBest(const std::vector<Match>& x, const std::vector<Match>& y) {
    std::map<int, int> best;
    for (auto* v : {&x, &y})
        for (auto& m : *v) {
            auto it = best.find(m.index);
            if (it == best.end() || m.dist < it->second) best[m.index] = m.dist;
        }
    std::vector<Match> out;
    for (auto& kv : best) out.push_back({kv.first, kv.second});
    std::sort(out.begin(), out.end());
    return out;
}

static void show(const App& a, const std::vector<Match>& r, size_t limit = 20) {
    if (r.empty()) { std::cout << "  (no matches)\n"; return; }
    for (size_t i = 0; i < r.size() && i < limit; ++i)
        std::cout << "  #" << r[i].index << " " << a.data.names[r[i].index] << "  dist=" << r[i].dist
                  << "  " << a.tree->read(r[i].index) << "\n";
    if (r.size() > limit) std::cout << "  ... and " << (r.size() - limit) << " more\n";
}

static void doRange(App& a) {
    std::string q; int tau;
    if (!readQuery(q) || !readInt("Max distance (tau): ", tau, 0)) return;
    bool rc = askRevComp();
    auto r = a.tree->rangeSearch(q, tau);
    long long calls = a.tree->lastDistanceCalls();
    if (rc) {
        auto r2 = a.tree->rangeSearch(reverseComplement(q), tau);
        calls += a.tree->lastDistanceCalls();
        r = mergeBest(r, r2);
    }
    std::cout << r.size() << " match(es), distance calls: " << calls << " (dataset size " << a.tree->size() << ")\n";
    show(a, r);
}

static void doKnn(App& a) {
    std::string q; int k;
    if (!readQuery(q) || !readInt("k: ", k, 1)) return;
    bool rc = askRevComp();
    auto r = a.tree->kNearest(q, k);
    long long calls = a.tree->lastDistanceCalls();
    if (rc) {
        auto r2 = a.tree->kNearest(reverseComplement(q), k);
        calls += a.tree->lastDistanceCalls();
        r = mergeBest(r, r2);
        if (r.size() > static_cast<size_t>(k)) r.resize(k);
    }
    std::cout << "distance calls: " << calls << " (dataset size " << a.tree->size() << ")\n";
    show(a, r);
}

static void doCompare(App& a) {
    std::string q; int tau;
    if (!readQuery(q) || !readInt("Max distance (tau): ", tau, 0)) return;
    auto t0 = std::chrono::steady_clock::now();
    auto fast = a.tree->rangeSearch(q, tau);
    auto t1 = std::chrono::steady_clock::now();
    auto slow = bruteRange(a.data.reads, q, tau, a.metric);
    auto t2 = std::chrono::steady_clock::now();
    bool same = fast.size() == slow.size() &&
                std::equal(fast.begin(), fast.end(), slow.begin(),
                           [](const Match& x, const Match& y) { return x.index == y.index && x.dist == y.dist; });
    std::cout << "VP-tree : " << fast.size() << " matches, " << a.tree->lastDistanceCalls() << " distance calls, "
              << std::chrono::duration<double, std::micro>(t1 - t0).count() << " us\n";
    std::cout << "Brute   : " << slow.size() << " matches, " << a.data.reads.size() << " distance calls, "
              << std::chrono::duration<double, std::micro>(t2 - t1).count() << " us\n";
    std::cout << "Results identical: " << (same ? "YES" : "NO") << "\n";
}

static void menu() {
    std::cout << "\n===== Genisys2 : DNA read matching with VP-Tree =====\n"
              << "1. Load reads file (FASTA / one-per-line)\n"
              << "2. Choose metric (hamming / levenshtein) and rebuild\n"
              << "3. Range search  (all reads within distance tau)\n"
              << "4. k-nearest reads\n"
              << "5. Compare VP-tree vs brute force\n"
              << "6. Tree statistics\n"
              << "0. Exit\n> ";
}

int main(int argc, char** argv) {
    App app;
    if (argc >= 3) setMetric(app, argv[2]);
    if (argc >= 2) {
        app.data = loadReads(argv[1]);
        if (!app.data.ok) std::cout << "Cannot open " << argv[1] << "\n";
        else {
            std::cout << "Loaded " << app.data.reads.size() << " reads (" << app.data.skipped << " skipped)\n";
            buildTree(app);
        }
    }
    int choice = -1;
    while (true) {
        menu();
        if (!(std::cin >> choice)) {
            if (std::cin.eof()) break;
            std::cin.clear(); std::cin.ignore(1000, '\n');
            std::cout << "Please enter a number.\n";
            continue;
        }
        if (choice == 0) break;
        if (choice == 1) {
            std::string path; std::cout << "File path: "; std::cin >> path;
            LoadResult nd = loadReads(path);
            if (!nd.ok) { std::cout << "Cannot open file.\n"; continue; }
            app.data = std::move(nd);
            std::cout << "Loaded " << app.data.reads.size() << " reads (" << app.data.skipped << " skipped)\n";
            buildTree(app);
        } else if (choice == 2) {
            std::string m; std::cout << "Metric (hamming/levenshtein): "; std::cin >> m;
            setMetric(app, m); buildTree(app);
        } else if (choice >= 3 && choice <= 6) {
            if (!app.tree) { std::cout << "Load reads first (option 1).\n"; continue; }
            if (choice == 3) doRange(app);
            else if (choice == 4) doKnn(app);
            else if (choice == 5) doCompare(app);
            else std::cout << "n=" << app.tree->size() << " height=" << app.tree->height()
                           << " build distance calls=" << app.tree->buildDistanceCalls()
                           << " metric=" << app.metricName << "\n";
        } else std::cout << "Unknown option.\n";
    }
    std::cout << "Bye!\n";
    return 0;
}
