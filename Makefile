# MEMBER 3 - build file.  Targets: all, test, clean
CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

all: genisys gen_reads

genisys: src/genisys_app_3.cpp src/vantage_tree_1.hpp src/seq_distance_2.hpp src/read_loader_2.hpp src/linear_scan_3.hpp
	$(CXX) $(CXXFLAGS) src/genisys_app_3.cpp -o genisys

gen_reads: tools/read_simulator_3.cpp src/read_loader_2.hpp
	$(CXX) $(CXXFLAGS) tools/read_simulator_3.cpp -o gen_reads

test_metrics: tests/distance_tests_2.cpp src/seq_distance_2.hpp src/read_loader_2.hpp
	$(CXX) $(CXXFLAGS) tests/distance_tests_2.cpp -o test_metrics

test_edge: tests/tree_validation_3.cpp src/vantage_tree_1.hpp src/seq_distance_2.hpp src/linear_scan_3.hpp
	$(CXX) $(CXXFLAGS) tests/tree_validation_3.cpp -o test_edge

test: test_metrics test_edge
	./test_metrics
	./test_edge

clean:
	rm -f genisys gen_reads test_metrics test_edge tmp_test_reads.fasta

.PHONY: all test clean
