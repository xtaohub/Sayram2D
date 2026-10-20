ROOT := ../..
CXX ?= g++-16
BUILD ?= $(ROOT)/build
TIMED_DIR ?= $(BUILD)/legacy_timed
INCLUDES ?= -I/opt/homebrew/include -I/Users/xtao/local/include
CXXFLAGS ?= -std=c++17 -O2 $(INCLUDES) -I$(ROOT)/source -I$(ROOT)/source/Cases

.PHONY: all reference benchmark run-reference run-benchmark clean

all: reference benchmark

reference: $(BUILD)/legacy_reference

benchmark: $(BUILD)/legacy_benchmark

$(BUILD)/legacy_reference: $(ROOT)/tests/legacy/legacy_reference.cc \
                           $(ROOT)/tests/legacy/legacy_case.h \
                           $(ROOT)/source/Mesh.cc $(ROOT)/source/Solver.cc
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(ROOT)/tests/legacy/legacy_reference.cc \
		$(ROOT)/source/Mesh.cc $(ROOT)/source/Solver.cc \
		-o $@

$(TIMED_DIR)/Solver.cc: $(ROOT)/source/Solver.cc \
                        $(ROOT)/tests/legacy/make_timed_solver.py
	@mkdir -p $(TIMED_DIR)
	cd $(ROOT) && python3 tests/legacy/make_timed_solver.py $(abspath $@)

$(BUILD)/legacy_benchmark: $(ROOT)/tests/legacy/benchmark_original.cc \
                           $(ROOT)/tests/legacy/legacy_case.h \
                           $(TIMED_DIR)/Solver.cc $(ROOT)/source/Mesh.cc
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) $(ROOT)/tests/legacy/benchmark_original.cc \
		$(ROOT)/source/Mesh.cc $(TIMED_DIR)/Solver.cc \
		-o $@

run-reference: reference
	cd $(ROOT) && ./build/legacy_reference tests/legacy/reference

run-benchmark: benchmark
	rm -f $(ROOT)/tests/legacy/results/benchmark_original.csv
	cd $(ROOT) && ./build/legacy_benchmark 50 50 3 tests/legacy/results/benchmark_original.csv
	cd $(ROOT) && ./build/legacy_benchmark 100 100 3 tests/legacy/results/benchmark_original.csv
	cd $(ROOT) && ./build/legacy_benchmark 200 200 3 tests/legacy/results/benchmark_original.csv
	cd $(ROOT) && ./build/legacy_benchmark 500 500 3 tests/legacy/results/benchmark_original.csv

clean:
	rm -rf $(BUILD)/legacy_reference $(BUILD)/legacy_benchmark $(TIMED_DIR)
