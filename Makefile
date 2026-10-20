CXX := g++-16
AR ?= ar
EIGEN_INCLUDE ?= /Users/xtao/local/include

CPPFLAGS ?= -Iinclude -I$(EIGEN_INCLUDE)
CXXFLAGS ?= -std=c++20 -O2 -DEIGEN_NO_DEBUG -Wall -Wextra -Wno-class-memaccess

BUILD_DIR ?= build/sayram2d
LIB := $(BUILD_DIR)/libsayram2d.a

SRCS := src/mesh.cpp src/boundary_conditions.cpp src/solver.cpp
OBJS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

TEST_BIN := $(BUILD_DIR)/test_sayram2d
REF_BIN := $(BUILD_DIR)/new_reference
REF_DIR := build/new_reference
BENCH_BIN := $(BUILD_DIR)/benchmark_new
TIMED_NEW := $(BUILD_DIR)/timed_new/Solver.cpp
AY_BIN := $(BUILD_DIR)/albert_young_reference
AY_DATA := data/albert_young/ALBERT_YOUNG_80x80.bin
AY_NEW := tests/albert_young/new_reference.csv

.PHONY: all lib test reference albert-young compare benchmark-new clean

all: lib

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(TEST_BIN): tests/test_sayram2d.cc $(LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(REF_BIN): tests/new_reference.cc $(LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(AY_BIN): tests/albert_young_reference.cc tests/albert_young_case.h $(LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

$(TIMED_NEW): src/solver.cpp tests/make_timed_solver_new.py
	@mkdir -p $(BUILD_DIR)/timed_new
	python3 tests/make_timed_solver_new.py $(TIMED_NEW)

$(BENCH_BIN): tests/benchmark_new.cc $(TIMED_NEW) src/mesh.cpp src/boundary_conditions.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) tests/benchmark_new.cc src/mesh.cpp \
		src/boundary_conditions.cpp $(TIMED_NEW) -o $@

reference: $(REF_BIN)
	$(REF_BIN) $(REF_DIR)

albert-young: $(AY_BIN)
	$(AY_BIN) $(AY_DATA) $(AY_NEW)
	python3 tests/compare_albert_young.py \
		tests/albert_young/legacy_reference.csv $(AY_NEW)

test: $(TEST_BIN) $(REF_BIN) $(AY_BIN)
	$(TEST_BIN)
	$(REF_BIN) $(REF_DIR)
	python3 tests/compare_legacy.py $(REF_DIR) tests/legacy/reference
	$(AY_BIN) $(AY_DATA) $(AY_NEW)
	python3 tests/compare_albert_young.py \
		tests/albert_young/legacy_reference.csv $(AY_NEW)

compare: test

benchmark-new: $(BENCH_BIN)
	rm -f tests/results/benchmark_new.csv
	$(BENCH_BIN) 50 50 3 tests/results/benchmark_new.csv
	$(BENCH_BIN) 100 100 3 tests/results/benchmark_new.csv
	$(BENCH_BIN) 200 200 3 tests/results/benchmark_new.csv
	$(BENCH_BIN) 500 500 3 tests/results/benchmark_new.csv

clean:
	rm -rf $(BUILD_DIR) $(REF_DIR) $(AY_NEW)
