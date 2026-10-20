# Build the Sayram2D library, the standalone `sayram2d` app, or the small
# library example. Tests and benchmarks have separate Makefiles:
#
#   make -f Makefile.tests
#   make -f Makefile.benchmark
#
# A one-command case workflow is in Makefile.case.

CXX := g++-16
AR ?= ar
EIGEN_INCLUDE ?= /Users/xtao/local/include

HDF5_PREFIX ?= $(shell brew --prefix hdf5 2>/dev/null)
HDF5_INCLUDE ?= $(if $(HDF5_PREFIX),$(HDF5_PREFIX)/include,/usr/local/include)
HDF5_LIBDIR ?= $(if $(HDF5_PREFIX),$(HDF5_PREFIX)/lib,/usr/local/lib)
HDF5_LIBS ?= -lhdf5

CPPFLAGS ?= -Isrc/include -I$(EIGEN_INCLUDE)
APP_CPPFLAGS ?= -Isrc/app -I$(HDF5_INCLUDE)
HDF5_LDFLAGS ?= -L$(HDF5_LIBDIR)
CXXFLAGS ?= -std=c++20 -O2 -DEIGEN_NO_DEBUG -Wall -Wextra -Wno-class-memaccess

BUILD_DIR ?= build/sayram2d
LIB := $(BUILD_DIR)/libsayram2d.a

SRCS := src/mesh.cpp src/boundary_conditions.cpp src/solver.cpp
OBJS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

EXAMPLE_BIN := $(BUILD_DIR)/example
APP_BIN := sayram2d

.PHONY: all lib example app app-run clean

all: lib

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(EXAMPLE_BIN): src/main.cc $(LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

example: $(EXAMPLE_BIN)
	$(EXAMPLE_BIN)

$(APP_BIN): src/app/main.cc src/app/config.cc src/app/io.cc src/app/albert_young.cc $(LIB)
	$(CXX) $(APP_CPPFLAGS) $(CPPFLAGS) $(CXXFLAGS) src/app/main.cc \
		src/app/config.cc src/app/io.cc src/app/albert_young.cc $(LIB) \
		$(HDF5_LDFLAGS) $(HDF5_LIBS) -o $@

app: $(APP_BIN)

app-run: $(APP_BIN)
	./$(APP_BIN) app.ini

clean:
	rm -rf $(BUILD_DIR) $(APP_BIN)
