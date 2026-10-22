# Simple build for the standalone `sayram2d` application.
#
#   make          build ./sayram2d
#   make app      same
#   make lib      build the static library (requires ar)
#   make example1  build and run the time-independent paper problem
#   make example2  build and run the time-dependent example
#   make clean     remove temporary Make build files
#   make distclean remove temporary files, the static library, and the app
#
# Edit the four paths below if your installation is in a different place.
# Tests are configured with CMake/CTest.

CXX ?= c++
AR ?= ar

EIGEN_INCLUDE ?= /Users/xtao/local/include
HDF5_INCLUDE ?= /opt/homebrew/opt/hdf5/include
HDF5_LIBDIR ?= /opt/homebrew/opt/hdf5/lib
HDF5_LIBS ?= -lhdf5

CPPFLAGS := -Isrc/include -I$(EIGEN_INCLUDE)
APP_CPPFLAGS := -Isrc/app -Isrc/examples/albert_young -I$(HDF5_INCLUDE)
CXXFLAGS ?= -std=c++20 -O2 -DEIGEN_NO_DEBUG -Wall -Wextra -Wno-class-memaccess
HDF5_LDFLAGS := -L$(HDF5_LIBDIR)

BUILD_DIR ?= build/sayram2d
BUILD_LIB := $(BUILD_DIR)/libsayram2d.a
LIB := lib/libsayram2d.a
APP_BIN := sayram2d
EXAMPLE1_BIN := $(BUILD_DIR)/example1
EXAMPLE2_BIN := $(BUILD_DIR)/example2

CORE_SRCS := src/mesh.cc src/boundary_conditions.cc src/solver.cc
APP_SRCS := src/examples/albert_young/main.cc \
            src/app/config.cc src/app/io.cc \
            src/examples/albert_young/albert_young.cc $(CORE_SRCS)
HEADERS := $(wildcard src/include/*.h) $(wildcard src/app/*.h) \
           $(wildcard src/examples/albert_young/*.h)
OBJS := $(patsubst src/%.cc,$(BUILD_DIR)/%.o,$(CORE_SRCS))

.PHONY: all app lib example1 example2 clean distclean

all: app

# The application is built directly from source, so `make` does not require ar.
$(APP_BIN): $(APP_SRCS) $(HEADERS)
	$(CXX) $(APP_CPPFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(APP_SRCS) \
		$(HDF5_LDFLAGS) $(HDF5_LIBS) -o $@

app: $(APP_BIN)

# Make examples link the archive built by Make itself. CMake writes its own
# archive to lib/, so mixing build systems cannot make examples link objects
# compiled by a different C++ runtime.
$(BUILD_LIB): $(OBJS)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $^

# Optional static library for external callers.
$(LIB): $(BUILD_LIB)
	@mkdir -p $(dir $@)
	cp $< $@

$(BUILD_DIR)/%.o: src/%.cc $(HEADERS)
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

lib: $(LIB)

$(EXAMPLE1_BIN): src/examples/example1.cc $(BUILD_LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD_LIB) -o $@

example1: $(EXAMPLE1_BIN)
	$(EXAMPLE1_BIN)

$(EXAMPLE2_BIN): src/examples/example2.cc $(BUILD_LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD_LIB) -o $@

example2: $(EXAMPLE2_BIN)
	$(EXAMPLE2_BIN)

clean:
	rm -rf $(BUILD_DIR)

distclean: clean
	rm -f $(LIB) $(APP_BIN)
