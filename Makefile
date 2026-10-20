# Simple build for the standalone `sayram2d` application.
#
#   make          build ./sayram2d
#   make app      same
#   make app-run  run ./sayram2d with app.ini
#   make lib      build the static library (requires ar)
#   make example  build/run the small library example
#
# Edit the four paths below if your installation is in a different place.
# Tests and benchmarks have separate Makefiles:
#   make -f Makefile.tests
#   make -f Makefile.benchmark

CXX := g++-16
AR ?= ar

EIGEN_INCLUDE ?= /Users/xtao/local/include
HDF5_INCLUDE ?= /opt/homebrew/opt/hdf5/include
HDF5_LIBDIR ?= /opt/homebrew/opt/hdf5/lib
HDF5_LIBS ?= -lhdf5

CPPFLAGS := -Isrc/include -I$(EIGEN_INCLUDE)
APP_CPPFLAGS := -Isrc/app -I$(HDF5_INCLUDE)
CXXFLAGS ?= -std=c++20 -O2 -DEIGEN_NO_DEBUG -Wall -Wextra -Wno-class-memaccess
HDF5_LDFLAGS := -L$(HDF5_LIBDIR)

BUILD_DIR ?= build/sayram2d
LIB := $(BUILD_DIR)/libsayram2d.a
EXAMPLE_BIN := $(BUILD_DIR)/example
APP_BIN := sayram2d

CORE_SRCS := src/mesh.cpp src/boundary_conditions.cpp src/solver.cpp
APP_SRCS := src/app/main.cc src/app/config.cc src/app/io.cc \
            src/app/albert_young.cc $(CORE_SRCS)
OBJS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(CORE_SRCS))

.PHONY: all app app-run lib example clean

all: app

# The application is built directly from source, so `make` does not require ar.
$(APP_BIN): $(APP_SRCS)
	$(CXX) $(APP_CPPFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(APP_SRCS) \
		$(HDF5_LDFLAGS) $(HDF5_LIBS) -o $@

app: $(APP_BIN)

app-run: $(APP_BIN)
	./$(APP_BIN) app.ini

# Optional static library for external callers.
$(LIB): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(AR) rcs $@ $^

$(BUILD_DIR)/%.o: src/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

lib: $(LIB)

$(EXAMPLE_BIN): src/main.cc $(LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LIB) -o $@

example: $(EXAMPLE_BIN)
	$(EXAMPLE_BIN)

clean:
	rm -rf $(BUILD_DIR) $(APP_BIN)
