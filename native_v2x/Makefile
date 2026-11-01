# AI_DNA_ANALYZER V2 native engine.  Build tree lives next to the sources (HDD, not the main disk).
CXX      ?= g++
BUILD    ?= build
include $(dir $(lastword $(MAKEFILE_LIST)))../mk/htslib.mk
# Exactness contract: never fuse a*b+c, never fast-math.  Baseline ISA = generic x86-64;
# AVX2/FMA kernels (Phase D) are added with per-function target attributes + runtime dispatch.
EXACT    := -ffp-contract=off -fno-fast-math -fno-associative-math -frounding-math
ARCH     ?=
CXXFLAGS ?= -O3 -std=c++20 -Wall -Wextra -fPIC $(EXACT) $(ARCH) -Iinclude
$(BUILD)/libdnav2_kernels.so: src/capi.cpp include/dnav2/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -shared -o $@ src/capi.cpp
.PHONY: kernels clean
kernels: $(BUILD)/libdnav2_kernels.so
clean:
	rm -f $(BUILD)/libdnav2_kernels.so

# HTS_INC / HTS_LIB / HTS_RPATH are resolved in ../mk/htslib.mk (see scripts/hts_env.sh)
ENGINE_SRC := src/engine.cpp src/extract.cpp src/frame.cpp src/features.cpp
$(BUILD)/dnav2: src/main.cpp $(ENGINE_SRC) include/dnav2/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -pthread -isystem $(HTS_INC) -o $@ src/main.cpp $(ENGINE_SRC) $(HTS_LIB) $(HTS_RPATH)
.PHONY: engine
engine: $(BUILD)/dnav2

$(BUILD)/dnav2_prof: src/main.cpp $(ENGINE_SRC) include/dnav2/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -DDNAV2_PROFILE -pthread -isystem $(HTS_INC) -o $@ src/main.cpp $(ENGINE_SRC) $(HTS_LIB) $(HTS_RPATH)
.PHONY: profile
profile: $(BUILD)/dnav2_prof
$(BUILD)/bamscan: tools/bamscan.cpp
	@mkdir -p $(BUILD)
	$(CXX) -O2 -std=c++20 -isystem $(HTS_INC) -o $@ tools/bamscan.cpp $(HTS_LIB) $(HTS_RPATH)
.PHONY: tools
tools: $(BUILD)/bamscan

# ---- tests ------------------------------------------------------------------------------------
# Unit / equivalence tests that need no big data:  make test
# Data gates (need the HG002 chr20 300x archive):  see scripts/run_gates.sh
.PHONY: test
test: kernels
	V2_KERNEL_LIB=$(abspath $(BUILD))/libdnav2_kernels.so python3 tests/test_kernels.py
	python3 tests/test_blake2b.py

$(BUILD)/test_frame: tests/test_frame.cpp src/frame.cpp include/dnav2/frame.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -isystem $(HTS_INC) -o $@ tests/test_frame.cpp src/frame.cpp $(HTS_LIB) $(HTS_RPATH)
.PHONY: test-frame
test-frame: $(BUILD)/test_frame
	$(BUILD)/test_frame
