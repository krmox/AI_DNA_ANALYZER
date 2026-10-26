# AI_DNA_ANALYZER V2 native engine.  Build tree lives next to the sources (HDD, not the main disk).
CXX      ?= g++
BUILD    ?= build
HTS_INC  := $(shell python3 -c 'import pysam; print(pysam.get_include()[1])')
# Exactness contract: never fuse a*b+c, never fast-math.  Baseline ISA = generic x86-64;
# AVX2/FMA kernels (Phase D) are added with per-function target attributes + runtime dispatch.
EXACT    := -ffp-contract=off -fno-fast-math -fno-associative-math -frounding-math
ARCH     ?=
CXXFLAGS ?= -O3 -std=c++20 -Wall -Wextra -fPIC $(EXACT) $(ARCH) -Iinclude -Itail
$(BUILD)/libdnav2_kernels.so: src/capi.cpp include/dnav2/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -shared -o $@ src/capi.cpp
.PHONY: kernels clean
kernels: $(BUILD)/libdnav2_kernels.so
clean:
	rm -f $(BUILD)/libdnav2_kernels.so

HTS_LIB  := -L/usr/lib64 -l:libhts.so.3
ENGINE_SRC := src/engine.cpp src/extract.cpp src/frame.cpp src/features.cpp tail/tail.cpp tail/calibrate.cpp
$(BUILD)/dnav2: src/main.cpp $(ENGINE_SRC) include/dnav2/*.hpp tail/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -pthread -isystem $(HTS_INC) -o $@ src/main.cpp $(ENGINE_SRC) $(HTS_LIB) -Wl,-rpath,/usr/lib64
.PHONY: engine
engine: $(BUILD)/dnav2

$(BUILD)/dnav2_prof: src/main.cpp $(ENGINE_SRC) include/dnav2/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -DDNAV2_PROFILE -pthread -isystem $(HTS_INC) -o $@ src/main.cpp $(ENGINE_SRC) $(HTS_LIB) -Wl,-rpath,/usr/lib64
.PHONY: profile
profile: $(BUILD)/dnav2_prof
$(BUILD)/bamscan: tools/bamscan.cpp
	@mkdir -p $(BUILD)
	$(CXX) -O2 -std=c++20 -isystem $(HTS_INC) -o $@ tools/bamscan.cpp $(HTS_LIB) -Wl,-rpath,/usr/lib64
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
	$(CXX) $(CXXFLAGS) -isystem $(HTS_INC) -o $@ tests/test_frame.cpp src/frame.cpp $(HTS_LIB) -Wl,-rpath,/usr/lib64
.PHONY: test-frame
test-frame: $(BUILD)/test_frame
	$(BUILD)/test_frame

# ---- R2 tail: unit tests and hardened builds ------------------------------------------------------------
STRICT   := -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wformat=2 -Wnull-dereference -Wdouble-promotion
$(BUILD)/test_tail: tests/test_tail.cpp tail/tail.cpp tail/calibrate.cpp src/frame.cpp include/dnav2/*.hpp tail/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) $(CXXFLAGS) -isystem $(HTS_INC) -o $@ tests/test_tail.cpp tail/tail.cpp tail/calibrate.cpp src/frame.cpp src/extract.cpp src/features.cpp $(HTS_LIB) -Wl,-rpath,/usr/lib64
.PHONY: test-tail
test-tail: $(BUILD)/test_tail
	$(BUILD)/test_tail

# strict-warning compile of the new code only (the inherited R1.1 sources keep their original warning level)
.PHONY: strict
strict:
	$(CXX) -std=c++20 -O2 $(STRICT) -Iinclude -Itail -isystem $(HTS_INC) -fsyntax-only tail/tail.cpp tail/calibrate.cpp

# sanitizer builds of the whole engine: make asan / make tsan   (outputs build/dnav2_asan, build/dnav2_tsan)
$(BUILD)/dnav2_asan: src/main.cpp $(ENGINE_SRC) include/dnav2/*.hpp tail/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) -std=c++20 -O1 -g -fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize-undefined-trap-on-error $(EXACT) -Iinclude -Itail -pthread -isystem $(HTS_INC) -o $@ src/main.cpp $(ENGINE_SRC) $(HTS_LIB) -Wl,-rpath,/usr/lib64
$(BUILD)/dnav2_tsan: src/main.cpp $(ENGINE_SRC) include/dnav2/*.hpp tail/*.hpp
	@mkdir -p $(BUILD)
	$(CXX) -std=c++20 -O1 -g -fsanitize=thread $(EXACT) -Iinclude -Itail -pthread -isystem $(HTS_INC) -o $@ src/main.cpp $(ENGINE_SRC) $(HTS_LIB) -Wl,-rpath,/usr/lib64
.PHONY: asan tsan
asan: $(BUILD)/dnav2_asan
tsan: $(BUILD)/dnav2_tsan
