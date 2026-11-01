# Shared htslib resolution for the native Makefiles.  Include with:
#     include $(dir $(lastword $(MAKEFILE_LIST)))../mk/htslib.mk
# The logic lives in scripts/hts_env.sh (shared with native/build.sh and verify.sh).
# Override any of HTS_PREFIX / HTS_INC / HTS_LIBDIR / HTS_LIB / HTS_RPATH on the make
# command line or in the environment.  `make hts-info` prints what was resolved.
_HTS_SH := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))../scripts/hts_env.sh)
HTS_INC    ?= $(shell bash $(_HTS_SH) inc)
HTS_LIBDIR ?= $(shell bash $(_HTS_SH) libdir)
HTS_LIB    ?= $(shell HTS_INC='$(HTS_INC)' HTS_LIBDIR='$(HTS_LIBDIR)' bash $(_HTS_SH) lib)
HTS_RPATH  ?= -Wl,-rpath,$(HTS_LIBDIR)

.PHONY: hts-info
hts-info:
	@HTS_INC='$(HTS_INC)' HTS_LIBDIR='$(HTS_LIBDIR)' HTS_LIB='$(HTS_LIB)' HTS_RPATH='$(HTS_RPATH)' bash $(_HTS_SH) report

# Fail early, with instructions, for any goal that compiles against htslib.
_HTS_NEEDED := $(filter-out clean kernels hts-info test help,$(MAKECMDGOALS))
ifneq ($(_HTS_NEEDED),)
  ifeq ($(wildcard $(HTS_INC)/htslib/sam.h),)
    $(info htslib headers not found at '$(HTS_INC)/htslib/sam.h'.)
    $(info Run 'make hts-info', or: conda env create -f environment.yml, or set HTS_PREFIX=/path/to/htslib)
    $(error cannot build $(_HTS_NEEDED) without htslib; see docs/PLATFORM_SUPPORT.md)
  endif
endif
