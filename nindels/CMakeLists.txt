cmake_minimum_required(VERSION 3.16)
project(nindels CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
add_compile_options(-Wall -Wextra -Wpedantic)

# Core algorithm: no external dependencies.
add_library(nindels_core STATIC
  src/types.cpp
  src/pattern.cpp
  src/similarity.cpp
  src/memory.cpp
  src/evidence.cpp
  src/likelihood.cpp
  src/reassess.cpp
  src/model.cpp
  src/engine.cpp
  src/table_io.cpp
  src/vcf_out.cpp)
target_include_directories(nindels_core PUBLIC include)

enable_testing()
add_executable(nindels-tests tests/test_core.cpp)
target_link_libraries(nindels-tests nindels_core)
add_test(NAME core COMMAND nindels-tests)

# BAM/VCF front end: htslib. Point HTSLIB_ROOT at an install prefix, or rely
# on pkg-config. Without htslib only the table-based subcommands are built.
find_package(PkgConfig QUIET)
if(HTSLIB_ROOT)
  find_path(HTSLIB_INCLUDE_DIR htslib/sam.h PATHS ${HTSLIB_ROOT}/include NO_DEFAULT_PATH)
  find_library(HTSLIB_LIBRARY hts PATHS ${HTSLIB_ROOT}/lib NO_DEFAULT_PATH)
elseif(PkgConfig_FOUND)
  pkg_check_modules(HTS QUIET htslib)
  if(HTS_FOUND)
    find_path(HTSLIB_INCLUDE_DIR htslib/sam.h HINTS ${HTS_INCLUDE_DIRS})
    find_library(HTSLIB_LIBRARY hts HINTS ${HTS_LIBRARY_DIRS})
  endif()
endif()

find_package(Threads REQUIRED)
set(CLI_SOURCES src/main.cpp)
if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  message(STATUS "htslib: ${HTSLIB_LIBRARY}")
  list(APPEND CLI_SOURCES src/bam_scan.cpp src/truth.cpp)
  add_executable(nindels-call ${CLI_SOURCES})
  target_compile_definitions(nindels-call PRIVATE NINDELS_HAVE_HTSLIB=1)
  target_include_directories(nindels-call PRIVATE ${HTSLIB_INCLUDE_DIR})
  target_link_libraries(nindels-call nindels_core ${HTSLIB_LIBRARY} Threads::Threads)
  if(HTSLIB_ROOT)
    set_target_properties(nindels-call PROPERTIES BUILD_RPATH ${HTSLIB_ROOT}/lib INSTALL_RPATH ${HTSLIB_ROOT}/lib)
  endif()
else()
  message(WARNING "htslib not found: nindels-call is built without BAM/VCF subcommands")
  add_executable(nindels-call ${CLI_SOURCES})
  target_link_libraries(nindels-call nindels_core Threads::Threads)
endif()
