cmake_minimum_required(VERSION 3.12)
project(FIMD C CXX ASM)

# add subprojects FIMD-CPU and FIMD-GPU
add_subdirectory(cpu)
add_subdirectory(gpu)

# include header files "fimd_cpu.h" and "fimd_gpu.h"
include_directories(cpu gpu)

# add dummy target "FIMD"
add_custom_target(${PROJECT_NAME})

# this target compiles shared library files "libfimd_cpu.so" and "libfimd_gpu.so"
add_dependencies(${PROJECT_NAME} fimd_cpu fimd_gpu)

# output shared library absolute paths
file(REAL_PATH "${PROJECT_SOURCE_DIR}/lib/libfimd_cpu.so" FIMD_CPU_LIB)
file(REAL_PATH "${PROJECT_SOURCE_DIR}/lib/libfimd_gpu.so" FIMD_GPU_LIB)

get_directory_property(NESTED_PROJECT PARENT_DIRECTORY)
if(NESTED_PROJECT)
    set(FIMD_CPU_LIB ${FIMD_CPU_LIB} PARENT_SCOPE)
    set(FIMD_GPU_LIB ${FIMD_GPU_LIB} PARENT_SCOPE)
endif()

message("[FIMD shared libraries]")
message("-- FIMD-CPU library: ${FIMD_CPU_LIB}")
message("-- FIMD-GPU library: ${FIMD_GPU_LIB}")

message("[FIMD CMakeLists completed]")
