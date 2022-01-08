cmake_minimum_required(VERSION 3.10)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "cmake")

project(simple-2d-constraint-solver)

set(CMAKE_CXX_STANDARD 11)

# ========================================================
# GTEST

include(FetchContent)
FetchContent_Declare(
    googletest
    URL
    https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set_property(TARGET gmock PROPERTY FOLDER "gtest")
set_property(TARGET gmock_main PROPERTY FOLDER "gtest")
set_property(TARGET gtest PROPERTY FOLDER "gtest")
set_property(TARGET gtest_main PROPERTY FOLDER "gtest")

# ========================================================

add_library(simple-2d-constraint-solver STATIC
    # Source files
    src/rigid_body_system.cpp
    src/matrix.cpp

    # Header files
    include/rigid_body_system.h
    include/matrix.h
)

# GTEST

enable_testing()

add_executable(simple-2d-constraint-solver-test
    test/sanity_tests.cpp
    test/matrix_tests.cpp
)

target_link_libraries(simple-2d-constraint-solver-test
    gtest_main
    simple-2d-constraint-solver
)

include(GoogleTest)
gtest_discover_tests(simple-2d-constraint-solver-test)
