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
    src/fixed_rotation_constraint.cpp
    src/rigid_body_system.cpp
    src/matrix.cpp
    src/system_state.cpp
    src/utilities.cpp
    src/ode_solver.cpp
    src/euler_ode_solver.cpp
    src/rk4_ode_solver.cpp
    src/constraint.cpp
    src/force_generator.cpp
    src/rigid_body.cpp
    src/sle_solver.cpp
    src/gauss_seidel_sle_solver.cpp
    src/fixed_position_constraint.cpp
    src/static_force_generator.cpp
    src/gravity_force_generator.cpp
    src/link_constraint.cpp
    src/gaussian_elimination_sle_solver.cpp
    src/spring.cpp
    src/line_constraint.cpp
    src/constant_speed_motor.cpp
    src/rolling_constraint.cpp
    src/sparse_matrix.cpp
    src/constant_rotation_constraint.cpp
    src/nsv_ode_solver.cpp
    src/generic_rigid_body_system.cpp
    src/optimized_nsv_rigid_body_system.cpp
    src/clutch_constraint.cpp
    src/conjugate_gradient_sle_solver.cpp
    src/rotation_friction_constraint.cpp
    src/simple_gear_constraint.cpp

    # Header files
    include/fixed_rotation_constraint.h
    include/rigid_body_system.h
    include/matrix.h
    include/utilities.h
    include/system_state.h
    include/ode_solver.h
    include/euler_ode_solver.h
    include/rk4_ode_solver.h
    include/constraint.h
    include/force_generator.h
    include/rigid_body.h
    include/sle_solver.h
    include/gauss_seidel_sle_solver.h
    include/fixed_position_constraint.h
    include/static_force_generator.h
    include/gravity_force_generator.h
    include/link_constraint.h
    include/gaussian_elimination_sle_solver.h
    include/spring.h
    include/line_constraint.h
    include/constant_speed_motor.h
    include/rolling_constraint.h
    include/sparse_matrix.h
    include/constant_rotation_constraint.h
    include/nsv_ode_solver.h
    include/scs.h
    include/generic_rigid_body_system.h
    include/optimized_nsv_rigid_body_system.h
    include/clutch_constraint.h
    include/conjugate_gradient_sle_solver.h
    include/rotation_friction_constraint.h
    include/simple_gear_constraint.h
)

# GTEST

enable_testing()

add_executable(simple-2d-constraint-solver-test
    # Source files
    test/sanity_tests.cpp
    test/matrix_tests.cpp
    test/euler_ode_solver_tests.cpp
    test/rk4_ode_solver_tests.cpp
    test/gauss_seidel_sle_solver_tests.cpp
    test/gaussian_elimination_sle_solver_tests.cpp
    test/rolling_constraint_tests.cpp
    test/utilities.cpp
    test/sparse_matrix_tests.cpp
    test/conjugate_gradient_sle_solver_tests.cpp

    # Header files
    test/utilities.h
)

target_link_libraries(simple-2d-constraint-solver-test
    gtest_main
    simple-2d-constraint-solver
)

include(GoogleTest)
gtest_discover_tests(simple-2d-constraint-solver-test)
