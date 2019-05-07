cmake_minimum_required(VERSION 3.1)

option(SRC_BUILD_OPENMP
  "Build the OpenMP configuration" ON)

option(SRC_BUILD_MPI
  "Build the MPI configuration (If enabled, overrides OpenMP option)" OFF)

if (SRC_BUILD_OPENMP AND (NOT SRC_BUILD_MPI))
  
  find_package(OpenMP REQUIRED)
  
  add_executable(openmp main.cxx)
  
  # based on https://cliutils.gitlab.io/modern-cmake/chapters/packages/OpenMP.html
  if(NOT TARGET OpenMP::OpenMP_CXX)
    # CMake < 3.9
    find_package(Threads REQUIRED)

    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)

    set_property(TARGET OpenMP::OpenMP_CXX
      PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})

    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
      PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

    # set OpenMP_CXX_FLAGS
    if (OPENMP_FOUND AND DEFINED OpenMP_CXX_FLAGS)
      set(OpenMP_CXX_FOUND OPENMP_FOUND)
    endif()
  endif()

  if(OpenMP_CXX_FOUND)
    add_definitions(-DUSE_OPENMP)
    target_link_libraries(openmp PUBLIC OpenMP::OpenMP_CXX)
  endif()
  
endif()

if (SRC_BUILD_MPI)
  
  find_package(MPI REQUIRED)

  add_executable(mpi main.cxx)

  # based on https://cliutils.gitlab.io/modern-cmake/chapters/packages/MPI.html
  if (NOT TARGET MPI::MPI_CXX)
    # CMake < 3.9
    cmake_policy(SET CMP0004 OLD)
    add_library(MPI::MPI_CXX IMPORTED INTERFACE)
    cmake_policy(SET CMP0004 NEW)
    
    set_property(TARGET MPI::MPI_CXX
      PROPERTY INTERFACE_COMPILE_OPTIONS ${MPI_CXX_COMPILE_FLAGS})
    
    set_property(TARGET MPI::MPI_CXX
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_PATH}")
    
    set_property(TARGET MPI::MPI_CXX
      PROPERTY INTERFACE_LINK_LIBRARIES ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
    
  endif()
  
  if(MPI_CXX_FOUND)
    add_definitions(-DUSE_MPI)
    target_link_libraries(mpi PUBLIC MPI::MPI_CXX)
  endif()
  
endif()
