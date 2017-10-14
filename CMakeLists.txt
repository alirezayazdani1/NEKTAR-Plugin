# paraview build location
# /intrepid-fs0/projects/fl/eureka/packages/PV3.14.1/PV_BUILD_3.14.1/paraview/src/paraview-build
# blas location
# /soft/apps/atlas-3.8.2-gcc-shared/lib

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
  CMAKE_POLICY(SET CMP0017 NEW)
ENDIF(COMMAND CMAKE_POLICY)

project(Nektar)

# Try building a test executable to see if all dependencies
# of the plugin library are satisfied.
option(BUILD_TEST_EXE "Build a test executable." OFF)
#mark_as_advanced(BUILD_TEST_EXE)
if(BUILD_TEST_EXE)
  add_definitions(-DBUILD_EXE)
endif()

if(NOT MPI_LIBRARY)
  find_package(MPI REQUIRED)
endif()

enable_language(Fortran)
if(NOT BLAS_LIBRARIES)
  find_package(BLAS REQUIRED)
endif()

#set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
message("the current blas libs are ${BLAS_LIBRARIES}")

if(NOT LAPACK_lapack_LIBRARY)
  find_PACKAGE(LAPACK)
endif()

message("the current lapack lib is ${LAPACK_lapack_LIBRARY}")

try_compile(BLAS_MANGLING_HAS_UNDERSCORE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/linkBLAS.cxx
  CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${BLAS_LIBRARIES}"
  COMPILE_DEFINITIONS "-DBLAS_MANGLING_HAS_UNDERSCORE"
  OUTPUT_VARIABLE compileResult
)

try_compile(NORMAL_BLAS_CALLS
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/linkBLAS.cxx
  CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${BLAS_LIBRARIES}"
  COMPILE_DEFINITIONS "-DNO_BLAS_MANGLING"
  OUTPUT_VARIABLE compileResult
)

if(BLAS_MANGLING_HAS_UNDERSCORE)
  message("calling blas routines from C/C++ requires underscore")
  add_definitions(-DBLAS_MANGLING_HAS_UNDERSCORE)
elseif(NORMAL_BLAS_CALLS)
  message("calling blas routines from C/C++ without modifications")
  add_definitions(-DNO_BLAS_MANGLING)
else()
  #message(FATAL_ERROR "don't know how to call blas routines from C/C++")
  message("don't know how to call blas routines from C/C++")
endif()

add_subdirectory(Veclib)
add_subdirectory(metis)
add_subdirectory(Hlib)
add_subdirectory(Utilities)
add_subdirectory(NektarObj)

# ParaView
if(NOT ParaView_SOURCE_DIR)
  find_package(ParaView 3.14 REQUIRED)
  message("package version is " ${PARAVIEW_VERSION_FULL})
  if(${PARAVIEW_VERSION_MAJOR} EQUAL 3 AND ${PARAVIEW_VERSION_MINOR} EQUAL 14)
    message("using pv 3.14")
    add_definitions(-DVTK5)
  else()
    message("using something after pv 3.14")
  endif()
  include(${PARAVIEW_USE_FILE})
endif()

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/NektarObj
  ${MPI_CXX_INCLUDE_PATH}
  )

add_paraview_plugin(pvNektarReader "2.0"
   SERVER_MANAGER_XML pvNektarReader.xml
   SERVER_MANAGER_SOURCES vtkNektarReader.cxx
   REQUIRED_ON_SERVER)

target_link_libraries(pvNektarReader LINK_PUBLIC
  nektar
  hybrid
  vec
  metis
  nektar
  nektarobj
  ${BLAS_LIBRARIES}
  ${LAPACK_lapack_LIBRARY}
  ${BLAS_LIBRARIES}
  )

if(BUILD_TEST_EXE)
  add_executable(TestReader TestReader.cxx)
  target_link_libraries(TestReader
    pvNektarReader
    ${LAPACK_lapack_LIBRARY}
    ${BLAS_LIBRARIES}
    )
endif()

