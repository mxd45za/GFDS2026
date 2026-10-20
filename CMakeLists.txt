cmake_minimum_required(VERSION 3.20)

# This build is for the NVHPC CUDA Fortran path. Set the compiler before
# project() so CMake detects the correct module format and dependency rules.
if(NOT DEFINED CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER nvfortran CACHE FILEPATH "Fortran compiler")
endif()

project(FDS LANGUAGES Fortran)

# Optional SUNDIALS dependency for CVODE detailed chemistry.
option(USE_SUNDIALS "Enable CVODE detailed chemistry through SUNDIALS" OFF)
option(SINGLE_PRECISION "Build with EB and QB in single precision" OFF)

# Source directory.
set(SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/Source")

# Store Fortran module files in the build directory.
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod")

set(FDS_SOURCES
  ${SRC_DIR}/prec.f90
  ${SRC_DIR}/imkl.f90
  ${SRC_DIR}/prop.f90
  ${SRC_DIR}/cons.f90
  ${SRC_DIR}/pois.f90
  ${SRC_DIR}/devc.f90
  ${SRC_DIR}/type.f90
  ${SRC_DIR}/mesh.f90
  ${SRC_DIR}/func.f90
  ${SRC_DIR}/soot.f90
  ${SRC_DIR}/gsmv.f90
  ${SRC_DIR}/rcal.f90
  ${SRC_DIR}/ctrl.f90
  ${SRC_DIR}/turb.f90
  ${SRC_DIR}/data.f90
  ${SRC_DIR}/smvv.f90
  ${SRC_DIR}/geom.f90
  ${SRC_DIR}/vege.f90
  ${SRC_DIR}/init.f90
  ${SRC_DIR}/hvac.f90
  ${SRC_DIR}/ccib.f90
  ${SRC_DIR}/fire.f90
  ${SRC_DIR}/chem.f90
  ${SRC_DIR}/radi.f90
  ${SRC_DIR}/read.f90
  ${SRC_DIR}/dump.f90
  ${SRC_DIR}/pres.f90
  ${SRC_DIR}/wall.f90
  ${SRC_DIR}/velo.f90
  ${SRC_DIR}/part.f90
  ${SRC_DIR}/divg.f90
  ${SRC_DIR}/mass.f90
  ${SRC_DIR}/g_cons.f90
  ${SRC_DIR}/g_erro.f90
  ${SRC_DIR}/g_mesh.f90
  ${SRC_DIR}/g_ht1d.f90
  ${SRC_DIR}/g_dump.f90
  ${SRC_DIR}/g_init.f90
  ${SRC_DIR}/g_func.f90
  ${SRC_DIR}/g_exch.f90
  ${SRC_DIR}/g_fire.f90
  ${SRC_DIR}/g_pois.f90
  ${SRC_DIR}/g_mass.f90
  ${SRC_DIR}/g_wall.f90
  ${SRC_DIR}/g_divg.f90
  ${SRC_DIR}/g_pres.f90
  ${SRC_DIR}/g_velo.f90
  ${SRC_DIR}/g_radi.f90
  ${SRC_DIR}/g_tools.f90
  ${SRC_DIR}/main.f90
)

add_executable(fds ${FDS_SOURCES})
set_target_properties(fds PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)
# NVHPC compiler options.
target_compile_options(fds PRIVATE
  -Mpreprocess
  -mp
  -cuda
  -O3

)

# Intel MKL interfaces.
target_compile_definitions(fds PRIVATE WITH_MKL)

if(SINGLE_PRECISION)
  target_compile_definitions(fds PRIVATE GFDS_SINGLE_PRECISION)
endif()

if(USE_SUNDIALS)
  find_package(SUNDIALS CONFIG REQUIRED)
  target_link_libraries(fds PRIVATE SUNDIALS::fcvode_mod SUNDIALS::fnvecserial_mod)
  target_compile_definitions(fds PRIVATE WITH_SUNDIALS)
  if(DEFINED SUNDIALS_VERSION)
    target_compile_definitions(fds PRIVATE SUNDIALS_PP="${SUNDIALS_VERSION}")
  endif()
endif()

# NVHPC + OpenMPI path. Prefer an explicit -DNVHPC_MPI_ROOT or environment
# variable, then fall back to the OpenMPI bundled with the selected NVHPC SDK.
set(NVHPC_MPI_ROOT "$ENV{NVHPC_MPI_ROOT}" CACHE PATH
  "Root of the OpenMPI installation compatible with NVHPC")
if(NOT NVHPC_MPI_ROOT)
  if(DEFINED ENV{NVHPC_ROOT})
    set(_NVHPC_ROOT "$ENV{NVHPC_ROOT}")
  else()
    get_filename_component(_NVFORTRAN_BIN "${CMAKE_Fortran_COMPILER}" DIRECTORY)
    get_filename_component(_NVHPC_COMPILERS "${_NVFORTRAN_BIN}" DIRECTORY)
    get_filename_component(_NVHPC_ROOT "${_NVHPC_COMPILERS}" DIRECTORY)
  endif()
  file(GLOB _NVHPC_MPI_CANDIDATES LIST_DIRECTORIES true
    "${_NVHPC_ROOT}/comm_libs/*/openmpi*/openmpi-*")
  if(_NVHPC_MPI_CANDIDATES)
    list(GET _NVHPC_MPI_CANDIDATES 0 NVHPC_MPI_ROOT)
    set(NVHPC_MPI_ROOT "${NVHPC_MPI_ROOT}" CACHE PATH
      "Root of the OpenMPI installation compatible with NVHPC" FORCE)
  endif()
endif()
if(NOT EXISTS "${NVHPC_MPI_ROOT}/include" OR NOT EXISTS "${NVHPC_MPI_ROOT}/lib")
  message(FATAL_ERROR
    "Set NVHPC_MPI_ROOT to the compatible OpenMPI root containing include/ and lib/.\n"
    "Example: cmake -DNVHPC_MPI_ROOT=/path/to/openmpi ...")
endif()

target_include_directories(fds PRIVATE
  ${NVHPC_MPI_ROOT}/include
  ${NVHPC_MPI_ROOT}/lib
)

target_link_directories(fds PRIVATE ${NVHPC_MPI_ROOT}/lib)

# MPI Fortran libraries.
target_link_libraries(fds PRIVATE
  mpi_usempif08
  mpi_usempi_ignore_tkr
  mpi_mpifh
  mpi
)

# Runtime library paths relative to the executable.
set(RPATH_FLAGS "-Wl,-rpath,$ORIGIN:$ORIGIN/../../lib:$ORIGIN/../../../../../compilers/lib:$ORIGIN/../../../../../../compilers/lib")

target_link_options(fds PRIVATE
  ${RPATH_FLAGS}
  -Wl,-rpath,${NVHPC_MPI_ROOT}/lib
  -cuda
  -Mmkl
)
