cmake_minimum_required(VERSION 3.21)
project(flexis CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

include(cmake/downloadDataset.cmake)
include(cmake/third_party_bliss.cmake)

set(third_party "${CMAKE_CURRENT_SOURCE_DIR}/third_party/" CACHE STRING "third_party library directory")
option(INSTRUMENT "instrumentation / profiling support" OFF)

find_package(CLI11 REQUIRED)
find_package(fmt REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph range)
find_package(OpenMP REQUIRED)
find_package(indicators REQUIRED)
find_package(nlohmann_json REQUIRED)

add_subdirectory(third_party)
add_subdirectory(src)

set(DATAFOLDER "data")
set(DATADIR "${CMAKE_CURRENT_SOURCE_DIR}/${DATAFOLDER}")
set(BINARYDATADIR "${CMAKE_CURRENT_BINARY_DIR}/${DATAFOLDER}")
set(BINDIR "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(BINDATADIR "${CMAKE_CURRENT_BINARY_DIR}/bin/${DATAFOLDER}")

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${DATADIR} ${BINARYDATADIR})
file(MAKE_DIRECTORY ${BINDIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${DATADIR} ${BINDATADIR})

