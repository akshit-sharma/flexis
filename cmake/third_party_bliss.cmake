include(FetchContent)

set(BLISS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/bliss)

FetchContent_Declare(
  bliss
  URL https://users.aalto.fi/~tjunttil/bliss/downloads/bliss-0.77.zip
  SOURCE_DIR ${BLISS_SOURCE_DIR}
  DOWNLOAD_EXTRACT_TIMESTAMP YES
)

FetchContent_GetProperties(bliss)
if (NOT bliss_POPULATED)
  FetchContent_Populate(bliss)

  message(STATUS "Patching bliss at ${bliss_SOURCE_DIR}")

  # Use sed to change 'protected:' to 'public:' at line 67
  execute_process(
    COMMAND sed -i "1s/cmake_minimum_required(VERSION 3.5)/cmake_minimum_required(VERSION 3.19)/"
    "${bliss_SOURCE_DIR}/CMakeLists.txt"
  )
  # execute_process(
  #   COMMAND sed -i "67s/protected:/public:/" "${bliss_SOURCE_DIR}/src/digraph.hh"
  # )

  # Create the include directory inside the bliss source directory
  file(MAKE_DIRECTORY ${bliss_SOURCE_DIR}/include)

  # Check if the symlink already exists
  set(SYMLINK_PATH "${bliss_SOURCE_DIR}/include/bliss")

  file(MAKE_DIRECTORY ${bliss_SOURCE_DIR}/include)
  set(SYMLINK_PATH "${bliss_SOURCE_DIR}/include/bliss")
  if(NOT EXISTS ${SYMLINK_PATH})
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E create_symlink ../src ${SYMLINK_PATH}
      WORKING_DIRECTORY ${bliss_SOURCE_DIR}
    )
  endif()

  add_subdirectory(${bliss_SOURCE_DIR} ${bliss_BINARY_DIR})
endif()
