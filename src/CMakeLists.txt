
function(compile targetName cppFiles extraDefs extraLibs compileOptions linkOptions)
  add_executable(${targetName} ${cppFiles})

  set_target_properties(${targetName}
    PROPERTIES
      CXX_STANDARD 20
      CXX_STANDARD_REQUIRED YES
      CXX_EXTENSIONS NO
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
      LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

  target_link_libraries(${targetName}
    PUBLIC
      OpenMP::OpenMP_CXX
      CLI11::CLI11
      fmt::fmt
      vf3lib
      Boost::graph
      Boost::range
      bliss
      nlohmann_json::nlohmann_json
      ${extraLibs}
    )

  target_include_directories(${targetName} SYSTEM
    PUBLIC
      ${BLISS_SOURCE_DIR}/include
      )

  if (INSTRUMENT)
    target_compile_definitions(${targetName} PUBLIC INSTRUMENTATION)
  endif()

  if(extraDefs)
    target_compile_definitions(${targetName} PUBLIC ${extraDefs})
  endif()

  if (compileOptions)
    target_compile_options(${targetName} PUBLIC ${compileOptions})
  endif()

  if (linkOptions)
    target_link_options(${targetName} PUBLIC ${linkOptions})
  endif()

endfunction()

compile(flexis mining.cpp "" "" "" "")
compile(flexisSplit miningTiming.cpp "" "" "" "")
compile(process process.cpp "" "" "" "")
