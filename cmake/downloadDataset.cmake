
#download some datasets
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data/downloaded/")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data/extracted/")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data/processed/")

function (downloadTxtGz url filename)
  set (downloadedFile "${CMAKE_CURRENT_SOURCE_DIR}/data/downloaded/${filename}.gz")
  set (extractedFile "${CMAKE_CURRENT_SOURCE_DIR}/data/extracted/${filename}")

  if (NOT EXISTS "${downloadedFile}")
    message(STATUS "Downloading ${filename}...")
    add_custom_command(
      OUTPUT "${downloadedFile}"
      COMMAND wget -q --show-progress -O "${downloadedFile}" "${url}"
      COMMENT "Downloading ${filename}..."
    )
  endif()

  add_custom_command(
    OUTPUT "${extractedFile}"
    COMMAND gzip -d -c "${downloadedFile}" > "${extractedFile}"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data/extracted/"
    DEPENDS "${downloadedFile}"
    COMMENT "Extracting ${filename}..."
  )

  add_custom_target(${filename} ALL DEPENDS ${extractedFile})
endfunction()

downloadTxtGz("https://snap.stanford.edu/data/cit-Patents.txt.gz" "cit-Patents.txt")
downloadTxtGz("https://snap.stanford.edu/data/wiki-Vote.txt.gz" "wiki-Vote.txt")

