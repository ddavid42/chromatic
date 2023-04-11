add_executable(gauss_matrixmarket mmio.cpp gauss.cpp ${HEADER_FILES})

target_compile_features(gauss_matrixmarket PRIVATE cxx_std_17)

# Add all headers files under the include directory
if (NOT WITHOUT_ORIGINS)
  target_include_directories(gauss_matrixmarket PRIVATE ../../../include)
  target_compile_definitions(gauss_matrixmarket PRIVATE _ORIGINS _ORIGINS_ERROR)
  target_compile_options(gauss_matrixmarket PRIVATE -include "origins.hpp")
  set(COMPILATION_OPTIONS "-D_ORIGINS_MAIN")
  if (WITH_MATRIX)
    string(APPEND COMPILATION_OPTIONS " -D_MATRIX")
  endif()
  if (WITH_INPUT_VECTOR)
    string(APPEND COMPILATION_OPTIONS " -D_INPUT_VECTOR")
  endif()
  set_source_files_properties(gauss.cpp PROPERTIES COMPILE_FLAGS ${COMPILATION_OPTIONS})
endif()

