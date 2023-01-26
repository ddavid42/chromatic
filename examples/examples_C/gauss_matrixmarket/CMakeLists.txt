add_executable(gauss_matrixmarket mmio.cpp gauss.cpp ${HEADER_FILES})

target_compile_features(gauss_matrixmarket PRIVATE cxx_std_17)

# Add all headers files under the include directory
target_include_directories(gauss_matrixmarket PRIVATE include)
