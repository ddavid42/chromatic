add_executable(gauss_simple gauss.cpp ${HEADER_FILES})

target_compile_features(gauss_simple PRIVATE cxx_std_17)

# Add all headers files under the include directory
target_include_directories(gauss_simple PRIVATE include)
