add_executable(minimize
  minimize.cpp
)
                
target_link_libraries(minimize
  lapack
  imaging2
)

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_BINARY_DIR})