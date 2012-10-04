add_executable(statistic
  statistic.cpp
)
                
target_link_libraries(statistic
  lapack
  imaging2
)

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_BINARY_DIR})
