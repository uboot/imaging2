add_executable(xml
  xml.cpp
)
                
target_link_libraries(xml
  xml2  
  imaging2
)

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_BINARY_DIR})
