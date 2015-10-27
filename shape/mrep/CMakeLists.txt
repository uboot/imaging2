add_executable(mrep
  mrep.cpp
)
                
target_link_libraries(mrep
  GL
  GLU
  glut
  pthread
  png 
  xml2 
  imaging2
  Magick++ 
)

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_BINARY_DIR})
