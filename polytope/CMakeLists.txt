add_executable(polytope
  polytope.cpp
)
                
target_link_libraries(polytope
  GL
  GLU
  glut
  pthread
  xml2 
  imaging2
  Magick++
)

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_BINARY_DIR})
