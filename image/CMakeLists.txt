add_executable(image
  image.cpp
)

target_link_libraries(image
  imaging2
  lapack
  glut
  Magick++
)

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${CMAKE_BINARY_DIR})
