cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


option(FORTRAN "Compile Fortran components") 

project(imaging2 CXX C Fortran)

if(UNIX)
  include_directories(SYSTEM /usr/include/GL)
  include_directories(SYSTEM /usr/include/libxml2)
endif(UNIX)

include_directories(SYSTEM ${CMAKE_SOURCE_DIR})

set(CORE_OBJECTS
  core/Cmessage.cxx
  core/distribution_utilities.cxx
  core/utilities.cxx
  core/vector_utilities.cxx
  core/matrix_utilities.cxx
  core/cio.cxx
  core/xmlio.cxx
  external/cg+/cgfam.f
  external/cg+/cgsearch.f
  external/gl2ps/gl2ps.c
  external/gpc/gpc.c
  external/triangle/triangle.c 
  external/aip/aip.c 
  external/liblbfgs/lbfgs.c
  external/itpack/dsrc2c.f 
  graphics/DummyGraphics.cxx
  graphics/GraphicsInterface.cxx
  fem/Assembler.cxx
  fem/ElementIntegrator.cxx
  fem/FemKernel.cxx
  fem/Image2Grid.cxx
  fem/Image2Grid_impl.cxx
  fem/ShapeFunction.cxx
  fem/Transform.cxx
  fem/triangle.c
  fem/utilities.cxx
  image/Color.cxx
  image/cio.cxx 
  image/gio.cxx 
  image/GrayValue.cxx
  image/utilities.cxx
  lapack/linear_algebra.cxx
  minimize/CovarianceMatrixAdaptation.cxx
  minimize/Lbfgs.cxx
  minimize/NlCg.cxx
  minimize/SteepestDescent.cxx
  polytope/gio.cxx
  polytope/Polygon.cxx
  polytope/SimplePolygon.cxx
  polytope/xmlio.cxx
  shape/BoundaryDiscretizer.cxx
  shape/BsplineShape.cxx
  shape/Circle.cxx
  shape/gio.cxx 
  shape/mrep/MrepSkeleton2d.cxx
  shape/mrep/MrepModel2d.cxx
  shape/mrep/PolygonModel2d.cxx
  shape/mrep/Position2d.cxx
  shape/mrep/xmlio.cxx
  shape/mrep/gio.cxx
  shape/xmlio.cxx
  shape/ShapeStatistics.cxx
  statistic/LinearPca.cxx
  statistic/utilities.cxx
  solver/CgSolver.cxx
  solver/utilities.cxx
  spline/gio.cxx
  spline/utilities.cxx
  xml/XmlReader.cxx
  xml/XmlWriter.cxx
)

set(SPOOLES_OBJECTS
  solver/BiCgStabSolver.cxx
  solver/LuSolver.cxx
  solver/spooles.c
)

set(OPENGL_OBJECTS
  graphics/Objects.cxx 
  graphics/OpenGlViewer.cxx 
)

add_library(imaging2
  ${CORE_OBJECTS}
  ${SPOOLES_OBJECTS}
  ${OPENGL_OBJECTS}  
)
                
set_target_properties(imaging2 PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(core)
add_subdirectory(fem)
add_subdirectory(image)
add_subdirectory(polytope)
add_subdirectory(minimize)
add_subdirectory(shape)
add_subdirectory(statistic)
add_subdirectory(xml)
