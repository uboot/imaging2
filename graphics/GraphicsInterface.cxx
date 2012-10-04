#include <graphics/GraphicsInterface.hpp>

#include <core/utilities.hpp>
#include <core/matrix_utilities.hpp>

namespace imaging
{
  const GraphicsInterface::FlushCommand GraphicsInterface::flush;
  const GraphicsInterface::ResetCommand GraphicsInterface::reset_group;
  
  void GraphicsInterface::line_segment(const ublas::fixed_vector<float_t, 2> & v_0, const ublas::fixed_vector<float_t, 2> & v_1)
  {
    std::vector< ublas::fixed_vector<float_t, 2> > vertices(2);

    vertices[0] = v_0;
    vertices[1] = v_1;

    polyline(vertices);
  }
    
  void GraphicsInterface::arrow(const ublas::fixed_vector<float_t, 2> & starting_point, const ublas::fixed_vector<float_t, 2> & direction)
  {
    ublas::fixed_vector<float_t, 2> head = starting_point + direction;
    line_segment(starting_point, head);
    
    ublas::fixed_vector<float_t, 2> head_segment = prod(rotation_matrix(PI / 4.0), 0.3 * direction);
    line_segment(head, head - head_segment);
    head_segment = prod(rotation_matrix(- PI / 4.0), 0.3 * direction);
    line_segment(head, head - head_segment);
  }
    
  void GraphicsInterface::quadrangle(const ublas::fixed_vector<float_t, 2> & v_0, const ublas::fixed_vector<float_t, 2> & v_1, const ublas::fixed_vector<float_t, 2> & v_2, const ublas::fixed_vector<float_t, 2> & v_3)
  {
    std::vector< ublas::fixed_vector<float_t, 2> > vertices(4);

    vertices[0] = v_0;
    vertices[1] = v_1;
    vertices[2] = v_2;
    vertices[3] = v_3;

    polygon(vertices);
  }
    
  void GraphicsInterface::triangle(const ublas::fixed_vector<float_t, 2> & v_0, const ublas::fixed_vector<float_t, 2> & v_1, const ublas::fixed_vector<float_t, 2> & v_2)
  {    
    std::vector< ublas::fixed_vector<float_t, 2> > vertices(3);

    vertices[0] = v_0;
    vertices[1] = v_1;
    vertices[2] = v_2;

    polygon(vertices);
  }
}
