#include <fem/Transform.hpp>

#include <core/utilities.hpp>

namespace imaging
{
  namespace transform_impl
  {
    void unit_normal(const ublas::fixed_vector<float_t, 2> & start, const ublas::fixed_vector<float_t, 2> & end, ublas::fixed_vector<float_t, 2> & out)
    {
      out(0) = end(1) - start(1);
      out(1) = - ( end(0) - start(0) );
      
      float_t norm = norm_2(out);
      if(norm != 0.0)
        out /= norm;
    }
    
    void unit_normal(const ublas::fixed_vector<float_t, 3> & start, const ublas::fixed_vector<float_t, 3> & end_1, const ublas::fixed_vector<float_t, 3> & end_2, ublas::fixed_vector<float_t, 3> & out)
    {
      out(0) = (end_1(1) - start(1)) * (end_2(2) - start(2)) - (end_1(2) - start(2)) * (end_2(1) - start(1));
      out(1) = (end_1(2) - start(2)) * (end_2(0) - start(0)) - (end_1(0) - start(0)) * (end_2(2) - start(2));
      out(2) = (end_1(0) - start(0)) * (end_2(1) - start(1)) - (end_1(1) - start(1)) * (end_2(0) - start(0));
      
      float_t norm = norm_2(out);
      if(norm != 0.0)
        out /= norm;
    }
  }
  size_t Square2dTransform::face_vertex(size_t face, size_t vertex)
  {
    size_t index;

    if( !(face == 3 && vertex == 1) )
      index = face + vertex;
    else
      index = 0;
    
#ifdef DEBUG
    if( face > 3 || vertex > 1)
     throw Exception("Exception: Too large argument 'face' oder 'vertex' in Square2dTransform::face_vertex().");
#endif

    return index;
  }

  ublas::fixed_matrix<float_t, 2, 1> & Square2dTransform::boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_matrix<float_t, 2, 1> & out) const
  {
    switch(face)
    {
    case 0:
      for(size_t i = 0; i < 2; ++i)
        out(i,0) = .5 * (_vertices(1)(i) - _vertices(0)(i));
      break;
    case 1:
      for(size_t i = 0; i < 2; ++i)
        out(i,0) = .5 * (_vertices(2)(i) - _vertices(1)(i));
      break;
    case 2:
      for(size_t i = 0; i < 2; ++i)
        out(i,0) = .5 * (_vertices(3)(i) - _vertices(2)(i));
      break;
    case 3:
      for(size_t i = 0; i < 2; ++i)
        out(i,0) = .5 * (_vertices(0)(i) - _vertices(3)(i));
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' in Square2dTransform::boundary_derivative().");
#endif
      break;
    }


    return out;
  }

  ublas::fixed_vector<float_t, 2> & Square2dTransform::value(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const
  {
    Bilinear2dShapeFunction shape_function;

    out.assign(0.0);

    for(size_t i = 0; i < 4; ++i)
      out += shape_function.value(i, in) * _vertices(i);

    return out;
  }

  ublas::fixed_matrix<float_t, 2, 2> & Square2dTransform::derivative(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 2, 2> & out) const
  {
    Bilinear2dShapeFunction shape_function;

    ublas::fixed_vector<float_t, 2> temp;

    out.assign(0.0);

    for(size_t i = 0; i < 4; ++i)
    {
      shape_function.gradient(i, in, temp);
      ublas::matrix_row< ublas::fixed_matrix<float_t, 2, 2> >(out, 0) += _vertices(i)(0) * temp;
      ublas::matrix_row< ublas::fixed_matrix<float_t, 2, 2> >(out, 1) += _vertices(i)(1) * temp;
    }

    return out;
  }
  
  ublas::fixed_vector<float_t, 2> & Square2dTransform::boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 2> & out) const
  {
    switch(face_index)
    {
    case 0:
      out.assign(in(0), -1.0);
      break;
    case 1:
      out.assign(1.0, in(0));
      break;
    case 2:
      out.assign(-in(0), 1.0);
      break;
    case 3:
      out.assign(-1.0, -in(0));
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face_index' Transform::boundary2element().");
#endif

    }

    return out;
  }
  
  ublas::fixed_vector<float_t, 2> & Square2dTransform::boundary_normal(size_t face, ublas::fixed_vector<float_t, 2> & out) const
  {
    switch(face)
    {
    case 0:
      transform_impl::unit_normal(_vertices(0), _vertices(1), out);
      break;
    case 1:
      transform_impl::unit_normal(_vertices(1), _vertices(2), out);
      break;
    case 2:
      transform_impl::unit_normal(_vertices(2), _vertices(3), out);
      break;
    case 3:
      transform_impl::unit_normal(_vertices(3), _vertices(0), out);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' Transform::boundary_normal().");
#endif

    }
    
    ublas::fixed_vector<float_t, 2> orientation;
    transform_impl::unit_normal(_vertices(0), _vertices(1), orientation);
    
    if(inner_prod(orientation, _vertices(3) - _vertices(0)) > 0)
      out *= -1;

    return out;

  }


  size_t Triangle2dTransform::face_vertex(size_t face, size_t vertex)
  {
    size_t index;

    if( !(face == 2 && vertex == 1) )
      index = face + vertex;
    else
      index = 0;

    return index;
  }

  ublas::fixed_matrix<float_t, 2, 1> & Triangle2dTransform::boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_matrix<float_t, 2, 1> & out) const
  {
    switch(face)
    {
    case 0:
      for(size_t i = 0; i < 2; ++i)
        out(i,0) = 0.5 * (_vertices(1)(i) - _vertices(0)(i));
      break;
    case 1:
      for(size_t i = 0; i < 2; ++i)
        out(i,0) = 0.5 * (_vertices(2)(i) - _vertices(1)(i));
      break;
    case 2:
      for(size_t i = 0; i < 2; ++i)
        out(i,0) = 0.5 * (_vertices(0)(i) - _vertices(2)(i));
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' in Transform<3,2>::boundary_derivative().");
#endif
      break;
    }


    return out;
  }

  ublas::fixed_vector<float_t, 2> & Triangle2dTransform::value(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const
  {
    out = _vertices(0) + in(0) * (_vertices(1) - _vertices(0)) + in(1) * (_vertices(2) - _vertices(0));

    return out;
  }

  ublas::fixed_matrix<float_t, 2, 2> & Triangle2dTransform::derivative(const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 2, 2> & out) const
  {
    out(0,0) = _vertices(1)(0) - _vertices(0)(0);
    out(0,1) = _vertices(2)(0) - _vertices(0)(0);
    out(1,0) = _vertices(1)(1) - _vertices(0)(1);
    out(1,1) = _vertices(2)(1) - _vertices(0)(1);

    return out;
  }
  
  ublas::fixed_vector<float_t, 2> & Triangle2dTransform::boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 2> & out) const
  {
    float_t temp;

    switch(face_index)
    {
    case 0:
      out.assign( (in(0) + 1.0) / 2.0, 0.0);
      break;
    case 1:
      temp = (in(0) + 1.0) / 2.0;
      out.assign( 1 - temp, temp );
      break;
    case 2:
      out.assign(0.0, 1 - (in(0) + 1.0) / 2.0);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face_index' Transform::boundary2element().");
#endif

    }

    return out;
  }
  
  ublas::fixed_vector<float_t, 2> & Triangle2dTransform::boundary_normal(size_t face, ublas::fixed_vector<float_t, 2> & out) const
  {
    switch(face)
    {
    case 0:
      transform_impl::unit_normal(_vertices(0), _vertices(1), out);
      break;
    case 1:
      transform_impl::unit_normal(_vertices(1), _vertices(2), out);
      break;
    case 2:
      transform_impl::unit_normal(_vertices(2), _vertices(0), out);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' Transform::boundary_normal().");
#endif

    }
    
    ublas::fixed_vector<float_t, 2> orientation;
    transform_impl::unit_normal(_vertices(0), _vertices(1), orientation);
    
    if(inner_prod(orientation, _vertices(2) - _vertices(0)) > 0)
      out *= -1;

    return out;

  }


  size_t Tetrahedra3dTransform::face_vertex(size_t face, size_t vertex)
  {
    return face_vertex_matrix[face][vertex];
  }

  ublas::fixed_matrix<float_t, 3, 2> & Tetrahedra3dTransform::boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 3, 2> & out) const
  {
    switch(face)
    {
    case 0:
      for(size_t i = 0; i < 3; ++i)
      {
        out(i,0) = (_vertices(1)(i) - _vertices(0)(i));
        out(i,1) = (_vertices(3)(i) - _vertices(0)(i));
      }
      break;
    case 1:
      for(size_t i = 0; i < 3; ++i)
      {
        out(i,0) = (_vertices(2)(i) - _vertices(1)(i));
	      out(i,1) = (_vertices(3)(i) - _vertices(1)(i));
      } 
      break;
    case 2:
      for(size_t i = 0; i < 3; ++i)
      {
        out(i,0) = (_vertices(2)(i) - _vertices(0)(i)); 
        out(i,1) = (_vertices(3)(i) - _vertices(0)(i));
      }
      break;
    case 3:
      for(size_t i = 0; i < 3; ++i)
      {
        out(i,0) = (_vertices(1)(i) - _vertices(0)(i));
        out(i,1) = (_vertices(2)(i) - _vertices(0)(i));
      }
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' in Transform<3,2>::boundary_derivative().");
#endif
      break;
    }


    return out;
  }

  ublas::fixed_vector<float_t, 3> & Tetrahedra3dTransform::value(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const
  {
    out = _vertices(0) + in(0) * (_vertices(1) - _vertices(0))  
                        + in(1) * (_vertices(2) - _vertices(0))  
                        + in(2) * (_vertices(3) - _vertices(0));
    return out;
  }

  ublas::fixed_matrix<float_t, 3, 3> & Tetrahedra3dTransform::derivative(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_matrix<float_t, 3, 3> & out) const
  {
    out(0,0) = _vertices(1)(0) - _vertices(0)(0);
    out(0,1) = _vertices(2)(0) - _vertices(0)(0);
    out(0,2) = _vertices(3)(0) - _vertices(0)(0);
    out(1,0) = _vertices(1)(1) - _vertices(0)(1);
    out(1,1) = _vertices(2)(1) - _vertices(0)(1);
    out(1,2) = _vertices(3)(1) - _vertices(0)(1);
    out(2,0) = _vertices(1)(2) - _vertices(0)(2);
    out(2,1) = _vertices(2)(2) - _vertices(0)(2);
    out(2,2) = _vertices(3)(2) - _vertices(0)(2);
    
    return out;
  }
  
  ublas::fixed_vector<float_t, 3> & Tetrahedra3dTransform::boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 3> & out) const
  {
    switch(face_index)
    {
    case 0:
      out.assign( in(0), 0.0, in(1) );
      break;
    case 1:
      out.assign( 1 - in(0) - in(1), in(0), in(1));
      break;
    case 2:
      out.assign(0.0, in(0), in(1) );
      break;
    case 3:
      out.assign(in(0), in(1), 0.0);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face_index' Transform::boundary2element().");
#endif

    }

    return out;
  }
  
  ublas::fixed_vector<float_t, 3> & Tetrahedra3dTransform::boundary_normal(size_t face, ublas::fixed_vector<float_t, 3> & out) const
  {
    switch(face)
    {
    case 0:
      transform_impl::unit_normal(_vertices(0), _vertices(1), _vertices(3), out);
      break;
    case 1:
      transform_impl::unit_normal(_vertices(1), _vertices(2), _vertices(3), out);
      break;
    case 2:
      transform_impl::unit_normal(_vertices(0), _vertices(3), _vertices(2), out);
      break;
    case 3:
      transform_impl::unit_normal(_vertices(0), _vertices(2), _vertices(1), out);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' Transform::boundary_normal().");
#endif

    }
    
    ublas::fixed_vector<float_t, 3> orientation;
    transform_impl::unit_normal(_vertices(0), _vertices(1), _vertices(2), orientation);
    
    if(inner_prod(orientation, _vertices(3) - _vertices(0)) < 0)
      out *= -1;

    return out;

  }
  
  const size_t Tetrahedra3dTransform::face_vertex_matrix[4][3] = {{0, 1, 3}, {1, 2, 3}, {0, 2, 3}, {0, 1, 2}};

  size_t Interval1dTransform::face_vertex(size_t face, size_t vertex)
  {
    return face;
  }

  ublas::fixed_matrix<float_t, 1, 0> & Interval1dTransform::boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 0> & in, ublas::fixed_matrix<float_t, 1, 0> & out) const
  {
    return out;
  }

  ublas::fixed_vector<float_t, 1> & Interval1dTransform::value(const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 1> & out) const
  {
    out = _vertices(0) + (1 + in(0)) / 2.0 * (_vertices(1) - _vertices(0));

    return out;
  }

  ublas::fixed_matrix<float_t, 1, 1> & Interval1dTransform::derivative(const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_matrix<float_t, 1, 1> & out) const
  {
    out(0,0) = ( _vertices(1)(0) - _vertices(0)(0) ) / 2.0;

    return out;
  }
  
  ublas::fixed_vector<float_t, 1> & Interval1dTransform::boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 0> & in, ublas::fixed_vector<float_t, 1> & out) const
  {
    switch(face_index)
    {
    case 0:
      out.assign(-1.0);
      break;
    case 1:
      out.assign(1.0);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' Transform::boundary2element().");
#endif

    }

    return out;
  }
  
  ublas::fixed_vector<float_t, 1> & Interval1dTransform::boundary_normal(size_t face, ublas::fixed_vector<float_t, 1> & out) const
  {
    switch(face)
    {
    case 0:
      out.assign(-1.0);
      break;
    case 1:
      out.assign(1.0);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' Transform::boundary_normal().");
#endif

    }

    return out;

  }


  size_t Cube3dTransform::face_vertex(size_t face, size_t vertex)
  {
    return face_vertex_matrix[face][vertex];
  }


  ublas::fixed_matrix<float_t, 3, 2> & Cube3dTransform::boundary_derivative(size_t face, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_matrix<float_t, 3, 2> & out) const
  {
     for(size_t i = 0; i < 3; ++i)
	{
        out(i,0) = .5 * (_vertices(face_vertex_matrix[face][2])(i) - _vertices(face_vertex_matrix[face][0])(i));
        out(i,1) = .5 * (_vertices(face_vertex_matrix[face][3])(i) - _vertices(face_vertex_matrix[face][0])(i));
 	}    

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' in Transform<4,2>::boundary_derivative().");
#endif
    return out;
  }

  ublas::fixed_vector<float_t, 3> & Cube3dTransform::value(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const
  {
    Trilinear3dShapeFunction shape_function;

    out.assign(0.0);

    for(size_t i = 0; i < 8; ++i)
      out += shape_function.value(i, in) * _vertices(i);
    return out;
  }


  ublas::fixed_matrix<float_t, 3, 3> & Cube3dTransform::derivative(const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_matrix<float_t, 3, 3> & out) const
  {
    Trilinear3dShapeFunction shape_function;

    ublas::fixed_vector<float_t, 3> temp;

    out.assign(0.0);

    for(size_t i = 0; i < 8; ++i)
    {
      shape_function.gradient(i, in, temp);
	  for(size_t j=0; j<3; ++j)
      ublas::matrix_row< ublas::fixed_matrix<float_t, 3, 3> >(out, j) += _vertices(i)(j) * temp;
    }

    return out;
  }

  
  ublas::fixed_vector<float_t, 3> & Cube3dTransform::boundary2element(size_t face_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 3> & out) const
  {
    switch(face_index)
    {
    case 0:
      out.assign(in(0), in(1),-1.0);
      break;
    case 1:
      out.assign(in(0), in(1), 1.0);
      break;
    case 2:
      out.assign(-1.0, in(0), in(1));
      break;
    case 3:
      out.assign(1.0, in(0), in(1) );
      break;
	case 4:
      out.assign(in(0), -1.0, in(1));
      break;
    case 5:
      out.assign(in(0), 1.0, in(1));
      break;
	  

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face_index' Transform::boundary2element().");
#endif

    }

    return out;
  }
  
  ublas::fixed_vector<float_t, 3> & Cube3dTransform::boundary_normal(size_t face, ublas::fixed_vector<float_t, 3> & out) const
  {
    switch(face)
    {
    case 0:
      transform_impl::unit_normal(_vertices(0), _vertices(2), _vertices(1), out);
      break;
    case 1:
      transform_impl::unit_normal(_vertices(4), _vertices(5), _vertices(6), out);
      break;
    case 2:
      transform_impl::unit_normal(_vertices(0), _vertices(4), _vertices(2), out);
      break;
    case 3:
      transform_impl::unit_normal(_vertices(1), _vertices(3), _vertices(5), out);
      break;
    case 4:
      transform_impl::unit_normal(_vertices(0), _vertices(1), _vertices(4), out);
      break;
    case 5:
      transform_impl::unit_normal(_vertices(2), _vertices(6), _vertices(3), out);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'face' Transform::boundary_normal().");
#endif

    }
    
    ublas::fixed_vector<float_t, 3> orientation;
    transform_impl::unit_normal(_vertices(0), _vertices(1), _vertices(3), orientation);
    
    if(inner_prod(orientation, _vertices(4) - _vertices(0)) < 0)
      out *= -1;

    return out;

  }

  /*face nummerierung 0 - unten - 1 - oben - 3 rechts - 4 links - 5 vorne - 6 hinten*/
  const size_t Cube3dTransform::face_vertex_matrix[6][4] = {{0, 1, 3, 2}, {4 ,5 ,7 ,6 },{0 , 2, 6, 4 },{1 ,3 ,7 ,5 },{0 ,1 ,5 ,4 }, {2 ,3 ,7 ,6 }};


}


