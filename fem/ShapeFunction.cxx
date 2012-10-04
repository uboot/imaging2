#include <fem/ShapeFunction.hpp>

#include <fem/Transform.hpp>


namespace imaging
{
  size_t Bilinear2dShapeFunction::face_node(size_t face, size_t node)
  {
    return Square2dTransform::face_vertex(face, node);
  }
  
  float_t Bilinear2dShapeFunction::value(size_t node, const ublas::fixed_vector<float_t, 2> & in) const
  {
    switch(node)
    {
    case 0:
      return float_t( 1.0/4.0 * (1.0-in(0)) * (1.0-in(1)) );
    case 1:
      return float_t( 1.0/4.0 * (1.0+in(0)) * (1.0-in(1)) );
    case 2:
      return float_t( 1.0/4.0 * (1.0+in(0)) * (1.0+in(1)) );
    case 3:
      return float_t( 1.0/4.0 * (1.0-in(0)) * (1.0+in(1)) );

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::value().");
#endif

    }
    return 0.0;
  }
  
  ublas::fixed_vector<float_t, 2> & Bilinear2dShapeFunction::gradient(size_t node_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const
  {
    switch(node_index)
    {
    case 0:
      out.assign(-(1.0-in(1)), -(1.0-in(0)));
      break;
    case 1:
      out.assign(1.0-in(1), -(1.0+in(0)));
      break;
    case 2:
      out.assign(1.0+in(1), 1.0+in(0));
      break;
    case 3:
      out.assign(-(1.0+in(1)),1.0-in(0));
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::gradient().");
#endif

    }

    out *= 1.0/4.0;

    return out;
  }
  
  size_t Trilinear3dShapeFunction::face_node(size_t face, size_t node)
  {
    return Cube3dTransform::face_vertex(face, node);
  }
  
  float_t Trilinear3dShapeFunction::value(size_t node, const ublas::fixed_vector<float_t, 3> & in) const
  {
    switch(node)
    {
    case 0:
      return 1.0/8.0 * (1.0-in(0)) * (1.0-in(1)) *(1.0-in(2));
    case 1:
      return 1.0/8.0 * (1.0+in(0)) * (1.0-in(1)) *(1.0-in(2));
    case 2:
      return 1.0/8.0 * (1.0-in(0)) * (1.0+in(1)) *(1.0-in(2));
    case 3:
      return 1.0/8.0 * (1.0+in(0)) * (1.0+in(1)) *(1.0-in(2));
    case 4:
      return 1.0/8.0 * (1.0-in(0)) * (1.0-in(1)) *(1.0+in(2));
    case 5:
      return 1.0/8.0 * (1.0+in(0)) * (1.0-in(1)) *(1.0+in(2));
    case 6:
      return 1.0/8.0 * (1.0-in(0)) * (1.0+in(1)) *(1.0+in(2));
    case 7:
      return 1.0/8.0 * (1.0+in(0)) * (1.0+in(1)) *(1.0+in(2));

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Trilinear3dShapeFunction::value().");
#endif

    }
    return 0.0;
  }

  ublas::fixed_vector<float_t, 3> & Trilinear3dShapeFunction::gradient(size_t node_index, const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const
  {
    switch(node_index)
    {

    case 0:
      out.assign( - (1.0-in(1)) * (1.0-in(2)),  - (1.0-in(0)) * (1.0-in(2)), - (1.0-in(0)) * (1.0-in(1)) );
   break;
    case 1:
      out.assign( (1.0-in(1)) *(1.0-in(2)), -(1.0+in(0)) *(1.0-in(2)), - (1.0+in(0)) *(1.0-in(1)) );
   break;
    case 2:
      out.assign( - (1.0+in(1)) * (1.0-in(2)), (1.0-in(0)) * (1.0-in(2)), - (1.0-in(0)) * (1.0+in(1)) );
   break;
    case 3:
      out.assign(  (1.0+in(1)) *(1.0-in(2)),(1.0+in(0)) *(1.0-in(2)),- (1.0+in(0)) *(1.0+in(1)) );
   break;
  case 4:
      out.assign(- (1.0-in(1)) * (1.0+in(2)),- (1.0-in(0)) * (1.0+in(2)),(1.0-in(0)) * (1.0-in(1)) );
   break;
    case 5:
      out.assign( (1.0-in(1)) * (1.0+in(2)),- (1.0+in(0)) * (1.0+in(2)),(1.0+in(0)) * (1.0-in(1)) );
   break;
    case 6:
      out.assign(-(1.0+in(1)) *(1.0+in(2)),(1.0-in(0)) *(1.0+in(2)),(1.0-in(0)) *(1.0+in(1)) );
   break;
  case 7:
      out.assign((1.0+in(1)) *(1.0+in(2)),(1.0+in(0)) *(1.0+in(2)), (1.0+in(0)) *(1.0+in(1)));
        break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Trilinear3dShapeFunction::gradient().");
#endif

    }

    out *= 1.0/8.0;

    return out;
  }

  size_t Linear2dShapeFunction::face_node(size_t face, size_t node)
  {
    return Triangle2dTransform::face_vertex(face, node);
  }
  
  float_t Linear2dShapeFunction::value(size_t node_index, const ublas::fixed_vector<float_t, 2> & in) const
  {
    switch(node_index)
    {
    case 0:
      return float_t( 1.0 - in(0) - in(1) );
    case 1:
      return float_t( in(0) );
    case 2:
      return float_t( in(1) );

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::value().");
#endif

    }

    return 0.0;
  }

  ublas::fixed_vector<float_t, 2> & Linear2dShapeFunction::gradient(size_t node_index, const ublas::fixed_vector<float_t, 2> & in, ublas::fixed_vector<float_t, 2> & out) const
  {
    switch(node_index)
    {
    case 0:
      out.assign(-1.0, -1.0);
      break;
    case 1:
      out.assign(1.0, 0.0);
      break;
    case 2:
      out.assign(0.0, 1.0);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::gradient().");
#endif

    }

    return out;
  }
  
  size_t Linear3dShapeFunction::face_node(size_t face, size_t node)
  {
    return Tetrahedra3dTransform::face_vertex(face, node);
  }
  
  float_t Linear3dShapeFunction::value(size_t node_index, const ublas::fixed_vector<float_t, 3> & in) const
  {
    switch(node_index)
    {
    case 0:
      return float_t( 1.0 - in(0) - in(1) - in(2) );
    case 1:
      return float_t( in(0) );
    case 2:
      return float_t( in(1) );
    case 3:
return float_t( in(2) );

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::value().");
#endif

    }

    return 0.0;
  }

  ublas::fixed_vector<float_t, 3> & Linear3dShapeFunction::gradient(size_t node_index, const ublas::fixed_vector<float_t, 3> & in, ublas::fixed_vector<float_t, 3> & out) const
  {
    switch(node_index)
    {
    case 0:
      out.assign(-1.0, -1.0, -1.0);
      break;
    case 1:
      out.assign(1.0, 0.0, 0.0);
      break;
    case 2:
      out.assign(0.0, 1.0, 0.0);
      break;
    case 3:
      out.assign(0.0, 0.0, 1.0);
      break;
#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::gradient().");
#endif

    }

    return out;
  }
  
  size_t Linear1dShapeFunction::face_node(size_t face, size_t node)
  {
    return Interval1dTransform::face_vertex(face, node);
  }
  
  float_t Linear1dShapeFunction::value(size_t node, const ublas::fixed_vector<float_t, 1> & in) const
  {
    switch(node)
    {
    case 0:
      return 1.0 - (in(0) + 1.0) / 2.0;
    case 1:
      return 0.0 + (in(0) + 1.0) / 2.0;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::value().");
#endif

    }

    return 0.0;
  }

  ublas::fixed_vector<float_t, 1> & Linear1dShapeFunction::gradient(size_t node, const ublas::fixed_vector<float_t, 1> & in, ublas::fixed_vector<float_t, 1> & out) const
  {
    switch(node)
    {
    case 0:
      out.assign(-0.5);
      break;
    case 1:
      out.assign(0.5);
      break;

#ifdef DEBUG
    default:
      throw Exception("Exception: Too large argument 'node' Bilinear2dShapeFunction::gradient().");
#endif

    }

    return out;
  }
}

