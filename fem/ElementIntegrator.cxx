#include <fem/ElementIntegrator.hpp>


namespace imaging
{
  template<>
  SquareIntegrator<4>::SquareIntegrator() : ElementIntegrator<2, 4>()
  {
    static const float_t GAUSS_NODE = 0.5773502694;

    _nodes(0).assign(GAUSS_NODE, GAUSS_NODE);
    _nodes(1).assign(-GAUSS_NODE, GAUSS_NODE);
    _nodes(2).assign(-GAUSS_NODE, -GAUSS_NODE);
    _nodes(3).assign(GAUSS_NODE, -GAUSS_NODE);

    _weights.assign(1.0);
  }

  template<>
  TriangleIntegrator<1>::TriangleIntegrator() : ElementIntegrator<2, 1>()
  {
    _nodes(0).assign(1.0/3.0, 1.0/3.0); // barycenter node
    _weights.assign(1.0/2.0);
  }

  template<>
  TriangleIntegrator<4>::TriangleIntegrator() : ElementIntegrator<2, 4>()
  {
    _nodes(0).assign(0.33333333, 0.33333333);
    _nodes(1).assign(0.73333333, 0.13333333);
    _nodes(2).assign(0.13333333, 0.73333333);
    _nodes(3).assign(0.13333333, 0.13333333);

    _weights.assign(0.28125000, 0.26041666, 0.26041666, 0.26041666);
  }

  template<>
  TriangleIntegrator<7>::TriangleIntegrator() : ElementIntegrator<2, 7>()
  {
    _nodes(0).assign(0.33333333, 0.33333333);
    _nodes(1).assign(0.05961587, 0.47014206);
    _nodes(2).assign(0.47014206, 0.05961587);
    _nodes(3).assign(0.47014206, 0.47014206);
    _nodes(4).assign(0.79742699, 0.10128651);
    _nodes(5).assign(0.10128651, 0.79742699);
    _nodes(6).assign(0.10128651, 0.10128651);

    _weights(0) = 0.22500000;
    _weights(1) = 0.13239415;
    _weights(2) = 0.13239415;
    _weights(3) = 0.13239415;
    _weights(4) = 0.12593918;
    _weights(5) = 0.12593918;
    _weights(6) = 0.12593918;
  }

 template<>
  CubeIntegrator<8>::CubeIntegrator() : ElementIntegrator<3, 8>()
  {
    static const float_t GAUSS_NODE = 0.5773502694;

    _nodes(0).assign(GAUSS_NODE, GAUSS_NODE, GAUSS_NODE);
    _nodes(1).assign(-GAUSS_NODE, GAUSS_NODE, GAUSS_NODE);
    _nodes(2).assign(GAUSS_NODE, -GAUSS_NODE, GAUSS_NODE);
    _nodes(3).assign(-GAUSS_NODE, -GAUSS_NODE, GAUSS_NODE);
	  _nodes(4).assign(GAUSS_NODE, GAUSS_NODE, -GAUSS_NODE);
	  _nodes(5).assign(-GAUSS_NODE, GAUSS_NODE, -GAUSS_NODE);
  	_nodes(6).assign(GAUSS_NODE, -GAUSS_NODE, -GAUSS_NODE);
	  _nodes(7).assign(-GAUSS_NODE, -GAUSS_NODE, -GAUSS_NODE);

    _weights.assign(1.0);
  }


  template<>
  TetrahedraIntegrator<1>::TetrahedraIntegrator() : ElementIntegrator<3, 1>()
  {
    _nodes(0).assign(0.25, 0.25, 0.25); // barycenter node

    _weights.assign(1.0/6.0);
  }




  template<>
  IntervalIntegrator<2>::IntervalIntegrator() : ElementIntegrator<1, 2>()
  {
    static const float_t GAUSS_NODE = 0.5773502694;

    _nodes(0).assign(-GAUSS_NODE);
    _nodes(1).assign(GAUSS_NODE);

    _weights.assign(1.0);
  }

  template<>
  IntervalIntegrator<1>::IntervalIntegrator() : ElementIntegrator<1, 1>()
  {
    _nodes(0).assign(0.0);

    _weights.assign(2.0);
  }
}

