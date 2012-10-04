// This file is part of the imaging2 class library.
//
// University of Innsbruck, Infmath Imaging, 2009.
// http://infmath.uibk.ac.at
//
// All rights reserved.


namespace imaging {}
namespace img = imaging;

#include <complex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <core/ublas/fixed_vector.hpp>
#include <core/ublas/fixed_matrix.hpp>

#include <core/Exception.hpp>

#include <core/float_types.hpp>

/** \ingroup core
    <tt>\#include <core/imaging2.hpp></tt>
    \brief Main library namespace (abbreviated as \em img).
    
    This is the main namespace of \em imaging2. All objects in the class library belong to this namespace with the only exception of the classes in tla. Instead of \em imaging also the short version \em img can be used.
*/
namespace imaging
{
  namespace ublas = boost::numeric::ublas;
  
  /** \ingroup core
      <tt>\#include <core/float_types.hpp></tt>
      
      Defines \em imaging::size_t to be \em std::size_t.
  */
  typedef std::size_t size_t;
  
  /** \ingroup core
      <tt>\#include <core/float_types.hpp></tt>
      
      Defines \em imaging::complex_t to be <em>std::complex<float_t></em>.
  */
  typedef std::complex<float_t> complex_t;
}
