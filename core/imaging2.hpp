/* 
*  Copyright 2009 University of Innsbruck, Infmath Imaging
*
*  This file is part of imaging2.
*
*  Imaging2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  Imaging2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with stromx-studio.  If not, see <http://www.gnu.org/licenses/>.
*/

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
