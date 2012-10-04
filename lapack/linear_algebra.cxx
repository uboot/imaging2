#include <lapack/linear_algebra.hpp>

#include <core/utilities.hpp>
#include <lapack/imaging_lapack.hpp>


namespace imaging
{
  void eigensystem(const ublas::matrix<float_t> & in, ublas::matrix<float_t> & eigenvectors, ublas::vector<float_t> & eigenvalues)
  {
    if (in.size1() == 0 || in.size2() == 0)
      throw Exception("Exception: empty matrix in eigensystem()."); 
      
    if (in.size1() != in.size2())
      throw Exception("Exception: non-square matrix in eigensystem()."); 
      
    std::size_t size = in.size1();
    
    char uplo = 'u';
    int n = int(size);
    float_t* e = new float_t[n-1];
    float_t* tau = new float_t[n-1];
    int lwork = n * 4;
    float_t* work = new float_t[lwork];
    int info;

    eigenvectors.resize(size, size);
    eigenvalues.resize(size);
    
    ublas::matrix<float_t> temp_eigenvectors(in);
    ublas::vector<float_t> temp_eigenvalues(eigenvalues.size());

    float_t *temp_eigenvectors_ptr, *temp_eigenvalues_ptr;

    temp_eigenvectors_ptr = &temp_eigenvectors(0,0);
    temp_eigenvalues_ptr = &temp_eigenvalues(0);

    IMAGING_DSYTRD(uplo, n, temp_eigenvectors_ptr, n, temp_eigenvalues_ptr, e,
            tau, work, lwork, &info);
            
    if(info != 0)
      throw Exception("Exception: Error in (lapack) dsytrd_().");

    IMAGING_DORGTR(uplo, n, temp_eigenvectors_ptr, n,
            tau, work, lwork, &info);
            
    if(info != 0)
      throw Exception("Exception: Error in (lapack) dorgtr_().");

    char compz = 'V';
    
    float_t* work_2 = new float_t[max(1,2*n-2)];

    IMAGING_DSTEQR(compz, n, temp_eigenvalues_ptr, e, temp_eigenvectors_ptr,
            n, work_2, &info);
    
    for(ublas::matrix<float_t>::iterator1 iter1 = eigenvectors.begin1(); iter1 != eigenvectors.end1(); ++iter1)
      for(ublas::matrix<float_t>::iterator2 iter2 = iter1.begin(); iter2!= iter1.end(); ++iter2)
        *iter2 = temp_eigenvectors(eigenvectors.size2() - iter2.index2() - 1, iter2.index1());
        
    for(ublas::vector<float_t>::iterator iter = eigenvalues.begin(); iter != eigenvalues.end(); ++iter)
      *iter = temp_eigenvalues[eigenvalues.size() - iter.index() - 1];
            
    delete [] e, tau, work, work_2;
  }
  
  void symmetric_square_root(const ublas::matrix<float_t> & in, ublas::matrix<float_t> & root)
  { 
    if (in.size1() == 0 || in.size2() == 0)
      throw Exception("Exception: empty matrix in symmetric_square_root()."); 
      
    if (in.size1() != in.size2())
      throw Exception("Exception: non-square matrix in symmetric_square_root()."); 
    
    const std::size_t n = in.size1();
    
    ublas::matrix<float_t> eigenvectors(n, n);
    ublas::vector<float_t> eigenvalues(n);
    
    eigensystem(in, eigenvectors, eigenvalues);
    
    root.resize(n, n);
    root = ublas::scalar_matrix<float_t>(n, n, 0.0);
    
    for(ublas::matrix<float_t>::iterator1 iter1 = root.begin1(); iter1 != root.end1(); ++iter1)
      for(ublas::matrix<float_t>::iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        *iter2 = eigenvectors(iter2.index1(), iter2.index2()) * 
                 sqrt(fabs(eigenvalues(iter2.index2())));
  }
  
  void square_root(const ublas::matrix<float_t> & in, ublas::matrix<float_t> & root)
  { 
    if (in.size1() == 0 || in.size2() == 0)
      throw Exception("Exception: empty matrix in square_root()."); 
      
    if (in.size1() != in.size2())
      throw Exception("Exception: non-square matrix in square_root()."); 
    
    const std::size_t n = in.size1();
    
    ublas::matrix<float_t> eigenvectors(n, n);
    ublas::vector<float_t> eigenvalues(n);
    
    eigensystem(in, eigenvectors, eigenvalues);
    
    root.resize(n, n);
    
    for(ublas::matrix<float_t>::iterator1 iter1 = root.begin1(); iter1 != root.end1(); ++iter1)
      for(ublas::matrix<float_t>::iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        *iter2 = eigenvectors(n - iter2.index1() - 1, iter2.index2()) * 
                 sqrt(fabs(eigenvalues(iter2.index2())));
                 
    root = prod(root, trans(eigenvectors));
  }
  
  void inverse_square_root(const ublas::matrix<float_t> & in, ublas::matrix<float_t> & root)
  { 
    if (in.size1() == 0 || in.size2() == 0)
      throw Exception("Exception: empty matrix in inverse_square_root()."); 
      
    if (in.size1() != in.size2())
      throw Exception("Exception: non-square matrix in inverse_square_root()."); 
    
    const std::size_t n = in.size1();
    
    ublas::matrix<float_t> eigenvectors(n, n);
    ublas::vector<float_t> eigenvalues(n);
    
    eigensystem(in, eigenvectors, eigenvalues);
    
    root.resize(n, n);
    
    for(ublas::matrix<float_t>::iterator1 iter1 = root.begin1(); iter1 != root.end1(); ++iter1)
      for(ublas::matrix<float_t>::iterator2 iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
        *iter2 = eigenvectors(iter2.index1(), iter2.index2()) / 
                 sqrt(fabs(eigenvalues(iter2.index2())));
  }
  
  void inverse(const ublas::matrix<float_t> & in, ublas::matrix<float_t> & out)
  {
    if (in.size1() == 0 || in.size2() == 0)
      throw Exception("Exception: empty matrix in inverse()."); 
      
    if (in.size1() != in.size2())
      throw Exception("Exception: non-square matrix in inverse()."); 
      
    int n = in.size1();
    
    out = in;
    
    int * ipiv = new int[n];
    int info;
    int lwork = n * 4;
    float_t* work = new float_t[lwork];
    float_t* out_ptr = &out(0, 0);

    IMAGING_DGETRF(n, n, out_ptr, n, ipiv, &info); // C style interface!

    if(info != 0)
      throw Exception("Exception: Error in (lapack) dgetrf_().");

    IMAGING_DGETRI(n, out_ptr, n, ipiv, work, lwork, &info);

    if(info != 0)
      throw Exception("Exception: Error in (lapack) dgetri_().");
    
    delete [] ipiv, work;
  }

}

