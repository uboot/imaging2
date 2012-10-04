#include <statistic/utilities.hpp>
#include <statistic/LinearPca.hpp>

#include <core/cio.hpp>

using namespace imaging;

int main ( int argc, char **argv )
{
  try
  {
    ublas::matrix<imaging::float_t> data(5, 3);
    
    data(0,0) = 1.0;
    data(0,1) = 0.0;
    data(0,2) = 0.0;
    
    data(1,0) = 2.1;
    data(1,1) = 0.9;
    data(1,2) = 0.1;
    
    data(2,0) = 3.3;
    data(2,1) = 2.1;
    data(2,2) = 0.2;
    
    data(3,0) = 0.2;
    data(3,1) = -1.1;
    data(3,2) = 0.35;
    
    data(4,0) = -1.3;
    data(4,1) = -2.5;
    data(4,2) = 0.45;
    
    
    std::cout << "Statistic utilities\n" << std::endl;
    
    ublas::vector<imaging::float_t> data_mean;
    ublas::vector<imaging::float_t> data_variance;
    
    imaging::mean(data, data_mean);
    imaging::var(data, data_variance);
    
    std::cout << "Mean: " << data_mean << std::endl;
    std::cout << "Variance: " << data_variance << std::endl;
    
    std::cout << "\n2D-PCA of a dataset in 3D\n" << std::endl;
    
    LinearPca pca(data, 2);
    
    std::cout << "Data:\n " << data << std::endl;
    
    std::cout << "Mean:\n " << pca.mean() << std::endl;
    
    std::cout << "Standard deviations:\n ";
    std::cout << pca.standard_deviations() << std::endl;
    
    ublas::vector<imaging::float_t> coefficients(pca.dimension());
    ublas::vector<imaging::float_t> vector(3);
    std::cout << "Principal components:\n ";
    
    coefficients(0) = 1.0;
    coefficients(1) = 0.0;
    pca.compute_vector(coefficients, vector);
    std::cout << vector - pca.mean();
    
    coefficients(0) = 0.0;
    coefficients(1) = 1.0;
    pca.compute_vector(coefficients, vector);
    std::cout << vector - pca.mean();

    std::cout << std::endl;
    
    vector(0) = 0.0;
    vector(1) = 0.0;
    vector(2) = 0.0;
    std::cout << "PCA coefficients of origin:\n ";
    pca.compute_coefficients(vector, coefficients);
    std::cout << coefficients << std::endl;
    
    std::cout << "Reconstruction of origin:\n ";
    pca.compute_vector(coefficients, vector);
    std::cout << vector << std::endl;
    
    std::cout << "Mahalanobis norm of origin:\n ";
    std::cout << pca.norm(vector) << std::endl;

  }

  catch ( Exception &exception )
  {
    std::cerr << exception.error_msg() << std::endl;

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
