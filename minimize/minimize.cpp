#include <core/utilities.hpp>
#include <core/Cmessage.hpp>
#include <minimize/SteepestDescent.hpp>
#include <minimize/Lbfgs.hpp>
#include <minimize/CovarianceMatrixAdaptation.hpp>
#include <minimize/NlCg.hpp>
#include <core/DifferentiableFunctionalInterface.hpp>
#include <minimize/DifferentiableFunctionalAdaptor.hpp>


using namespace imaging;

class QuadraticFunctional : public img::DifferentiableFunctionalInterface
{
public:
  img::float_t operator()(const ublas::vector<img::float_t> & x)
  {
    return 3.0 * x(0) * x(0) * x(0) * x(0) + x(1) * x(1) * x(1) * x(1);
  }
  
  img::float_t operator()(const ublas::vector<img::float_t> & x, ublas::vector<img::float_t> & gradient)
  {
    gradient.resize(2);
    gradient(0) = 12.0 * x(0) * x(0) * x(0);
    gradient(1) = 4.0 * x(1) * x(1) * x(1);
    
    return operator()(x);
  }
  
  img::size_t dimension() const
  {
    return 2;
  }
};

void reset_energy(DifferentiableEnergyInterface & energy)
{
  energy.current_argument()(0) = 1.0;
  energy.current_argument()(1) = 2.0;
  
  energy.set_argument_with_gradient();
}


int main ( int argc, char **argv )
{
  try
  {
    Cmessage::out.set_verbose_level ( Cmessage::ALL_MESSAGES );
    
    QuadraticFunctional my_functional;
    DifferentiableFunctionalAdaptor<QuadraticFunctional> functional_adaptor(my_functional);
    
    img::float_t epsilon = 1e-2;
    img::size_t n_max_steps = 100;
    
    if(argc > 1)
      epsilon = boost::lexical_cast<img::float_t>(argv[1]);
    
    if(argc > 2)
      n_max_steps = boost::lexical_cast<img::size_t>(argv[2]);
      
    SteepestDescent steepest_descent(functional_adaptor, epsilon, 0.1);
    Lbfgs lbfgs(functional_adaptor, epsilon);
    CovarianceMatrixAdaptation cma(functional_adaptor, 1.0, epsilon);
    NlCg nlcg(functional_adaptor, epsilon);
    
    std::vector<std::string> labels(4);
    std::vector<MinimizerInterface *> minimizers(4);
    
    labels[0] = "Steepest Descent (gradient based)";
    minimizers[0] = & steepest_descent;
    
    labels[1] = "LBFGS (quasi-Newton method)";
    minimizers[1] = & lbfgs;
    
    labels[2] = "CMA (evolutionary algorithm)";
    minimizers[2] = & cma;
    
    labels[3] = "NL CG (nonlinear conjugated gradient, Fletcher-Reeves)";
    minimizers[3] = & nlcg;
    
    for(std::size_t i = 0; i < 4; ++i)
    {
      bool terminated;
      img::size_t n_actual_steps;
      reset_energy(functional_adaptor);
      
      MessageInterface::out(labels[i], MessageInterface::IMPORTANT);
      terminated = minimizers[i]->minimize(n_max_steps, n_actual_steps);
      MessageInterface::out("terminated: " + boost::lexical_cast<std::string>(terminated), MessageInterface::IMPORTANT, +1);
      MessageInterface::out("n_actual_steps: " + boost::lexical_cast<std::string>(n_actual_steps), MessageInterface::IMPORTANT, +1);
      MessageInterface::out("energy: " + boost::lexical_cast<std::string>(functional_adaptor.current_energy()), MessageInterface::IMPORTANT, +1);
    }
  }

  catch ( Exception &exception )
  {
    std::cerr << exception.error_msg() << std::endl;

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
