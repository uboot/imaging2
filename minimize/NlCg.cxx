#include <minimize/NlCg.hpp>

#include <core/MessageInterface.hpp>

#include <boost/lexical_cast.hpp>

extern "C" 
{
  void cgfam_( int * N, double * X, double * F, double * G, double * D, double * GOLD, int * IPRINT, double * EPS, double * W, int * IFLAG, int * IREST, int * METHOD, int * FINISH );
}
  
  

namespace imaging
{
  NlCg::NlCg(DifferentiableEnergyInterface & energy, float_t convergence_constant, size_t method) : 
    _energy(energy),
    _convergence_constant(convergence_constant),
    _d(energy.dimension()),
    _g_old(energy.dimension()),
    _w(energy.dimension()),
    _iflag(0),
    _finish(false)
  {
    switch(method)
    {
    case FLETCHER_REEVES:
      _method = 1;
      break;
    case POLAK_RIBIERE:
      _method = 2;
      break;
    case POSITIVE_POLAK_RIBIERE:
      _method = 3;
      break;
    default: 
      throw Exception("Exception: Invalid method specified in NlCg::minimize().");
    }
  }
    
  bool NlCg::minimize(size_t n_max_steps,
    size_t & n_steps)
  {
    if(_finish)
    {
      n_steps = 0;
      return true;
    }
      
    int iprint[2];
    int irest = 0;
    int n = _energy.dimension();
    
    iprint[0] = -1;
    iprint[1] = 0;
    
    float_t f, gnorm;
    ublas::vector<double> g(_energy.dimension());
    ublas::vector<double> x(_energy.dimension());
    
    x = _energy.current_argument();
    
    
    for(n_steps = 0; n_steps < n_max_steps; ++n_steps)
    {
      MessageInterface::out("Nonlinear Conjugated Gradients: Step " + boost::lexical_cast<std::string>(n_steps + 1), MessageInterface::DEBUG_ONLY);

      _energy.set_argument_with_gradient();
      f = _energy.current_energy();
      g = _energy.current_gradient();
      
      gnorm = norm_2(g);
      
      MessageInterface::out("Functional value: " + boost::lexical_cast<std::string>(f),
               MessageInterface::DEBUG_ONLY, +1);
      MessageInterface::out("Norm of gradient: " + boost::lexical_cast<std::string>(gnorm),
               MessageInterface::DEBUG_ONLY, +1);
      
      for( ; ; )
      {
        cgfam_(& n, & x(0), & f, & g(0), & _d(0), & _g_old(0), iprint, & _convergence_constant, & _w(0), & _iflag, & irest, & _method, & _finish);
        
        if(_iflag == 1 || _iflag <= 0 || _iflag > 2)
          break;
        
        if(gnorm <= _convergence_constant)
          _finish = true;
      }
      
      if(_iflag <= 0 || _iflag > 2)
        break;
        
      _energy.current_argument() = x;
      
    }
      
    _energy.set_argument_with_gradient();
    
    MessageInterface::out("NL-CG optimization terminated with status code " + boost::lexical_cast<std::string>(_iflag) + ".", MessageInterface::DEBUG_ONLY);
    
    return _finish;
  }
  
}

