#include <minimize/Lbfgs.hpp>

#include <core/utilities.hpp>
#include <core/MessageInterface.hpp>

#include <boost/lexical_cast.hpp>


namespace imaging
{
  Lbfgs::Lbfgs(DifferentiableEnergyInterface & energy,
    float_t epsilon,
    float_t line_search_f_tolerance,
    float_t line_search_g_tolerance,
    size_t n_correction_steps,
    size_t n_max_line_search_steps) : 
    _energy(energy),
    _terminated(false),
    _n_steps(0)
  {
    lbfgs_parameter_init(& _lbfgs_parameters);
    
    _lbfgs_parameters.epsilon = epsilon;
    _lbfgs_parameters.ftol = line_search_f_tolerance;
    _lbfgs_parameters.gtol = line_search_g_tolerance;
    _lbfgs_parameters.m = n_correction_steps;
    _lbfgs_parameters.max_linesearch = n_max_line_search_steps;
  }
    
  bool Lbfgs::minimize(size_t n_max_steps,
    size_t & n_steps)
  {
    if(_terminated)
    {
      n_steps = 0;
      return true;
    }
      
    if(n_max_steps == 0)
    {
      n_steps = 0;
      return true;
    }
      
    _lbfgs_parameters.max_iterations = n_max_steps;
    
    size_t n = 0;
    
    // Todo: multiples of 4 necessary???
    while(n < _energy.dimension())
      n += 4;
     
    n = _energy.dimension();
    
    int ret;
    lbfgsfloatval_t x[n], fx;
    
    for(size_t i = 0; i < _energy.dimension(); ++i)
      x[i] = _energy.current_argument()(i);
    
    ret = lbfgs(n, x, &fx, evaluate, progress, this, & _lbfgs_parameters);

    /* Report the result. */
    MessageInterface::out("L-BFGS optimization terminated with status code " + boost::lexical_cast<std::string>(ret) +
      ", fx = " +  boost::lexical_cast<std::string>(float_t(fx)), MessageInterface::DEBUG_ONLY);
    
    n_steps = _n_steps;
    return _terminated;
  }
  
  lbfgsfloatval_t Lbfgs::evaluate(void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step)
  {
    Lbfgs * object = static_cast<Lbfgs*>(instance);
    if(n < object->_energy.dimension())
      throw Exception("Exception: n too small in Lbfgs::evaluate().");
      
    for(size_t i = 0; i < object->_energy.dimension(); ++i)
      object->_energy.current_argument()(i) = float_t(x[i]);
    
    object->_energy.set_argument_with_gradient();
    
    for(size_t i = 0; i < object->_energy.dimension(); ++i)
      g[i] = lbfgsfloatval_t(object->_energy.current_gradient()(i));
      
    for(size_t i = object->_energy.dimension(); i < n; ++i)
      g[i] = 0.0;
    
    return object->_energy.current_energy();
  }
  
  
  int Lbfgs::progress(void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls)
  {
    
    Lbfgs * object = static_cast<Lbfgs*>(instance);
    
    if(gnorm < object->_lbfgs_parameters.epsilon * max(1.0, xnorm))
      object->_terminated = true;
      
    ++object->_n_steps;
    MessageInterface::out("NL-BFGS: Step " + boost::lexical_cast<std::string>(object->_n_steps), MessageInterface::DEBUG_ONLY);
    MessageInterface::out("Functional value: " + boost::lexical_cast<std::string>(fx),
            MessageInterface::DEBUG_ONLY, +1);
    MessageInterface::out("Norm of gradient: " + boost::lexical_cast<std::string>(gnorm),
            MessageInterface::DEBUG_ONLY, +1);
    
    return 0;
  }
  
}

