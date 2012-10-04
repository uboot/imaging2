#include <minimize/SteepestDescent.hpp>

#include <core/MessageInterface.hpp>

#include <boost/lexical_cast.hpp>

namespace imaging
{
  SteepestDescent::SteepestDescent(DifferentiableEnergyInterface & energy, float_t min_gradient_norm, float_t step_size) :
    _energy(energy), _min_gradient_norm(min_gradient_norm), _step_size(step_size), _terminated(false) {}
    
  bool SteepestDescent::minimize(size_t n_max_steps, size_t & n_steps)
  {
    if(_terminated)
    {
      n_steps = 0;
      return true;
    }
      
    for(n_steps = 0; n_steps < n_max_steps; ++n_steps)
    {
      MessageInterface::out("SteepestDescent: Step " + boost::lexical_cast<std::string>(n_steps + 1), MessageInterface::DEBUG_ONLY);

      if(_energy.current_argument().size() != _energy.current_gradient().size())
        throw Exception("Exception: Sizes of current_argument() and current_gradient() do not agree in SteepestDescent::minimize().");
        
      _energy.current_argument() -= _step_size * _energy.current_gradient();
      _energy.set_argument_with_gradient();
      
      float_t gradient_abs = norm_2(_energy.current_gradient());

      MessageInterface::out("Current functional value: " + boost::lexical_cast<std::string>(_energy.current_energy()),
               MessageInterface::DEBUG_ONLY, +1);
      MessageInterface::out("Current gradient (abs): " + boost::lexical_cast<std::string>(gradient_abs),
               MessageInterface::DEBUG_ONLY, +1);

      if(gradient_abs < _min_gradient_norm)
      {
        _terminated = true;
        return true;
      }
    }
    
    return false;
  }

}

