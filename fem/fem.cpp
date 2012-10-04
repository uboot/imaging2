//Including all neccessary header files

#include <core/imaging2.hpp>

#include <image/Image.hpp>
#include <image/Color.hpp>
#include <image/CastAccessor.hpp>
#include <fem/Grid.hpp>
#include <fem/fem_2d_square_types.hpp>
#include <fem/equation/DiffusionStep.hpp>
#include <fem/equation/McmStep.hpp>
#include <fem/equation/TvFlowStep.hpp>
#include <fem/equation/SimpleEquationAdaptor.hpp>
#include <fem/Image2Grid.hpp>
#include <fem/SimpleAssembler.hpp>
#include <solver/CgSolver.hpp>

#include <graphics/OpenGlViewer.hpp>
#include <image/gio.hpp>
#include <fem/gio.hpp>


using namespace img;


void filter(std::string file_name, img::float_t diffusion_step_size, img::float_t mcm_step_size, img::float_t tvf_step_size, img::size_t n_steps)
{
    typedef img::float_t float_t;
    Image<2, GrayValue> input_image(file_name);
    CastAccessor< GrayImage2d, float_t> input_accessor(input_image);

    //Preparations for finite element method;
    //Construction of the grid and transformation of the image to a vector
    typedef fem_2d_square_types fem_types;

    Grid<fem_types> grid;
    ublas::compressed_matrix<float_t> stiffness_matrix;
    boost::shared_ptr< ublas::vector<float_t> > input_1(new ublas::vector<float_t>);
    boost::shared_ptr< ublas::vector<float_t> > input_2(new ublas::vector<float_t>);
    boost::shared_ptr< ublas::vector<float_t> > input_3(new ublas::vector<float_t>);

    Image2Grid<fem_types> image2grid(input_image.size());
    image2grid.construct_grid(grid);
    image2grid.stiffness_matrix_prototype(stiffness_matrix);
    image2grid.image2vector(input_accessor, *input_1);
    image2grid.image2vector(input_accessor, *input_2);
    image2grid.image2vector(input_accessor, *input_3);

    boost::shared_ptr< ublas::vector<float_t> > diffusion(new ublas::vector<float_t>);
    *diffusion = ublas::scalar_vector<float_t>(input_3->size(), 1.0);


    //Setting up different types of equations and initialization of parameters
    //mean curvature motion
    McmStep<fem_types> mcm_equation;
    mcm_equation.set_step_size(mcm_step_size);
    mcm_equation.set_input(input_1);
    mcm_equation.set_epsilon(1.0e-4);
    //total variation flow 
    TvFlowStep<fem_types> tvf_equation;
    tvf_equation.set_step_size(tvf_step_size);
    tvf_equation.set_input(input_2);
    tvf_equation.set_epsilon(1.0e-4);
    //diffusion equation 
    DiffusionStep<fem_types> diffusion_equation;
    diffusion_equation.set_step_size(diffusion_step_size);
    diffusion_equation.set_diffusion(diffusion);
    diffusion_equation.set_input(input_3);


    //Definition of data types;
    //assembler -> ??
    SimpleAssembler assembler;
    //force vector -> right hand side of the equation
    ublas::vector<float_t> force_vector;
    //solver -> type of solver 
    CgSolver solver;


    //Creating output accessors 
    //mean curvature motion    
    GrayImage2d mcm_output_image(input_image.size());
    CastAccessor<GrayImage2d, float_t> mcm_output_accessor(mcm_output_image);
    //total variation flow 
    GrayImage2d tvf_output_image(input_image.size());
    CastAccessor<GrayImage2d, float_t> tvf_output_accessor(tvf_output_image);
    //total variation flow 
    GrayImage2d diffusion_output_image(input_image.size());
    CastAccessor<GrayImage2d, float_t> diffusion_output_accessor(diffusion_output_image);


    for (std::size_t i=1; i <= n_steps; ++i)
    {
      std::cout << "Iteration step: " << i << std::endl;

      /*Assembling of stiffness matrix, force vector of the corresponding 
        equation and solving the system of equations*/
      //mean curvature motion
      assembler.assemble(mcm_equation, grid, stiffness_matrix, force_vector);
      solver.solve(stiffness_matrix, force_vector, *input_1);

      //total variation flow 
      assembler.assemble(tvf_equation, grid, stiffness_matrix, force_vector);
      solver.solve(stiffness_matrix, force_vector, *input_2);

      //diffusion equation
      assembler.assemble(diffusion_equation, grid, stiffness_matrix, force_vector);
      solver.solve(stiffness_matrix, force_vector, *input_3);


      //Converting the vector of the solution back to an image 
      //mean curvature motion 
      image2grid.vector2image(*input_1, mcm_output_accessor);
      //total variation flow 
      image2grid.vector2image(*input_2, tvf_output_accessor);
      //diffusion equation 
      image2grid.vector2image(*input_3, diffusion_output_accessor);

      //Graphic output
      OpenGlViewer::out << GraphicsInterface::reset_group;
      GraphicsInterface::StreamStatus old_status = OpenGlViewer::out.get_stream_status();
      OpenGlViewer::out << diffusion_output_image;
      OpenGlViewer::out << GraphicsInterface::translate(ublas::fixed_vector< float_t, 2 >(0, float_t(input_image.size()[1])));
      OpenGlViewer::out << mcm_output_image;
      OpenGlViewer::out << GraphicsInterface::translate(ublas::fixed_vector< float_t, 2 >(0, float_t(mcm_output_image.size()[1])));
      OpenGlViewer::out << tvf_output_image;
      OpenGlViewer::out << GraphicsInterface::flush;
      OpenGlViewer::out << old_status;
    }
    
//     OpenGlViewer::out.write_image("test.eps");

}


int main(int argc, char ** argv)
{
  typedef img::float_t float_t;

  try
  {
    OpenGlViewer::out.init(argc, argv, ublas::fixed_vector<float_t, 2> (0.0, 0.0),
    ublas::fixed_vector<float_t, 2>(100.0, 150.0), 400, 600);

    pthread_t display_thread;
    pthread_create(&display_thread, NULL, &OpenGlViewer::out.display, 0x0);

    OpenGlViewer::out << GraphicsInterface::set_group(0);
    filter("lenna.png", 0.5, 1.0, 0.1, 5);
    OpenGlViewer::out << GraphicsInterface::translate(ublas::fixed_vector< float_t, 2 >(50,0));
    OpenGlViewer::out << GraphicsInterface::set_group(1);
    filter("square.png", 0.1, 5.0, 0.05, 50);

    std::cout << "Successful!" << std::endl;

    int * ret;
    pthread_join(display_thread, (void**)&ret);

  }
  catch(Exception & e)
  {
    std::cerr << e.error_msg() << std::endl;
  }
  return 0;
}
