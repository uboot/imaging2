#include <graphics/OpenGlViewer.hpp>
#include <core/Cmessage.hpp>
#include <shape/mrep/MrepModel2d.hpp>
#include <shape/mrep/PolygonModel2d.hpp>
#include <shape/mrep/gio.hpp>
#include <shape/mrep/xmlio.hpp>


using namespace imaging;


int main ( int argc, char **argv )
{
  try
  {
    typedef OpenGlViewer gr;
    
    img::size_t verbose_level = Cmessage::ALL_MESSAGES;
    img::size_t main_argument_base_index = 1;
    
    if(argc == 1)
      throw Exception("No command line parameters provided!");
      
    if(argc == 3)
    {
      verbose_level = boost::lexical_cast<img::size_t>(argv[1]);
      main_argument_base_index = 2;
    }
    
    Cmessage::out.set_verbose_level(verbose_level);
    gr::out.init (argc, argv, ublas::fixed_vector<imaging::float_t, 2>(0.0, 0.0), ublas::fixed_vector<imaging::float_t, 2> ( 600.0, 300.0 ), 512, 256);
    
    pthread_t display_thread;
    if(pthread_create(&display_thread, NULL, &img::OpenGlViewer::out.display, 0x0) != 0 )
      throw img::Exception("Exception: Could not create OpenGL thread.");
      
    XmlReader xml_in(argv[main_argument_base_index], "imaging");
    std::vector<MrepModel2d> mrep_models;
    std::vector<PolygonModel2d> polygon_models;
    xml_in >> mrep_models;
    xml_in >> polygon_models;
    
    for(std::size_t i = 0; i < mrep_models.size(); ++i)
      gr::out << mrep_models[i];
    
    for(std::size_t i = 0; i < polygon_models.size(); ++i)
      gr::out << polygon_models[i];

    gr::out << gr::flush;
    
    int *ret;
    pthread_join(display_thread, (void**)&ret);
  }

  catch ( img::Exception &exception )
  {
    std::cerr << exception.error_msg() << std::endl;

    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

