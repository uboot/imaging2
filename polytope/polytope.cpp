
#include <core/Cmessage.hpp>
#include <graphics/OpenGlViewer.hpp>
#include <polytope/SimplePolygon.hpp>
#include <polytope/gio.hpp>
#include <polytope/xmlio.hpp>


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
    gr::out.init (argc, argv, ublas::fixed_vector<imaging::float_t, 2>(-1.0, -1.0), ublas::fixed_vector<imaging::float_t, 2> ( 3.0, 3.0 ), 400, 400);
    
    pthread_t display_thread;
    if(pthread_create(&display_thread, NULL, &img::OpenGlViewer::out.display, 0x0) != 0 )
      throw img::Exception("Exception: Could not create OpenGL thread.");
      
    XmlReader xml_in(argv[main_argument_base_index], "imaging");
    std::vector<Polygon> polygons;
    Polygon result;
    xml_in >> polygons;
    
    for(std::size_t i = 0; i < polygons.size(); ++i)
    {
      gr::out << polygons[i];
      polygon_union(result, polygons[i], result);
    }
    
    gr::out << gr::offset_z_layer(+1) << gr::set_color(Color::RED);
    gr::out << result;
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

