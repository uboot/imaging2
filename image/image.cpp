#include <image/Image.hpp>
#include <image/Color.hpp>

#include <graphics/OpenGlViewer.hpp>
#include <image/gio.hpp>



int main(int argc, char ** argv)
{
  using namespace img;
  
  typedef img::float_t float_t;
  typedef img::size_t size_t;
  
  try
  {
    OpenGlViewer::out.init(argc, argv, ublas::fixed_vector<float_t, 2> (0.0, 0.0),
    ublas::fixed_vector<float_t, 2>(200.0, 100.0), 800, 400);

    pthread_t display_thread;
    pthread_create(&display_thread, NULL, &OpenGlViewer::out.display, 0x0);

    ColorImage2d image(ublas::fixed_vector<size_t, 2>(200, 100));
    
    size_t i = 0;
    for(boost::multi_array<Color, 2>::iterator iter1 = image.begin(); iter1 != image.end(); ++iter1)
      for(boost::multi_array<Color, 1>::iterator iter2 = iter1->begin(); iter2 != iter1->end(); ++iter2, ++i)
        *iter2 = Color(i, i * i, i * i * i);
    
    OpenGlViewer::out << image;
    
    int * ret;
    pthread_join(display_thread, (void**)&ret);

  }
  catch(Exception & e)
  {
    std::cerr << e.error_msg() << std::endl;
  }
  return 0;
}
