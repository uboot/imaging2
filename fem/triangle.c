#define VOID int
#define REAL double
#define NULL 0x0

#include <string.h>

#include <external/triangle/triangle.h>

void triangle_triangulate(unsigned int n_in_points,
                   const double* in_points,
                   double max_triangle_area,
                   unsigned int *n_out_points,
                   double **out_points,
                   unsigned int *n_out_elements,
                   unsigned int **out_elements,
                   unsigned int *n_out_boundary_elements,
                   unsigned int **out_boundary_elements)
{
  unsigned int i, j;
  
  struct triangulateio in, out;

  /* Define input points. */

  in.numberofpoints = n_in_points;
  in.numberofpointattributes = 0;
  in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

  for(i = 0; i < 2 * n_in_points; ++i)
    in.pointlist[i] = in_points[i];
    
  in.numberofsegments = n_in_points;
  in.segmentlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(int));
  in.segmentmarkerlist = (REAL *) malloc(in.numberofpoints * sizeof(int));
  
  for(i = 0; i < n_in_points; ++i)
  {
    in.segmentlist[2*i] = i;
    in.segmentlist[2*i+1] = (i+1) % n_in_points;
    in.segmentmarkerlist[i] = 1;
  }

  in.pointattributelist = (REAL *) NULL;
  in.pointmarkerlist = (int *) NULL;

  in.numberofholes = 0;
  in.numberofregions = 0;
  in.regionlist = (REAL *) NULL;
  
  /* Make necessary initializations so that Triangle can return a */
  /*   triangulation in `out' and a voronoi diagram in `vorout'.  */

  out.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  /* Not needed if -N switch used or number of point attributes is zero: */
  out.pointattributelist = (REAL *) NULL;
  out.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  out.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  /* Not needed if -E switch used or number of triangle attributes is zero: */
  out.triangleattributelist = (REAL *) NULL;
  out.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  /* Needed only if segments are output (-p or -c) and -P not used: */
  out.segmentlist = (int *) NULL;
  /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  out.segmentmarkerlist = (int *) NULL;
  out.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  out.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */
  
  /* Triangulate the points.  Switches are chosen to read and write a  */
  /*   PSLG (p), preserve the convex hull (c), number everything from  */
  /*   zero (z), assign a regional attribute to each element (A), and  */
  /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  /*   neighbor list (n).                                              */

  char max_triangle_area_string[10];
  char attribute_list[17] = "QpqYeza";
  
  snprintf(max_triangle_area_string, 10, "%f", max_triangle_area);
  strcat(attribute_list, max_triangle_area_string);
  
  triangulate(attribute_list, &in, &out, NULL);
  
  *n_out_points = out.numberofpoints;
  *out_points = (double*)malloc(out.numberofpoints * 2 * sizeof(double));
  
  for(i = 0; i < 2 * out.numberofpoints; ++i)
    (*out_points)[i] = out.pointlist[i];
    
  *n_out_elements = out.numberoftriangles;
  *out_elements = (unsigned int*)malloc(out.numberoftriangles * 3 * sizeof(unsigned int));
  
  for(i = 0; i < 3 * out.numberoftriangles; ++i)
    (*out_elements)[i] = (unsigned int)out.trianglelist[i];
    
  *n_out_boundary_elements = 0;  
  for(i = 0; i < out.numberofedges; ++i)
    if(out.edgemarkerlist[i] == 1)
      (*n_out_boundary_elements)++;
      
  *out_boundary_elements = (unsigned int*)malloc(*n_out_boundary_elements * 2 * sizeof(unsigned int));
  
  for(i = 0, j = 0; i < out.numberofedges; ++i)
    if(out.edgemarkerlist[i] == 1)
    {
      (*out_boundary_elements)[2 * j] = (unsigned int)out.edgelist[2 * i];
      (*out_boundary_elements)[2 * j + 1] = (unsigned int)out.edgelist[2 * i + 1];
      ++j;
    }
  
  
  
  free(in.pointlist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  /* free(in.pointattributelist);
  free(in.pointmarkerlist);
  free(in.regionlist); */
  free(out.pointlist);
  // free(out.pointattributelist);
  free(out.pointmarkerlist);
  free(out.trianglelist);
  free(out.edgelist);
  free(out.edgemarkerlist);
  // free(out.triangleattributelist);
  /* free(out.trianglearealist);
  free(out.neighborlist); */
  free(out.segmentlist);
  free(out.segmentmarkerlist);
  
  
}

