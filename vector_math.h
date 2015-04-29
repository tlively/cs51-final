/********************************************************************
 * the header file for some handy vector math
 ********************************************************************/
#ifndef VECTOR_INCLUDED
#define VECTOR_INCLUDED

/* A generic struct for handing two dimensional vectors. Used in many
 * situations with different meanings. */
typedef struct po_vector {
  float x;
  float y;
} po_vector;

/* forms a line vector from two coordinate points */
po_vector vect_from_points(po_vector origin, po_vector destination);

/* gets the right hand normal vector */
po_vector vect_normal(po_vector vect);

/* calculate the dot product of two vectors */
float vect_dot_prod(po_vector a, po_vector b);

/* the projection of vector a onto vector p */
po_vector vect_project(po_vector a, po_vector p);

/* get corresponding axis from two points; vertex1 should be the orgin */
po_vector vect_axis(po_vector vertex1, po_vector vertex2);

#endif