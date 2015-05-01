/********************************************************************
 * the header file for some handy vector math
 ********************************************************************/
#ifndef VECTOR_INCLUDED
#define VECTOR_INCLUDED

#include <math.h>

/* A generic struct for handing two dimensional vectors. Used in many
 * situations with different meanings. */
typedef struct po_vector {
  float x;
  float y;
} po_vector;

/* a 2x2 matrix */
typedef struct matrix {
  po_vector rows[2];
} matrix;

/* forms a line vector from two coordinate points */
po_vector vect_from_points(po_vector origin, po_vector destination);

/* gets the slope from two points */
float vect_slope(po_vector p1, po_vector p2);

/* gets the right hand normal vector */
po_vector vect_normal(po_vector vect);

/* turn this into a vector of length 1 */
po_vector vect_unit(po_vector vect);

/* calculate the dot product of two vectors */
float vect_dot_prod(po_vector a, po_vector b);

/* add two vectors together */
po_vector vect_add(po_vector a, po_vector b);

/* take the reciprocal of a vector  */
po_vector vect_recip(po_vector a);

/* get the maginitude of a vector squared*/
float vect_mag_squared(po_vector vect);

/* the projection of vector a onto vector p */
po_vector vect_project(po_vector a, po_vector p);

/* get corresponding axis from two points; vertex1 should be the orgin */
po_vector vect_axis(po_vector vertex1, po_vector vertex2);

/* makes 2x2 a matrix of the form 
 * {a,b}
 * {c,d} */
matrix vect_create_matrix(float a, float b, float c, float d);

/* multiplies a vector with a 2x2 matrix
 * matrix[n] for n is an int, will grab the n row of the matrix. 
 * matrix[n].x corresponds to the first column in that row, 
 * matrix[n].y to the second column in that row */
po_vector vect_matrix_mult(po_vector vert, matrix a_matrix);

/* multiplies a vector by a scalar*/
po_vector vect_scaled(po_vector vect, float a);

/* adds a scalar to a vector */
po_vector vect_add_scalar(po_vector vect, float a);


/* a function for crossing a vector with a scalar */
po_vector vect_cross_scalar(po_vector vect, float a);

/* a function for crossing a scalar with a vector
 * yes, this is different than the above function */
po_vector vect_scalar_cross(float a, po_vector vect);

/* cross product of two vectors in 2D */
float vect_cross_prod(po_vector vect1, po_vector vect2);

/* subtract second vector from first */
po_vector vect_minus(po_vector vect1,po_vector vect2);

#endif
